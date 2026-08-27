#!/usr/bin/env python3
"""Differential run against STAR on the nf-core/rnaseq test dataset.

A middle-sized fixture: 50 000 paired reads against *S. cerevisiae* chrI plus
the GFP transgene, which is small enough to fetch and index in under a minute
and large enough to move the numbers the unit tests cannot reach — mapping
rates, the unmapped-reason buckets, and the multimapper depth histogram.

Everything is fetched from public URLs; nothing is vendored.

    python3 test/nfcore_diff.py --rustar ./target/release/rustar-aligner \\
        --work /tmp/nfcore --star STAR

Exit status is 0 when every threshold holds, 1 otherwise, so it can gate a CI
job. `--report-only` always exits 0 and just prints the comparison, which is
what to use while a difference is being investigated.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import urllib.request
from pathlib import Path

REF = "626c8fab639062eade4b10747e919341cbf9b41a"
FILES = {
    "genome.fasta": f"https://raw.githubusercontent.com/nf-core/test-datasets/{REF}/reference/genome.fasta",
    "genes.gtf.gz": f"https://raw.githubusercontent.com/nf-core/test-datasets/{REF}/reference/genes_with_empty_tid.gtf.gz",
    "reads_1.fastq.gz": "https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357072_1.fastq.gz",
    "reads_2.fastq.gz": "https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357072_2.fastq.gz",
}

# How far each figure may drift from STAR before the run fails. These are not
# aspirations: they are the measured gaps plus headroom, so a regression trips
# them and today's state does not.
THRESHOLDS = {
    "uniquely_mapped_frac": 0.005,   # fraction of input reads
    "multi_mapped_frac": 0.005,
    "unmapped_short_frac": 0.010,
    "unmapped_other_frac": 0.010,
    "max_nh_ratio": 2.0,             # deepest multimapper, rustar / STAR
}


def fetch(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    for name, url in FILES.items():
        dest = work / name
        if dest.exists() and dest.stat().st_size > 0:
            continue
        print(f"fetching {name}", flush=True)
        urllib.request.urlretrieve(url, dest)  # noqa: S310 - fixed public URLs

    # Decompress the reads. `--readFilesCommand zcat` is deliberately not used:
    # on macOS zcat silently yields nothing for both aligners, and a fixture
    # that reads as "0 input reads" is worse than a slower one.
    import gzip
    import shutil

    for n in (1, 2):
        plain = work / f"reads_{n}.fastq"
        if not plain.exists():
            with gzip.open(work / f"reads_{n}.fastq.gz", "rb") as fin, open(plain, "wb") as fout:
                shutil.copyfileobj(fin, fout)


def run(cmd: list[str], log: Path) -> None:
    with open(log, "w") as f:
        proc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        sys.exit(f"command failed ({proc.returncode}): {' '.join(cmd)}\nsee {log}")


def align(exe: str, work: Path, tag: str, is_star: bool) -> Path:
    idx = work / f"{tag}_idx"
    idx.mkdir(exist_ok=True)
    prefix = str(work / f"{tag}_")
    gen = [exe]
    if not is_star:
        gen += ["--runMode", "genomeGenerate"]
    else:
        gen += ["--runMode", "genomeGenerate"]
    gen += [
        "--genomeDir", str(idx),
        "--genomeFastaFiles", str(work / "genome.fasta"),
        "--genomeSAindexNbases", "9",
        "--outFileNamePrefix", prefix + "idx_",
    ]
    run(gen, work / f"{tag}_index.log")

    aln = [exe]
    if not is_star:
        aln += ["--runMode", "alignReads"]
    aln += [
        "--genomeDir", str(idx),
        "--readFilesIn", str(work / "reads_1.fastq"), str(work / "reads_2.fastq"),
        "--outFilterMultimapNmax", "20",
        "--outSAMtype", "SAM",
        "--runThreadN", "4",
        "--outFileNamePrefix", prefix,
    ]
    run(aln, work / f"{tag}_align.log")
    return Path(prefix)


def final_log(prefix: Path) -> dict[str, int]:
    rows: dict[str, int] = {}
    for line in open(f"{prefix}Log.final.out"):
        if "|" not in line:
            continue
        label, value = line.split("|", 1)
        value = value.strip()
        if value.endswith("%"):
            continue
        try:
            rows[label.strip()] = int(value)
        except ValueError:
            pass
    return rows


def nh_histogram(prefix: Path) -> dict[int, int]:
    hist: dict[int, int] = {}
    with open(f"{prefix}Aligned.out.sam") as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            flag = int(fields[1])
            if flag & 0x100 or flag & 0x800 or flag & 0x4:
                continue
            for tag in fields[11:]:
                if tag.startswith("NH:i:"):
                    nh = int(tag[5:])
                    hist[nh] = hist.get(nh, 0) + 1
                    break
    return hist


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rustar", default="./target/release/rustar-aligner")
    ap.add_argument("--star", default="STAR")
    ap.add_argument("--work", default="/tmp/nfcore-diff")
    ap.add_argument("--report-only", action="store_true")
    ap.add_argument("--json", help="write the measured figures here")
    args = ap.parse_args()

    work = Path(args.work)
    fetch(work)

    star_prefix = align(args.star, work, "star", is_star=True)
    rustar_prefix = align(args.rustar, work, "rustar", is_star=False)

    s, r = final_log(star_prefix), final_log(rustar_prefix)
    s_nh, r_nh = nh_histogram(star_prefix), nh_histogram(rustar_prefix)
    n_input = s["Number of input reads"]

    fields = [
        ("uniquely_mapped_frac", "Uniquely mapped reads number"),
        ("multi_mapped_frac", "Number of reads mapped to multiple loci"),
        ("unmapped_short_frac", "Number of reads unmapped: too short"),
        ("unmapped_other_frac", "Number of reads unmapped: other"),
    ]

    measured: dict[str, float] = {}
    failures: list[str] = []

    print(f"\n{'metric':<40} {'STAR':>10} {'rustar':>10} {'delta':>10}")
    for key, label in fields:
        sv, rv = s.get(label, 0), r.get(label, 0)
        delta = abs(sv - rv) / n_input
        measured[key] = delta
        print(f"{label:<40} {sv:>10} {rv:>10} {delta:>9.3%}")
        if delta > THRESHOLDS[key]:
            failures.append(f"{label}: |{sv} - {rv}| / {n_input} = {delta:.3%} > {THRESHOLDS[key]:.3%}")

    s_max, r_max = max(s_nh, default=1), max(r_nh, default=1)
    ratio = r_max / s_max if s_max else float("inf")
    measured["max_nh_ratio"] = ratio
    print(f"{'deepest multimapper (NH)':<40} {s_max:>10} {r_max:>10} {ratio:>9.2f}x")
    if ratio > THRESHOLDS["max_nh_ratio"]:
        failures.append(f"deepest NH: {r_max} against STAR's {s_max} ({ratio:.2f}x > {THRESHOLDS['max_nh_ratio']}x)")

    if args.json:
        Path(args.json).write_text(json.dumps({"measured": measured, "star": s, "rustar": r}, indent=2))

    if failures:
        print("\nFAILED:")
        for f in failures:
            print(f"  {f}")
        return 0 if args.report_only else 1

    print("\nNFCORE DIFF PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
