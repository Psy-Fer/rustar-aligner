#!/usr/bin/env python3
"""Wall time and peak memory for rustar-aligner against STAR, end to end.

Uses the same nf-core/rnaseq fixture as `nfcore_diff.py` (50 000 paired reads,
S. cerevisiae chrI plus the GFP transgene), timing indexing and alignment
separately because they are different questions: indexing is dominated by
suffix-array construction, alignment by the seed and stitch loops.

    python3 test/speed_bench.py --rustar ./target/release/rustar-aligner

Prints a Markdown table, so the output can go straight into a job summary. No
thresholds and no non-zero exit: a wall time from a shared CI runner is a
data point, not a gate. Compare runs on the same machine.
"""

from __future__ import annotations

import argparse
import json
import resource
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from nfcore_diff import fetch  # noqa: E402  - shared fixture download


def timed(cmd: list[str], log: Path) -> tuple[float, float]:
    """Run `cmd`, returning (wall seconds, peak RSS in MB) for the child."""
    before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    start = time.monotonic()
    with open(log, "w") as f:
        proc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    wall = time.monotonic() - start
    after = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if proc.returncode != 0:
        sys.exit(f"command failed ({proc.returncode}): {' '.join(cmd)}\nsee {log}")
    # ru_maxrss is kilobytes on Linux, bytes on macOS.
    peak = max(after, before)
    peak_mb = peak / 1024 if sys.platform != "darwin" else peak / (1024 * 1024)
    return wall, peak_mb


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rustar", default="./target/release/rustar-aligner")
    ap.add_argument("--star", default="STAR")
    ap.add_argument("--work", default="/tmp/rustar-speed")
    ap.add_argument("--threads", default="4")
    ap.add_argument("--json")
    args = ap.parse_args()

    work = Path(args.work)
    fetch(work)

    results: dict[str, dict[str, float]] = {}
    for tag, exe, run_mode in (("STAR", args.star, False), ("rustar", args.rustar, True)):
        idx = work / f"{tag}_idx"
        idx.mkdir(parents=True, exist_ok=True)
        prefix = str(work / f"{tag}_")

        gen = [exe, "--runMode", "genomeGenerate",
               "--genomeDir", str(idx),
               "--genomeFastaFiles", str(work / "genome.fasta"),
               "--genomeSAindexNbases", "9",
               "--runThreadN", args.threads,
               "--outFileNamePrefix", prefix + "idx_"]
        index_wall, index_rss = timed(gen, work / f"{tag}_index.log")

        aln = [exe]
        if run_mode:
            aln += ["--runMode", "alignReads"]
        aln += ["--genomeDir", str(idx),
                "--readFilesIn", str(work / "reads_1.fastq"), str(work / "reads_2.fastq"),
                "--outSAMtype", "SAM",
                "--runThreadN", args.threads,
                "--outFileNamePrefix", prefix]
        align_wall, align_rss = timed(aln, work / f"{tag}_align.log")

        results[tag] = {
            "index_seconds": index_wall,
            "index_peak_mb": index_rss,
            "align_seconds": align_wall,
            "align_peak_mb": align_rss,
        }

    print("### End-to-end timings\n")
    print(f"Fixture: nf-core/rnaseq test data, 50 000 pairs, {args.threads} threads.\n")
    print("| stage | STAR | rustar | ratio |")
    print("|---|---|---|---|")
    for stage, unit in (("index_seconds", "s"), ("align_seconds", "s")):
        s, r = results["STAR"][stage], results["rustar"][stage]
        ratio = r / s if s else float("inf")
        print(f"| {stage.replace('_', ' ')} | {s:.1f}{unit} | {r:.1f}{unit} | {ratio:.2f}x |")
    print()
    print("Peak RSS is measured across the whole child process group, so it is")
    print("only meaningful when the stages are run one at a time, as they are here.")

    if args.json:
        Path(args.json).write_text(json.dumps(results, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
