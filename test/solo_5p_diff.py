#!/usr/bin/env python3
"""Differential test for the 5' paired-end solo path (`--soloBarcodeMate 1`).

Synthesizes a small 5' 10x-style dataset where mate 1 carries the barcode+UMI
followed by cDNA (mate 1 = CB+UMI+spacer+cDNA), and mate 2 is the paired cDNA
read, then runs BOTH real STAR and rustar-aligner with the cellgeni 5' flag set
(`--soloBarcodeMate 1 --clip5pNbases 39 0` + the CellRanger-matching UMI/CB flags)
and compares the raw Gene and GeneFull count matrices.

Usage:
    python3 test/solo_5p_diff.py --star /path/to/STAR --rustar /path/to/rustar-aligner

Exit 0 = matrices match, 1 = mismatch / error.
"""
import argparse
import os
import random
import shutil
import subprocess
import sys
import tempfile

# CellRanger-matching flags used by cellgeni's 5' PE command (no CellRanger4 adapter
# clip — the 5' path uses clip5pNbases instead).
CR_FLAGS = [
    "--outFilterScoreMin", "30",
    "--soloCBmatchWLtype", "1MM_multi_Nbase_pseudocounts",
    "--soloUMIfiltering", "MultiGeneUMI_CR",
    "--soloUMIdedup", "1MM_CR",
]
# Barcode + spacer clipped off mate 1's 5' before alignment (cellgeni: 39).
CLIP5P_M1 = 39
CB_LEN = 16
UMI_LEN = 12
SPACER = 11               # CB(16)+UMI(12)+spacer(11) = 39 = CLIP5P_M1
CDNA_LEN = 80
BASES = "ACGT"
_COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

GENE_A_START = 10000
GENE_B_START = 30000


def rand_seq(rng, n):
    return "".join(rng.choice(BASES) for _ in range(n))


def revcomp(s):
    return "".join(_COMP[b] for b in reversed(s))


def _plant_gene(g, s, rng):
    g[s : s + 150] = list(rand_seq(rng, 150))          # exon1
    g[s + 150 : s + 400] = list(rand_seq(rng, 250))    # intron
    g[s + 150], g[s + 151] = "G", "T"                  # donor
    g[s + 398], g[s + 399] = "A", "G"                  # acceptor
    g[s + 400 : s + 550] = list(rand_seq(rng, 150))    # exon2


def build_genome(rng, length=50000):
    g = list(rand_seq(rng, length))
    _plant_gene(g, GENE_A_START, rng)
    _plant_gene(g, GENE_B_START, rng)
    return "".join(g)


def gene_pair(genome, exon_start):
    """A concordant FR cDNA pair inside exon1: mate1 forward at [a, a+CDNA_LEN),
    mate2 = revcomp of [a+20, a+20+CDNA_LEN) (both inside the 150 bp exon)."""
    a = exon_start + 20
    m1 = genome[a : a + CDNA_LEN]
    m2 = revcomp(genome[a + 20 : a + 20 + CDNA_LEN])
    return m1, m2


def write_files(d, genome, rng):
    fa = os.path.join(d, "genome.fa")
    with open(fa, "w") as f:
        f.write(">chr1\n")
        for i in range(0, len(genome), 70):
            f.write(genome[i : i + 70] + "\n")
    gtf = os.path.join(d, "genes.gtf")
    with open(gtf, "w") as f:
        f.write('chr1\tsrc\texon\t10001\t10150\t.\t+\t.\tgene_id "GENEA"; transcript_id "GENEA.1"; gene_name "GeneA";\n')
        f.write('chr1\tsrc\texon\t10401\t10550\t.\t+\t.\tgene_id "GENEA"; transcript_id "GENEA.1"; gene_name "GeneA";\n')
        f.write('chr1\tsrc\texon\t30001\t30150\t.\t+\t.\tgene_id "GENEB"; transcript_id "GENEB.1"; gene_name "GeneB";\n')
        f.write('chr1\tsrc\texon\t30401\t30550\t.\t+\t.\tgene_id "GENEB"; transcript_id "GENEB.1"; gene_name "GeneB";\n')
    wl = os.path.join(d, "whitelist.txt")
    cbs = ["AAAACCCCGGGGTTTT", "ACACACACGTGTGTGT", "TTTTGGGGCCCCAAAA", "GTGTGTGTACACACAC"]
    with open(wl, "w") as f:
        f.write("\n".join(cbs) + "\n")

    a1, a2 = gene_pair(genome, GENE_A_START)
    b1, b2 = gene_pair(genome, GENE_B_START)
    # (cell, (mate1_cdna, mate2_cdna), umi, n_reads)
    plan = [
        (cbs[0], (a1, a2), "ACGTACGTACGT", 5),
        (cbs[0], (a1, a2), "ACGTACGTACGA", 1),   # 1MM_CR neighbor -> collapses
        (cbs[0], (a1, a2), "TGCATGCATGCA", 3),
        (cbs[0], (b1, b2), "GGGGTTTTAACC", 2),
        (cbs[1], (a1, a2), "CATGCATGCATG", 4),
    ]
    r1 = os.path.join(d, "R1.fq")  # mate 1: CB+UMI+spacer+cDNA
    r2 = os.path.join(d, "R2.fq")  # mate 2: cDNA
    ci = 0
    with open(r1, "w") as f1, open(r2, "w") as f2:
        for (cb, (m1c, m2c), umi, n) in plan:
            for _ in range(n):
                name = f"read{ci}"; ci += 1
                spacer = rand_seq(rng, SPACER)
                m1 = cb + umi + spacer + m1c
                f1.write(f"@{name}\n{m1}\n+\n{'I' * len(m1)}\n")
                f2.write(f"@{name}\n{m2c}\n+\n{'I' * len(m2c)}\n")
    return fa, gtf, wl, r1, r2


def run(cmd, **kw):
    print("  $", " ".join(str(c) for c in cmd))
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        print(r.stdout[-2000:]); print(r.stderr[-4000:])
        raise SystemExit(f"command failed ({r.returncode}): {cmd[0]}")
    return r


def _solo_geometry(wl):
    return [
        "--soloType", "CB_UMI_Simple", "--soloCBwhitelist", wl,
        "--soloCBstart", "1", "--soloCBlen", str(CB_LEN),
        "--soloUMIstart", str(CB_LEN + 1), "--soloUMIlen", str(UMI_LEN),
        "--soloBarcodeMate", "1", "--clip5pNbases", str(CLIP5P_M1), "0",
        "--soloStrand", "Forward",
        "--soloFeatures", "Gene", "GeneFull",
    ]


def run_star(star, d, fa, gtf, wl, r1, r2):
    idx = os.path.join(d, "star_index"); os.makedirs(idx, exist_ok=True)
    run([star, "--runMode", "genomeGenerate", "--genomeDir", idx,
         "--genomeFastaFiles", fa, "--sjdbGTFfile", gtf,
         "--genomeSAindexNbases", "7", "--sjdbOverhang", "89"])
    gp = os.path.join(idx, "genomeParameters.txt")
    lines = open(gp).read().splitlines()
    with open(gp, "w") as f:
        for ln in lines:
            f.write("sjdbGTFfile\t-\n" if ln.startswith("sjdbGTFfile\t") else ln + "\n")
    out = os.path.join(d, "star_out") + os.sep
    run([star, "--genomeDir", idx, "--readFilesIn", r1, r2]
        + _solo_geometry(wl) + ["--outSAMtype", "None", "--outFileNamePrefix", out] + CR_FLAGS)
    log = os.path.join(out, "Log.final.out")
    if os.path.exists(log):
        for ln in open(log):
            if "Number of input reads" in ln and ln.strip().endswith("0"):
                raise SystemExit("STAR read 0 reads (broken binary on this host)")
    return os.path.join(out, "Solo.out")


def run_rustar(rustar, d, fa, gtf, wl, r1, r2):
    idx = os.path.join(d, "rustar_index"); os.makedirs(idx, exist_ok=True)
    run([rustar, "--runMode", "genomeGenerate", "--genomeDir", idx,
         "--genomeFastaFiles", fa, "--sjdbGTFfile", gtf,
         "--genomeSAindexNbases", "7", "--sjdbOverhang", "89"])
    out = os.path.join(d, "rustar_out") + os.sep
    run([rustar, "--genomeDir", idx, "--readFilesIn", r1, r2, "--sjdbGTFfile", gtf]
        + _solo_geometry(wl) + ["--outSAMtype", "None", "--outFileNamePrefix", out] + CR_FLAGS)
    return os.path.join(out, "Solo.out")


def decode_matrix(raw_dir):
    feats, barcodes = [], []
    with open(os.path.join(raw_dir, "features.tsv")) as f:
        feats = [l.rstrip("\n").split("\t")[0] for l in f]
    with open(os.path.join(raw_dir, "barcodes.tsv")) as f:
        barcodes = [l.strip() for l in f]
    out = {}
    with open(os.path.join(raw_dir, "matrix.mtx")) as f:
        lines = [l for l in f if not l.startswith("%")]
    for entry in lines[1:]:
        p = entry.split()
        if len(p) >= 3:
            out[(barcodes[int(p[1]) - 1], feats[int(p[0]) - 1])] = int(float(p[2]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--star", default=shutil.which("STAR"))
    ap.add_argument("--rustar", required=True)
    ap.add_argument("--keep", action="store_true")
    ap.add_argument("--seed", type=int, default=20260710)
    args = ap.parse_args()
    if not (args.star and os.path.exists(args.star)):
        raise SystemExit(f"STAR not found: {args.star}")
    if not os.path.exists(args.rustar):
        raise SystemExit(f"rustar not found: {args.rustar}")

    d = tempfile.mkdtemp(prefix="solo5p_")
    print(f"workdir: {d}\nSTAR: {args.star}\nrustar: {args.rustar}")
    rng = random.Random(args.seed)
    try:
        genome = build_genome(rng)
        fa, gtf, wl, r1, r2 = write_files(d, genome, rng)

        print("\n== rustar (5' PE) ==")
        rustar_solo = run_rustar(args.rustar, d, fa, gtf, wl, r1, r2)
        print("\n== STAR (5' PE) ==")
        star_solo = run_star(args.star, d, fa, gtf, wl, r1, r2)

        ok = True
        for feat in ("Gene", "GeneFull"):
            rm = decode_matrix(os.path.join(rustar_solo, feat, "raw"))
            sm = decode_matrix(os.path.join(star_solo, feat, "raw"))
            print(f"\n== {feat}: STAR={len(sm)} entries, rustar={len(rm)} entries ==")
            if sm == rm:
                print(f"  PASS: {feat} raw matrix byte-identical to STARsolo ({len(sm)} entries)")
            else:
                ok = False
                print(f"  FAIL: {feat} mismatch:")
                for k in sorted(set(sm) | set(rm)):
                    if sm.get(k) != rm.get(k):
                        print(f"     {k}: STAR={sm.get(k)} rustar={rm.get(k)}")
        return 0 if ok else 1
    finally:
        if args.keep:
            print(f"(kept {d})")
        else:
            shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
