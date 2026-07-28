#!/usr/bin/env python3
"""Cross-compare CellRanger-style summary metrics across rustar / STARsolo /
CellRanger.

rustar and STARsolo emit `Solo.out/<feature>/Summary.csv` (key,value with
fractions in [0,1]); CellRanger emits `metrics_summary.csv` (one header row + one
value row, percentages like "53.5%" and comma-grouped integers). This pulls the
shared metrics into one table — genome/exon/intron/intergenic mapping rates plus
per-cell UMI/gene stats.

Usage:
  solo_summary_compare.py \
    --rustar    <rustar Solo.out/GeneFull/Summary.csv> \
    --starsolo  <STARsolo Solo.out/GeneFull/Summary.csv> \
    --cellranger <outs/metrics_summary.csv>
"""
import argparse
import csv
import sys


def load_summary_csv(path):
    """rustar/STARsolo Summary.csv -> {key: float-or-int}."""
    d = {}
    with open(path) as fh:
        for line in fh:
            if "," not in line:
                continue
            k, v = line.rstrip("\n").split(",", 1)
            try:
                d[k] = float(v) if "." in v else int(v)
            except ValueError:
                d[k] = v
    return d


def load_cr_metrics(path):
    """CellRanger metrics_summary.csv -> {key: float} (percents -> fraction)."""
    with open(path) as fh:
        rows = list(csv.reader(fh))
    keys, vals = rows[0], rows[1]
    out = {}
    for k, v in zip(keys, vals):
        v = v.strip()
        if v.endswith("%"):
            out[k] = float(v[:-1]) / 100.0
        else:
            try:
                out[k] = float(v.replace(",", ""))
            except ValueError:
                out[k] = v
    return out


def fmt_pct(x):
    return f"{x*100:.1f}%" if isinstance(x, (int, float)) else str(x)


def fmt_int(x):
    return f"{int(x):,}" if isinstance(x, (int, float)) else str(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rustar", required=True, help="rustar GeneFull Summary.csv")
    ap.add_argument("--starsolo", required=True, help="STARsolo GeneFull Summary.csv")
    ap.add_argument("--cellranger", required=True, help="CellRanger metrics_summary.csv")
    ap.add_argument("--feature", default="GeneFull")
    a = ap.parse_args()

    R = load_summary_csv(a.rustar)
    S = load_summary_csv(a.starsolo)
    C = load_cr_metrics(a.cellranger)
    f = a.feature

    # (label, rustar key, starsolo key, cellranger key, formatter)
    pct = fmt_pct
    intg = fmt_int
    rows = [
        ("Valid barcodes", "Reads With Valid Barcodes", "Reads With Valid Barcodes", "Valid Barcodes", pct),
        ("Sequencing saturation", "Sequencing Saturation", "Sequencing Saturation", "Sequencing Saturation", pct),
        ("Reads mapped to genome (U+M)", "Reads Mapped to Genome: Unique+Multiple", "Reads Mapped to Genome: Unique+Multiple", "Reads Mapped to Genome", pct),
        ("  ... exonic", "Reads Mapped Confidently to Exonic Regions", None, "Reads Mapped Confidently to Exonic Regions", pct),
        ("  ... intronic", "Reads Mapped Confidently to Intronic Regions", None, "Reads Mapped Confidently to Intronic Regions", pct),
        ("  ... intergenic", "Reads Mapped Confidently to Intergenic Regions", None, "Reads Mapped Confidently to Intergenic Regions", pct),
        ("Reads antisense to gene", "Reads Mapped Antisense to Gene", None, "Reads Mapped Antisense to Gene", pct),
        ("Estimated number of cells", "Estimated Number of Cells", "Estimated Number of Cells", "Estimated Number of Cells", intg),
        ("Mean reads / cell", "Mean Reads per Cell", "Mean Reads per Cell", "Mean Reads per Cell", intg),
        (f"Median genes / cell", f"Median {f} per Cell", f"Median {f} per Cell", "Median Genes per Cell", intg),
        ("Median UMI / cell", "Median UMI per Cell", "Median UMI per Cell", "Median UMI Counts per Cell", intg),
        ("Total genes detected", f"Total {f} Detected", f"Total {f} Detected", "Total Genes Detected", intg),
        ("Fraction reads in cells", "Fraction of Unique Reads in Cells", "Fraction of Unique Reads in Cells", "Fraction Reads in Cells", pct),
    ]

    w = 34
    print(f"\nCross-tool summary ({f} for rustar/STARsolo; CellRanger raw is intron-inclusive)\n")
    print(f"{'metric':<{w}}{'rustar':>14}{'STARsolo':>14}{'CellRanger':>14}")
    print("-" * (w + 42))
    for label, rk, sk, ck, fn in rows:
        rv = fn(R.get(rk)) if rk and rk in R else "—"
        sv = fn(S.get(sk)) if sk and sk in S else "—"
        cv = fn(C.get(ck)) if ck and ck in C else "—"
        print(f"{label:<{w}}{rv:>14}{sv:>14}{cv:>14}")
    print()


if __name__ == "__main__":
    sys.exit(main())
