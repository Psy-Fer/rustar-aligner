#!/usr/bin/env python3
"""Compare GeneFull (intron-inclusive) quantification across rustar / STARsolo /
CellRanger, plus the EmptyDrops-filtered cell sets.

Part A (raw count parity): load each tool's raw matrix, report total UMIs, genes
detected, cells with >0 UMI, and the per-cell UMI-total correlation between
rustar-GeneFull and STARsolo-GeneFull (they should match closely) and each vs
CellRanger (whose default raw matrix is intron-inclusive).

Part B (filtered h5 parity): given EmptyDrops-called barcode lists for each tool
(from the `emptydrops` Rust binary) and CellRanger's own filtered barcodes,
report cell-set overlap (Jaccard) and per-cell UMI agreement on shared cells.

Usage:
  solo_genefull_compare.py \
    --rustar  <rustar genefull raw dir> \
    --starsolo <starsolo genefull raw dir> \
    --cellranger <cellranger raw dir> \
    [--rustar-cells f.txt --starsolo-cells f.txt --cr-cells f.txt] \
    --out compare_genefull.json
"""
import argparse
import gzip
import json
import os
import sys


def _open(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def _find(d, base):
    for c in (base, base + ".gz"):
        p = os.path.join(d, c)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"{base}[.gz] not in {d}")


def load_raw(d):
    """Return (barcodes list, dict cell_idx->total_umi, n_genes, total_umi).

    Barcodes are normalized (10x '-1' gem-group suffix stripped) so they are
    comparable across tools (CellRanger appends '-1', STARsolo/rustar do not)."""
    bcs = [l.split("\t")[0].split("-")[0].strip() for l in _open(_find(d, "barcodes.tsv"))]
    genes = [l.split("\t")[0].strip() for l in _open(_find(d, "features.tsv"))]
    totals = {}
    total_umi = 0
    genes_seen = set()
    with _open(_find(d, "matrix.mtx")) as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            break  # first non-% line is the dims header; skip it
        for line in fh:
            g, c, v = line.split()[:3]
            v = int(float(v))
            if v == 0:
                continue
            ci = int(c) - 1
            totals[ci] = totals.get(ci, 0) + v
            total_umi += v
            genes_seen.add(int(g) - 1)
    return bcs, totals, len(genes), total_umi, len(genes_seen)


def summarize(name, d):
    bcs, totals, n_genes, total_umi, genes_detected = load_raw(d)
    cells = sum(1 for v in totals.values() if v > 0)
    print(f"[{name}] cells>0={cells:,}  total_UMI={total_umi:,}  "
          f"genes_detected={genes_detected:,}/{n_genes:,}")
    return {"name": name, "barcodes": bcs, "totals": totals,
            "cells_gt0": cells, "total_umi": total_umi,
            "genes_detected": genes_detected, "n_genes": n_genes}


def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx == 0 or syy == 0:
        return float("nan")
    return sxy / (sxx ** 0.5 * syy ** 0.5)


def per_cell_corr(a, b):
    """Per-cell UMI-total correlation over the shared barcode set."""
    a_by_bc = {a["barcodes"][i]: t for i, t in a["totals"].items()}
    b_by_bc = {b["barcodes"][i]: t for i, t in b["totals"].items()}
    shared = sorted(set(a_by_bc) & set(b_by_bc))
    xs = [a_by_bc[bc] for bc in shared]
    ys = [b_by_bc[bc] for bc in shared]
    r = pearson(xs, ys)
    exact = sum(1 for x, y in zip(xs, ys) if x == y)
    return {"shared_cells": len(shared), "pearson_r": r,
            "exact_total_match": exact,
            "exact_frac": exact / len(shared) if shared else float("nan")}


def read_cells(p):
    if not p or not os.path.exists(p):
        return None
    return set(l.split("\t")[0].strip() for l in _open(p))


def jaccard(a, b):
    if a is None or b is None:
        return None
    inter = len(a & b)
    union = len(a | b)
    return {"a": len(a), "b": len(b), "intersection": inter,
            "jaccard": inter / union if union else float("nan"),
            "a_only": len(a - b), "b_only": len(b - a)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rustar", required=True)
    ap.add_argument("--starsolo", required=True)
    ap.add_argument("--cellranger", required=True)
    ap.add_argument("--rustar-cells")
    ap.add_argument("--starsolo-cells")
    ap.add_argument("--cr-cells")
    ap.add_argument("--out", default="compare_genefull.json")
    a = ap.parse_args()

    print("=== Part A: GeneFull raw count parity ===")
    R = summarize("rustar-GeneFull", a.rustar)
    S = summarize("STARsolo-GeneFull", a.starsolo)
    C = summarize("CellRanger-raw", a.cellranger)

    print("\n=== per-cell UMI-total correlation ===")
    rs = per_cell_corr(R, S)
    print(f"rustar vs STARsolo : shared={rs['shared_cells']:,}  r={rs['pearson_r']:.6f}  "
          f"exact_total={rs['exact_frac']:.4%}")
    rc = per_cell_corr(R, C)
    print(f"rustar vs CellRgr  : shared={rc['shared_cells']:,}  r={rc['pearson_r']:.6f}")
    sc = per_cell_corr(S, C)
    print(f"STAR   vs CellRgr  : shared={sc['shared_cells']:,}  r={sc['pearson_r']:.6f}")

    out = {
        "raw": {k: {kk: v[kk] for kk in ("cells_gt0", "total_umi",
                                          "genes_detected", "n_genes")}
                for k, v in (("rustar", R), ("starsolo", S), ("cellranger", C))},
        "corr": {"rustar_vs_starsolo": rs, "rustar_vs_cr": rc, "starsolo_vs_cr": sc},
    }

    rcells = read_cells(a.rustar_cells)
    scells = read_cells(a.starsolo_cells)
    ccells = read_cells(a.cr_cells)
    if rcells or ccells:
        print("\n=== Part B: filtered cell-set overlap (EmptyDrops / CellRanger) ===")
        out["filtered"] = {}
        if rcells and ccells:
            j = jaccard(rcells, ccells)
            print(f"rustar-ED vs CR-filtered : rustar={j['a']:,} CR={j['b']:,} "
                  f"shared={j['intersection']:,} jaccard={j['jaccard']:.4f}")
            out["filtered"]["rustar_vs_cr"] = j
        if scells and ccells:
            j = jaccard(scells, ccells)
            print(f"STAR-ED   vs CR-filtered : star={j['a']:,} CR={j['b']:,} "
                  f"shared={j['intersection']:,} jaccard={j['jaccard']:.4f}")
            out["filtered"]["starsolo_vs_cr"] = j
        if rcells and scells:
            j = jaccard(rcells, scells)
            print(f"rustar-ED vs STAR-ED     : jaccard={j['jaccard']:.4f}")
            out["filtered"]["rustar_vs_starsolo"] = j

    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=2, default=str)
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    sys.exit(main())
