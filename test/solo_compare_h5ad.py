#!/usr/bin/env python3
"""Knee-call + compare CellRanger / STARsolo / rustar-aligner raw matrices.

For a fair comparison that isolates *counting* differences from *cell-calling*
differences, the SAME knee filter (CellRanger 2.2 — STARsolo's default
--soloCellFilter) is applied to each tool's RAW matrix. Each filtered result is
written as an .h5ad (AnnData, cells x genes) and the three are compared:
n cells, median UMI/genes per cell, barcode overlap, per-cell UMI correlation on
shared barcodes, and gene-level pseudobulk correlation.

Usage:
    .venv/bin/python test/solo_compare_h5ad.py \
        --cellranger <outs/raw_feature_bc_matrix> \
        --starsolo   <Solo.out/Gene/raw> \
        --rustar     <Solo.out/Gene/raw> \
        --out <dir>
"""
import argparse
import gzip
import json
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


def _find(d, base):
    for c in (base, base + ".gz"):
        p = os.path.join(d, c)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"{base}[.gz] not found in {d}")


def _open_text(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def load_raw(d):
    """Load a 10x/STARsolo raw matrix dir -> (X cells x genes CSR, barcodes, gene_ids)."""
    mp = _find(d, "matrix.mtx")
    handle = gzip.open(mp, "rb") if mp.endswith(".gz") else open(mp, "rb")
    with handle:
        m = scipy.io.mmread(handle)  # features x barcodes
    X = sp.csr_matrix(m).T.tocsr()  # -> barcodes (cells) x features (genes)
    barcodes = np.array([l.split("\t")[0].strip() for l in _open_text(_find(d, "barcodes.tsv"))])
    genes = np.array([l.split("\t")[0].strip() for l in _open_text(_find(d, "features.tsv"))])
    return X, barcodes, genes


def norm_bc(bc):
    """Strip 10x '-1' gem-group suffix so barcodes are comparable across tools."""
    return np.array([b.split("-")[0] for b in bc])


def revcomp(s):
    t = str.maketrans("ACGT", "TGCA")
    return s.translate(t)[::-1]


def knee_cr22(totals, n_expected=3000, max_pct=0.99, max_min_ratio=10):
    """CellRanger-2.2 knee threshold on per-barcode totals (STARsolo default)."""
    counts = np.sort(totals[totals > 0])[::-1]
    if counts.size == 0:
        return 0.0
    idx = min(int(round(n_expected * (1 - max_pct))), counts.size - 1)
    robust_max = counts[idx]
    return robust_max / max_min_ratio


def load_cell_set(path):
    """Load an EmptyDrops cells.txt (one barcode/line) -> normalized set, or None."""
    if not path or not os.path.exists(path):
        return None
    with _open_text(path) as fh:
        return set(l.split("\t")[0].split("-")[0].strip() for l in fh if l.strip())


def build_filtered(name, raw_dir, rc_barcodes=False, cells=None):
    """Filter a raw matrix to called cells. If `cells` (a normalized barcode set,
    e.g. from EmptyDrops) is given, keep exactly those; otherwise CR2.2 knee."""
    X, bc, genes = load_raw(raw_dir)
    bc = norm_bc(bc)
    if rc_barcodes:
        bc = np.array([revcomp(b) for b in bc])
    totals = np.asarray(X.sum(axis=1)).ravel()
    if cells is not None:
        thr = -1.0
        keep = np.array([b in cells for b in bc])
    else:
        thr = knee_cr22(totals)
        keep = totals >= thr
    Xf = X[keep]
    bcf = bc[keep]
    A = ad.AnnData(X=Xf, obs=pd.DataFrame(index=bcf), var=pd.DataFrame(index=genes))
    A.obs["n_umi"] = np.asarray(Xf.sum(axis=1)).ravel()
    A.obs["n_genes"] = np.asarray((Xf > 0).sum(axis=1)).ravel()
    return A, thr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cellranger", required=True)
    ap.add_argument("--starsolo", required=True)
    ap.add_argument("--rustar", required=True)
    ap.add_argument("--out", required=True)
    # Optional EmptyDrops cells.txt per tool; when given, filter by these calls
    # instead of the CR2.2 knee (CellRanger uses its own filtered barcodes).
    ap.add_argument("--rustar-cells")
    ap.add_argument("--starsolo-cells")
    ap.add_argument("--cellranger-cells")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    r_cells = load_cell_set(args.rustar_cells)
    s_cells = load_cell_set(args.starsolo_cells)
    c_cells = load_cell_set(args.cellranger_cells)

    # Build STARsolo / rustar first; detect whether CellRanger barcodes need RC
    # (some 5' chemistries report the reverse complement).
    star, star_thr = build_filtered("STARsolo", args.starsolo, cells=s_cells)
    rust, rust_thr = build_filtered("rustar", args.rustar, cells=r_cells)

    cr_plain, _ = build_filtered("CellRanger", args.cellranger, rc_barcodes=False, cells=c_cells)
    ov_plain = len(set(cr_plain.obs_names) & set(star.obs_names))
    cr_rc, _ = build_filtered("CellRanger", args.cellranger, rc_barcodes=True, cells=c_cells)
    ov_rc = len(set(cr_rc.obs_names) & set(star.obs_names))
    cr = cr_rc if ov_rc > ov_plain else cr_plain
    cr_orient = "reverse-complement" if ov_rc > ov_plain else "as-reported"

    objs = {"CellRanger": cr, "STARsolo": star, "rustar-aligner": rust}
    for name, A in objs.items():
        path = os.path.join(args.out, f"{name.replace('-aligner','')}.filtered.h5ad")
        A.write_h5ad(path)

    print(f"\nCellRanger barcode orientation vs STARsolo: {cr_orient} "
          f"(overlap as-reported={ov_plain}, rc={ov_rc})")

    # ---- per-tool summary ----
    print("\n================ filtered (CR2.2 knee) summary ================")
    hdr = f"{'tool':<16}{'cells':>8}{'median UMI/cell':>17}{'median genes/cell':>19}{'genes detected':>16}{'total UMI':>12}"
    print(hdr); print("-" * len(hdr))
    rows = {}
    for name, A in objs.items():
        med_umi = int(np.median(A.obs["n_umi"])) if A.n_obs else 0
        med_g = int(np.median(A.obs["n_genes"])) if A.n_obs else 0
        genes_det = int((np.asarray(A.X.sum(axis=0)).ravel() > 0).sum())
        tot = int(A.X.sum())
        rows[name] = dict(cells=A.n_obs, median_umi=med_umi, median_genes=med_g,
                          genes_detected=genes_det, total_umi=tot)
        print(f"{name:<16}{A.n_obs:>8}{med_umi:>17}{med_g:>19}{genes_det:>16}{tot:>12}")

    # ---- barcode overlap (called-cell sets) ----
    sets = {n: set(A.obs_names) for n, A in objs.items()}
    names = list(objs)
    print("\n================ called-cell barcode overlap ================")
    allc = sets[names[0]] & sets[names[1]] & sets[names[2]]
    print(f"shared by all 3: {len(allc)}")
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            inter = len(sets[a] & sets[b]); uni = len(sets[a] | sets[b])
            print(f"  {a} ∩ {b}: {inter}  (Jaccard {inter/uni:.3f})")

    # ---- correlations on shared cells & genes ----
    print("\n================ agreement on shared cells/genes ================")
    shared_genes = list(set(cr.var_names) & set(star.var_names) & set(rust.var_names))
    common_cells = sorted(allc)
    corr = {}
    if common_cells and shared_genes:
        # per-cell total UMI vectors (aligned to common cells)
        def cell_totals(A):
            idx = [A.obs_names.get_loc(c) for c in common_cells]
            return np.asarray(A[idx].X.sum(axis=1)).ravel()
        tot = {n: cell_totals(A) for n, A in objs.items()}
        # pseudobulk per gene (sum over shared cells), aligned to shared genes
        def pseudobulk(A):
            idx = [A.obs_names.get_loc(c) for c in common_cells]
            gi = [A.var_names.get_loc(g) for g in shared_genes]
            return np.asarray(A[idx][:, gi].X.sum(axis=0)).ravel()
        pb = {n: pseudobulk(A) for n, A in objs.items()}
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                a, b = names[i], names[j]
                rc_cell = np.corrcoef(tot[a], tot[b])[0, 1]
                rc_gene = np.corrcoef(pb[a], pb[b])[0, 1]
                corr[f"{a} vs {b}"] = dict(per_cell_umi_r=round(float(rc_cell), 4),
                                           pseudobulk_gene_r=round(float(rc_gene), 4))
                print(f"  {a} vs {b}: per-cell UMI r={rc_cell:.4f}, gene pseudobulk r={rc_gene:.4f}  "
                      f"(n_cells={len(common_cells)}, n_genes={len(shared_genes)})")

    out = dict(threshold=dict(STARsolo=star_thr, rustar=rust_thr),
               cellranger_orientation=cr_orient, summary=rows, correlations=corr,
               shared_all3_cells=len(allc))
    with open(os.path.join(args.out, "compare.json"), "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote {len(objs)} h5ad files + compare.json to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
