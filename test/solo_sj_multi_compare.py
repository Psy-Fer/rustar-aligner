#!/usr/bin/env python3
"""Diff rustar vs STARsolo SJ-feature and --soloMultiMappers matrices.

Both tools index barcodes by the same sorted whitelist (columns align directly)
and genes by the same GTF order (Gene rows align). SJ junctions differ per tool
(each has its own SJ.out.tab), so SJ rows are matched by (chr,start,end).

Reports, per matrix: shared rows/cols, total counts, Pearson r over shared
entries, and the fraction of shared entries that match exactly.

Usage:
  solo_sj_multi_compare.py --rustar <Solo.out> --starsolo <Solo.out>
"""
import argparse
import gzip
import os
import sys

import numpy as np
import scipy.io
import scipy.sparse as sp


def _open(p):
    return gzip.open(p, "rb") if p.endswith(".gz") else open(p, "rb")


def _find(d, base):
    for c in (base, base + ".gz"):
        p = os.path.join(d, c)
        if os.path.exists(p):
            return p
    return None


def load_mtx(d, name="matrix.mtx"):
    p = _find(d, name)
    if p is None:
        return None
    with _open(p) as fh:
        m = scipy.io.mmread(fh).tocsr()  # features x barcodes
    return m


def load_features_keys(d):
    """SJ features.tsv → list of (chr, start, end) per row. STARsolo symlinks
    features.tsv → SJ.out.tab (run root); fall back to that if the symlink is
    broken (it points at the in-container path)."""
    p = _find(d, "features.tsv")
    if p is None or not os.path.exists(p):
        # d = .../Solo.out/SJ/raw → run root is three levels up.
        alt = os.path.join(d, "..", "..", "..", "SJ.out.tab")
        p = alt if os.path.exists(alt) else None
    keys = []
    op = gzip.open(p, "rt") if p.endswith(".gz") else open(p)
    with op as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            keys.append((f[0], f[1], f[2]))
    return keys


def compare_aligned(name, A, B):
    """A, B are features×barcodes with identical row+col indexing."""
    if A is None or B is None:
        print(f"[{name}] missing matrix"); return
    r = min(A.shape[0], B.shape[0])
    c = min(A.shape[1], B.shape[1])
    A = A[:r, :c]
    B = B[:r, :c]
    da = np.asarray(A.sum()); db = np.asarray(B.sum())
    # union of nonzero coords
    U = (A != 0).astype(np.int8) + (B != 0).astype(np.int8)
    coo = U.tocoo()
    av = np.asarray(A[coo.row, coo.col]).ravel()
    bv = np.asarray(B[coo.row, coo.col]).ravel()
    rr = np.corrcoef(av, bv)[0, 1] if len(av) > 1 else float("nan")
    exact = np.mean(np.isclose(av, bv, atol=1e-4)) if len(av) else float("nan")
    print(f"[{name}] rustar_total={float(da):,.1f} star_total={float(db):,.1f} "
          f"shared_entries={len(av):,} r={rr:.5f} exact={exact:.4%}")


def compare_sj(rdir, sdir):
    ra = load_mtx(rdir); sa = load_mtx(sdir)
    if ra is None or sa is None:
        print("[SJ] missing matrix"); return
    rk = load_features_keys(rdir); sk = load_features_keys(sdir)
    print(f"[SJ] rustar junctions={len(rk):,} star junctions={len(sk):,}")
    sidx = {k: i for i, k in enumerate(sk)}
    shared = [(i, sidx[k]) for i, k in enumerate(rk) if k in sidx]
    print(f"[SJ] shared junctions (by chr/start/end) = {len(shared):,} "
          f"({len(shared)/max(len(rk),1):.1%} of rustar)")
    if not shared:
        return
    rrows = [i for i, _ in shared]
    srows = [j for _, j in shared]
    c = min(ra.shape[1], sa.shape[1])
    Rm = ra[rrows, :c]
    Sm = sa[srows, :c]
    U = (Rm != 0).astype(np.int8) + (Sm != 0).astype(np.int8)
    coo = U.tocoo()
    av = np.asarray(Rm[coo.row, coo.col]).ravel()
    bv = np.asarray(Sm[coo.row, coo.col]).ravel()
    rr = np.corrcoef(av, bv)[0, 1] if len(av) > 1 else float("nan")
    exact = np.mean(av == bv) if len(av) else float("nan")
    print(f"[SJ] on shared junctions: rustar_total={float(Rm.sum()):,} "
          f"star_total={float(Sm.sum()):,} shared_entries={len(av):,} "
          f"r={rr:.5f} exact={exact:.4%}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rustar", required=True, help="rustar Solo.out dir")
    ap.add_argument("--starsolo", required=True, help="STARsolo Solo.out dir")
    a = ap.parse_args()

    rg = os.path.join(a.rustar, "Gene", "raw")
    sg = os.path.join(a.starsolo, "Gene", "raw")
    print("=== Gene (unique) matrix sanity ===")
    compare_aligned("Gene", load_mtx(rg), load_mtx(sg))

    print("\n=== UniqueAndMult (--soloMultiMappers) ===")
    for method in ("Uniform", "PropUnique", "Rescue", "EM"):
        fn = f"UniqueAndMult-{method}.mtx"
        compare_aligned(method, load_mtx(rg, fn), load_mtx(sg, fn))

    print("\n=== SJ feature ===")
    compare_sj(os.path.join(a.rustar, "SJ", "raw"),
               os.path.join(a.starsolo, "SJ", "raw"))


if __name__ == "__main__":
    sys.exit(main())
