#!/usr/bin/env python3
"""GeneFull intron-gap + EmptyDrops-filtered h5 comparison (rustar vs CellRanger).

Loads matrices one at a time (memory-careful), reports:
  A. intron effect — rustar Gene vs GeneFull total UMI (same cells);
  B. raw-count parity — rustar GeneFull vs CellRanger raw total UMI / genes;
  C. cell-set agreement — rustar EmptyDrops cells vs CellRanger native filtered,
     and rustar-ED vs CellRanger-raw+same-EmptyDrops (isolates algorithm);
  D. per-cell UMI correlation on the shared filtered cells;
  writes rustar.GeneFull.filtered.h5ad + CellRanger.filtered.h5ad.
"""
import argparse, gzip, json, os, sys
import numpy as np, scipy.io, scipy.sparse as sp, anndata as ad, pandas as pd


def _find(d, base):
    for c in (base, base + ".gz"):
        p = os.path.join(d, c)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"{base}[.gz] in {d}")


def _open(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def load(d):
    mp = _find(d, "matrix.mtx")
    h = gzip.open(mp, "rb") if mp.endswith(".gz") else open(mp, "rb")
    with h:
        X = sp.csr_matrix(scipy.io.mmread(h)).T.tocsr()  # cells x genes
    bc = np.array([l.split("\t")[0].split("-")[0].strip() for l in _open(_find(d, "barcodes.tsv"))])
    genes = np.array([l.split("\t")[0].strip() for l in _open(_find(d, "features.tsv"))])
    return X, bc, genes


def cellset(p):
    with _open(p) as fh:
        return set(l.split("\t")[0].split("-")[0].strip() for l in fh if l.strip())


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rustar-gene", required=True)
    ap.add_argument("--rustar-genefull", required=True)
    ap.add_argument("--cellranger-raw", required=True)
    ap.add_argument("--rustar-ed-cells", required=True)
    ap.add_argument("--cr-ed-cells", required=True)
    ap.add_argument("--cr-native-cells", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    out = {}

    # ---- A. intron effect: rustar Gene vs GeneFull ----
    Xg, bcg, _ = load(a.rustar_gene)
    tot_gene = int(Xg.sum())
    g_by_bc = dict(zip(bcg, np.asarray(Xg.sum(1)).ravel()))
    del Xg
    Xf, bcf, genes_f = load(a.rustar_genefull)
    tot_genefull = int(Xf.sum())
    f_by_bc = dict(zip(bcf, np.asarray(Xf.sum(1)).ravel()))
    print(f"[A] intron effect (all barcodes):")
    print(f"    rustar Gene     total UMI = {tot_gene:,}")
    print(f"    rustar GeneFull total UMI = {tot_genefull:,}  "
          f"(+{100*(tot_genefull-tot_gene)/tot_gene:.1f}%)")
    out["intron_effect"] = {"gene_total_umi": tot_gene, "genefull_total_umi": tot_genefull,
                            "pct_increase": round(100*(tot_genefull-tot_gene)/tot_gene, 2)}

    # ---- B. raw parity: rustar GeneFull vs CellRanger raw ----
    Xc, bcc, genes_c = load(a.cellranger_raw)
    # 5' chemistry: CellRanger may report RC barcodes — detect against rustar.
    rust_set = set(bcf)
    ov_plain = len(set(bcc) & rust_set)
    bcc_rc = np.array([revcomp(b) for b in bcc])
    ov_rc = len(set(bcc_rc) & rust_set)
    if ov_rc > ov_plain:
        bcc = bcc_rc
        cr_orient = "reverse-complement"
    else:
        cr_orient = "as-reported"
    print(f"\n[B] CellRanger barcode orientation vs rustar: {cr_orient} "
          f"(overlap plain={ov_plain:,} rc={ov_rc:,})")
    tot_cr = int(Xc.sum())
    print(f"    rustar GeneFull raw total UMI = {tot_genefull:,}, genes={ (np.asarray(Xf.sum(0)).ravel()>0).sum():,}")
    print(f"    CellRanger      raw total UMI = {tot_cr:,}, genes={ (np.asarray(Xc.sum(0)).ravel()>0).sum():,}")
    c_by_bc = dict(zip(bcc, np.asarray(Xc.sum(1)).ravel()))
    out["raw_parity"] = {"rustar_genefull_total_umi": tot_genefull, "cellranger_total_umi": tot_cr,
                         "cr_orientation": cr_orient}

    # ---- C. cell-set agreement ----
    r_ed = cellset(a.rustar_ed_cells)
    cr_ed = cellset(a.cr_ed_cells)
    cr_nat = cellset(a.cr_native_cells)
    if cr_orient == "reverse-complement":
        cr_ed = {revcomp(b) for b in cr_ed}
        cr_nat = {revcomp(b) for b in cr_nat}

    def jac(x, y):
        i, u = len(x & y), len(x | y)
        return {"a": len(x), "b": len(y), "shared": i, "jaccard": round(i/u, 4) if u else None,
                "a_only": len(x - y), "b_only": len(y - x)}

    print("\n[C] cell-set agreement:")
    out["cell_sets"] = {}
    for label, x, y in [("rustar-ED vs CR-raw-ED (same algo)", r_ed, cr_ed),
                        ("rustar-ED vs CR-native-filtered", r_ed, cr_nat),
                        ("CR-raw-ED vs CR-native-filtered", cr_ed, cr_nat)]:
        j = jac(x, y)
        print(f"    {label:<38}: a={j['a']:,} b={j['b']:,} shared={j['shared']:,} "
              f"jaccard={j['jaccard']}")
        out["cell_sets"][label] = j

    # ---- D. per-cell UMI correlation on shared (rustar-ED ∩ CR-native) ----
    shared = sorted(r_ed & cr_nat)
    xs = [f_by_bc.get(b, 0) for b in shared]
    ys = [c_by_bc.get(b, 0) for b in shared]
    if len(shared) > 2:
        r = float(np.corrcoef(xs, ys)[0, 1])
        print(f"\n[D] per-cell UMI corr (rustar GeneFull vs CR raw) on {len(shared):,} shared "
              f"filtered cells: r={r:.4f}")
        out["per_cell_corr"] = {"shared_cells": len(shared), "pearson_r": round(r, 4)}

    # ---- write filtered h5ad ----
    def write_h5ad(name, X, bc, genes, keep_set):
        keep = np.array([b in keep_set for b in bc])
        Xk = X[keep]
        A = ad.AnnData(X=Xk, obs=pd.DataFrame(index=bc[keep]), var=pd.DataFrame(index=genes))
        A.obs["n_umi"] = np.asarray(Xk.sum(1)).ravel()
        A.obs["n_genes"] = np.asarray((Xk > 0).sum(1)).ravel()
        p = os.path.join(a.out, f"{name}.h5ad")
        A.write_h5ad(p)
        print(f"    wrote {p}  ({A.n_obs:,} cells)")
        return A.n_obs

    print("\n[E] writing EmptyDrops-filtered h5ad:")
    out["h5ad"] = {
        "rustar_genefull_ed": write_h5ad("rustar.GeneFull.emptydrops", Xf, bcf, genes_f, r_ed),
        "cellranger_native": write_h5ad("CellRanger.filtered", Xc, bcc, genes_c, cr_nat),
    }

    with open(os.path.join(a.out, "genefull_h5_compare.json"), "w") as fh:
        json.dump(out, fh, indent=2, default=str)
    print(f"\nwrote {a.out}/genefull_h5_compare.json")


if __name__ == "__main__":
    sys.exit(main())
