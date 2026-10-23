"""Jin et al. 2020 (Science; GEO GSE157977): in vivo Perturb-seq in developing mouse cortex. Lentivirus at E12.5 into
dividing neural progenitors, so the cells of one perturbation descend from a limited number of infected progenitors
(clones). Both gRNAs of a gene sit in one vector, so the independent replicate here is another library: each of the 18
libraries pools BFP+ cells of 1-3 litters (4-20 embryos) harvested on one day, i.e. separate transductions.
Cells: dial-out counts (per library), assigned when the top perturbation barcode has >= 2 UMIs and >= 2x the second;
barcodes mapped to genes with Table S5 (barcodes.csv; 'GFP' = GFP-targeting control). RNA: raw CellRanger h5,
cells with >= MINUMI UMIs; log1p CP10K. Classes by marker means (argmax, >= 0.5): excitatory neurons, interneurons,
astrocytes, oligodendrocyte lineage.
Effects: every cell of a perturbation minus the mean of the OTHER assigned cells of its class in its library (removes
library batch; leaving the perturbation's own cells out of the reference avoids a negative bias between its halves,
which the first version had: with small libraries the noise of one half entered the other half's reference).
Per perturbation and class (>= MINC cells):
  random   cells split by a random bit; cross-half signal sum_g x_A,g x_B,g, averaged over NSPLIT random splits
           (sex-linked genes excluded: embryos of both sexes are pooled); SE by bootstrap over perturbations
  library  libraries split into two random groups (NSPLIT draws; >= 10 cells in each); same cross-half signal
The difference random - library is the part of the reproducible signal that only cells of the same libraries share:
clones, plus differences between the embryos pooled within a library (no per-embryo label, so this is an upper bound
on the clone part). Scaled by the summed per-cell variance of the class (comparable to vc_expr).
Usage: jin_clone.py <dir with GEO files + barcodes.csv> <out.json>"""
import os, sys, glob, gzip, json, re
import numpy as np, pandas as pd, h5py
from scipy.sparse import csc_matrix

MINUMI, MINC, NSPLIT, SEED = 1000, 40, 50, 0
MARKERS = {"ExN": ["Neurod6", "Slc17a7", "Satb2", "Tbr1", "Neurod2"], "InN": ["Gad1", "Gad2", "Dlx1", "Dlx2", "Dlx6os1"],
           "Astro": ["Aqp4", "Gja1", "Aldh1l1", "Slc1a3"], "Oligo": ["Olig1", "Olig2", "Pdgfra", "Sox10", "Mbp"]}


def load_library(d, s, key):
    dl = glob.glob(os.path.join(d, "*dialout.sample%d.UMI*" % s))[0]
    h5 = glob.glob(os.path.join(d, "*UMI.counts.sample%d.h5" % s))[0]
    D = pd.read_csv(dl, index_col=0)
    D = D[[c for c in D.columns if c in key]]
    v = D.to_numpy(); o = np.sort(v, 1)[:, ::-1]
    ok = (o[:, 0] >= 2) & (o[:, 0] >= 2 * o[:, 1])
    pert = pd.Series(np.array([key[c] for c in D.columns])[v.argmax(1)], index=D.index)[ok]
    with h5py.File(h5, "r") as f:
        g = f["mm10"]; bc = np.array([x.decode().split("-")[0] for x in g["barcodes"][:]])
        genes = [x.decode() for x in g["gene_names"][:]]
        sh = tuple(g["shape"][:])
        M = csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=sh)   # genes x droplets
    pos = {b: i for i, b in enumerate(bc)}
    cells = [c for c in pert.index if c in pos]
    X = M[:, [pos[c] for c in cells]].T.tocsr().astype(np.float32)
    umi = np.asarray(X.sum(1)).ravel(); keep = umi >= MINUMI
    X = X[keep]; umi = umi[keep]
    X.data = np.log1p(X.data / np.repeat(umi, np.diff(X.indptr)) * 1e4)
    return genes, X, pert.loc[np.array(cells)[keep]].to_numpy()


def main():
    d, outp = sys.argv[1:3]; rng = np.random.default_rng(SEED)
    key = dict((r.barcode, r.gene) for r in pd.read_csv(os.path.join(d, "barcodes.csv")).itertuples())
    libs = sorted(int(re.search(r"sample(\d+)\.h5", p).group(1)) for p in glob.glob(os.path.join(d, "*UMI.counts.sample*.h5")))
    Xs, P, L = [], [], []
    for s in libs:
        genes, X, p = load_library(d, s, key)
        Xs.append(X); P.append(p); L.append(np.full(len(p), s))
        print("library", s, "cells", X.shape[0], flush=True)
    from scipy.sparse import vstack
    X = vstack(Xs).tocsr(); P = np.concatenate(P); L = np.concatenate(L)
    gi = {g: j for j, g in enumerate(genes)}
    score = np.stack([np.asarray(X[:, [gi[g] for g in m if g in gi]].mean(1)).ravel() for m in MARKERS.values()], 1)
    cls = np.where(score.max(1) >= 0.5, score.argmax(1), -1); names = list(MARKERS)
    out = {"rules": __doc__, "n_cells": int(X.shape[0]), "libraries": libs,
           "cells_per_class": {names[c]: int((cls == c).sum()) for c in range(len(names))}, "unassigned": int((cls < 0).sum()),
           "cells_per_perturbation": {k: int(v) for k, v in pd.Series(P).value_counts().items()}, "classes": {}}
    print(out["cells_per_class"], "unassigned", out["unassigned"], flush=True)
    for c, cname in enumerate(names):
        m = cls == c
        if m.sum() < 500:
            continue
        Xc = X[m].toarray(); Pc, Lc = P[m], L[m]
        SEX = {"Xist", "Tsix", "Ddx3y", "Eif2s3y", "Kdm5d", "Uty"}     # embryos of both sexes are pooled
        gm = (Xc.mean(0) >= 0.1) & np.array([g not in SEX for g in genes])
        Xc = Xc[:, gm]
        libs_c = np.unique(Lc)
        SL = {lib: Xc[Lc == lib].sum(0) for lib in libs_c}; NL = {lib: int((Lc == lib).sum()) for lib in libs_c}
        var = float(sum(((Xc[Lc == lib] - SL[lib] / NL[lib]) ** 2).sum(0) for lib in libs_c).sum() / len(Xc))
        res = {}
        X0 = Xc
        for p in np.unique(Pc):
            idx = np.where(Pc == p)[0]
            if len(idx) < MINC:
                continue
            # reference: the other cells of the same library and class (leave the perturbation out, so that the noise
            # of one half does not enter the reference of the other half)
            Xc = X0.copy()
            for lib in np.unique(Lc[idx]):
                own = idx[Lc[idx] == lib]; n_o = NL[lib] - len(own)
                Xc[own] -= (SL[lib] - X0[own].sum(0)) / max(n_o, 1)
            crs = []
            for _ in range(NSPLIT):                                   # average over random splits
                bit = rng.integers(0, 2, len(idx)); a, b = idx[bit == 0], idx[bit == 1]
                if len(a) >= 10 and len(b) >= 10:
                    crs.append(float(Xc[a].mean(0) @ Xc[b].mean(0)))
            cr = float(np.mean(crs)) if crs else np.nan
            lib_u = np.unique(Lc[idx]); cl = []
            for _ in range(NSPLIT):
                grp = set(rng.permutation(lib_u)[: len(lib_u) // 2])
                ga = idx[np.isin(Lc[idx], list(grp))]; gb = idx[~np.isin(Lc[idx], list(grp))]
                if len(ga) >= 10 and len(gb) >= 10:
                    cl.append(float(Xc[ga].mean(0) @ Xc[gb].mean(0)))
            res[p] = dict(n=int(len(idx)), n_libraries=int(len(lib_u)), random=cr, library=float(np.mean(cl)) if cl else np.nan,
                          n_splits=len(cl))
        tg = {p: r for p, r in res.items() if p != "GFP" and np.isfinite(r["random"]) and np.isfinite(r["library"])}
        rnd = np.array([r["random"] for r in tg.values()]); lib = np.array([r["library"] for r in tg.values()])
        summ = dict(n_cells=int(m.sum()), n_genes=int(gm.sum()), per_cell_var_sum=float(var), n_perturbations=len(tg),
                    mean_random=float(rnd.mean()), mean_library=float(lib.mean()),
                    clone_part_share=float(1 - lib.mean() / rnd.mean()),
                    clone_strength=float((rnd - lib).mean() / var),
                    clone_strength_se=float(np.std([rng.choice(rnd - lib, len(rnd)).mean() for _ in range(2000)]) / var),
                    positive_random_minus_library=float(((rnd - lib) > 0).mean()))
        if "GFP" in res:
            g = res["GFP"]; summ["GFP"] = dict(g, random_scaled=g["random"] / var, library_scaled=g["library"] / var,
                                               ratio_to_targeting_random=g["random"] / rnd.mean())
        out["classes"][cname] = dict(summary=summ, perturbations=res)
        print(cname, json.dumps(summ), flush=True)
    json.dump(out, open(outp, "w"), indent=1)


if __name__ == "__main__":
    main()
