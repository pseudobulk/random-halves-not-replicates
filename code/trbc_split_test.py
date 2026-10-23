"""TRBC halves against random halves of the same cells (output of cd4_trbc_split.py).
Each cell's deviation is taken against the control cells of its own TRBC class (removes any expression difference
between TRBC1 and TRBC2 cells that is unrelated to the perturbation); genes: control mean >= 0.1, TCR genes excluded.
Per target with >= MINH classed cells in each class and in each random half:
  trbc     class 0 cells against class 1 cells (no shared clone)
  random   the same classed cells split by r1 (and, as a second draw, by r2); both halves mix the classes
  cross-half signal sum_g x_A,g x_B,g
Targets: inert (own-gene control mean < 0.01) and expressed (>= 0.5).
Tests across targets: share positive; mean signal; paired sign test and Wilcoxon signed-rank test of random - trbc
(one-sided, random > trbc). Expected if the shared component is clonal: inert targets positive in about half of the
cases under TRBC halves, and random > trbc for most targets.
Usage: trbc_split_test.py <dir> <label> <out.json>"""
import sys, json
import numpy as np
from scipy.stats import binomtest, wilcoxon

MINH = 10


def main():
    d, label, outp = sys.argv[1:4]
    meta = json.load(open(d + "/meta.json")); genes, groups = meta["genes"], meta["groups"]
    S = np.load(d + "/sum8.npy", mmap_mode="r"); n8 = np.load(d + "/n8.npy").astype(float)
    ci = len(groups) - 1
    Sc = np.asarray(S[ci], dtype=np.float64); nc = n8[ci]
    cm = Sc.sum(0) / nc.sum()
    tcr = {j for j, g in enumerate(genes) if g.startswith(("TRBC", "TRAC", "TRBV", "TRAV", "TRGC", "TRDC"))}
    gk = np.array([j for j in np.where(cm >= 0.1)[0] if j not in tcr])
    gpos = {g: j for j, g in enumerate(genes)}
    mu = [Sc[0:4][:, gk].sum(0) / nc[0:4].sum(), Sc[4:8][:, gk].sum(0) / nc[4:8].sum()]
    C0, C1 = [0, 1, 2, 3], [4, 5, 6, 7]
    sides = {"r1": ([0, 1, 4, 5], [2, 3, 6, 7]), "r2": ([0, 2, 4, 6], [1, 3, 5, 7])}

    def eff(St, nt, codes):
        x = np.zeros(len(gk)); n = 0.0
        for c in codes:
            mc = mu[0] if c < 4 else mu[1]
            x += St[c] - nt[c] * mc; n += nt[c]
        return x / n

    rows = []
    for t in range(ci):
        g = groups[t]
        if g not in gpos:
            continue
        nt = n8[t]
        if nt[C0].sum() < MINH or nt[C1].sum() < MINH or min(nt[s].sum() for sd in sides.values() for s in sd) < MINH:
            continue
        St = np.asarray(S[t][:, gk], dtype=np.float64)
        r = dict(target=g, expr=float(cm[gpos[g]]), n0=float(nt[C0].sum()), n1=float(nt[C1].sum()),
                 trbc=float(eff(St, nt, C0) @ eff(St, nt, C1)))
        for k, (a, b) in sides.items():
            r[k] = float(eff(St, nt, a) @ eff(St, nt, b))
        rows.append(r)
    out = {"rules": __doc__, "label": label, "n_genes": int(len(gk)), "classed_cells": [meta["class0"], meta["class1"]],
           "control_class_difference_energy": float(((mu[0] - mu[1]) ** 2).sum()), "groups": {}}
    for name, sel in (("inert", lambda r: r["expr"] < 0.01), ("expressed", lambda r: r["expr"] >= 0.5)):
        R = [r for r in rows if sel(r)]
        if not R:
            continue
        tr = np.array([r["trbc"] for r in R]); r1 = np.array([r["r1"] for r in R]); r2 = np.array([r["r2"] for r in R])
        d1 = r1 - tr
        out["groups"][name] = dict(
            n_targets=len(R), median_cells_class0=float(np.median([r["n0"] for r in R])), median_cells_class1=float(np.median([r["n1"] for r in R])),
            share_pos_trbc=float((tr > 0).mean()), share_pos_random_r1=float((r1 > 0).mean()), share_pos_random_r2=float((r2 > 0).mean()),
            mean_trbc=float(tr.mean()), mean_random_r1=float(r1.mean()), mean_random_r2=float(r2.mean()),
            share_random_gt_trbc=float((d1 > 0).mean()),
            sign_test_p=float(binomtest(int((d1 > 0).sum()), len(d1), 0.5, alternative="greater").pvalue),
            wilcoxon_p=float(wilcoxon(r1, tr, alternative="greater").pvalue),
            sign_test_trbc_pos_vs_half_p=float(binomtest(int((tr > 0).sum()), len(tr), 0.5).pvalue))
        print(label, name, json.dumps(out["groups"][name]), flush=True)
    out["targets"] = rows
    json.dump(out, open(outp, "w"), indent=1)


if __name__ == "__main__":
    main()
