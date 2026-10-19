"""Direct clone tests on CD4 single cells (preprint review 2.1 and 2.3; A100; output of cd4_sc_extract_cells.py).
Per state (Rest, Stim48hr), on targets whose gene is measured:
DEPTH  Designs from the 8-way aggregates (code = gp * 4 + r1 * 2 + r2):
         random_half  r1 = 0 vs 1            random_quarter  (r1 = 0, r2 = 0) vs (r1 = 1, r2 = 0)
         guide_half   gp = 0 vs 1            guide_quarter   (gp = 0, r2 = 0) vs (gp = 1, r2 = 0)
       Effect = subset mean - control mean of the same subset (log1p CP10K; genes with control mean >= 0.1).
       Same targets for all designs (>= MINQ cells in every quarter side). For inert (control mean < 0.01) and
       expressed (>= 0.5) targets: cross energy per target (referenced to the inert mean, Z1 an), frac_pos,
       mean Pearson between sides and its Spearman-Brown value.
       Expected under the unit model (Z1 ap): random designs - inert cross energy about equal at quarter and half
       depth (the unit term does not shrink), Pearson / frac_pos rise with depth; guide designs - about 0 at both.
TRBC   Per cell TRBC1 fraction f = TRBC1 / (TRBC1 + TRBC2) (cells with >= MINTR counts). A clone carries one TCR beta
       rearrangement, so cells of one clone share f; a perturbation drawn from few clones has a mean f away from the
       population mean.
         overdispersion   var over targets of z = (mean f - f0) / (sd_f / sqrt(n)); 1 under random assignment;
                          compared with a label permutation of the same cells
         reproducibility  Pearson over targets of (mean f - f0) between sides, random halves vs guide halves
         link             Spearman over inert targets between the TRBC skew |z| and the expression cross energy of
                          the random halves (and the TRBC cross product d_A d_B vs expression cross energy)
Usage: clone_tests_cd4.py <cells_dir_state1> <label1> [<cells_dir_state2> <label2>] <out_dir>. Exploratory."""
import os, sys, json
import numpy as np
from scipy.stats import spearmanr

MINQ, MINTR, MINF = 8, 2, 15
DESIGNS = {"random_half": ([0, 1, 4, 5], [2, 3, 6, 7]), "random_quarter": ([0, 4], [2, 6]),
           "guide_half": ([0, 1, 2, 3], [4, 5, 6, 7]), "guide_quarter": ([0, 2], [4, 6])}


def sb(r):
    return float(2 * r / (1 + r)) if r > 0 else float("nan")


def pearson_rows(A, B):
    a = A - A.mean(1, keepdims=True); b = B - B.mean(1, keepdims=True)
    return (a * b).sum(1) / np.sqrt((a ** 2).sum(1) * (b ** 2).sum(1))


def run_state(d, label):
    meta = json.load(open(d + "/meta.json")); genes, groups = meta["genes"], meta["groups"]
    S = np.load(d + "/sum8.npy", mmap_mode="r"); n8 = np.load(d + "/n8.npy")
    ci = len(groups) - 1
    cm = np.asarray(S[ci].sum(0)) / n8[ci].sum()
    trg = {j for j, g in enumerate(genes) if g.startswith(("TRBC", "TRAC", "TRBV", "TRAV", "TRGC", "TRDC"))}
    gk = np.array([j for j in np.where(cm >= 0.1)[0] if j not in trg])   # TCR genes excluded (10-18 audit: circular link)
    gpos = {g: j for j, g in enumerate(genes)}; gkpos = {j: k for k, j in enumerate(gk)}
    texp = {t: cm[gpos[t]] for t in groups[:-1] if t in gpos}
    quarters = [0, 4, 2, 6]
    tids = [i for i, t in enumerate(groups[:-1]) if t in texp and n8[i, quarters].min() >= MINQ]
    tnames = [groups[i] for i in tids]
    inert = np.array([texp[t] < 0.01 for t in tnames]); expd = np.array([texp[t] >= 0.5 for t in tnames])
    res = {"label": label, "n_targets": len(tids), "n_inert": int(inert.sum()), "n_expressed": int(expd.sum()), "depth": {}}
    Sg = np.asarray(S[tids][:, :, gk], dtype=np.float64); Sc = np.asarray(S[ci][:, gk], dtype=np.float64)
    nt = n8[tids].astype(float); nc = n8[ci].astype(float)
    eff = {}
    for nm, (a, b) in DESIGNS.items():
        sides = []
        for codes in (a, b):
            m = Sg[:, codes].sum(1) / nt[:, codes].sum(1)[:, None] - Sc[codes].sum(0) / nc[codes].sum()
            for k, i in enumerate(tids):
                j = gkpos.get(gpos.get(groups[i], -1))
                if j is not None:
                    m[k, j] = 0
            sides.append(m)
        eff[nm] = sides
        A, B = sides
        Ai, Bi = A - A[inert].mean(0), B - B[inert].mean(0)
        cross = (Ai * Bi).sum(1); pr = pearson_rows(A, B)
        r = {"median_cells_per_side": [float(np.median(nt[:, a].sum(1))), float(np.median(nt[:, b].sum(1)))]}
        for sname, mm in (("inert", inert), ("expressed", expd)):
            r[sname] = dict(cross_per_target=float(cross[mm].mean()), frac_pos=float((cross[mm] > 0).mean()),
                            pearson_mean=float(np.nanmean(pr[mm])), pearson_sb=sb(np.nanmean(pr[mm])))
        res["depth"][nm] = r
        print(label, nm, json.dumps(r), flush=True)
    # TRBC
    c = np.load(d + "/cells.npz")
    grp, gp, r1 = c["grp"], c["gp"], c["r1"]; t1, t2 = c["trbc1"], c["trbc2"]
    tot = t1 + t2; ok = (grp >= 0) & (c["umi"] > 0) & (tot >= MINTR)
    f = np.where(tot > 0, t1 / np.maximum(tot, 1), np.nan)
    res["trbc"] = {"cells_valid": int(ok.sum()), "cells_total": int(((grp >= 0) & (c["umi"] > 0)).sum()),
                   "frac_only_trbc1": float(((t1 > 0) & (t2 == 0))[ok].mean()), "frac_only_trbc2": float(((t2 > 0) & (t1 == 0))[ok].mean()),
                   "f_mean": float(np.nanmean(f[ok])), "f_sd": float(np.nanstd(f[ok]))}
    f0, sdf = float(np.nanmean(f[ok])), float(np.nanstd(f[ok]))
    G = len(groups)
    def group_means(mask):
        s = np.bincount(grp[mask], weights=f[mask], minlength=G); n = np.bincount(grp[mask], minlength=G)
        return s / np.maximum(n, 1), n
    mf, nf = group_means(ok)
    tset = np.array([i for i in range(G - 1) if nf[i] >= 2 * MINF])
    z = (mf[tset] - f0) / (sdf / np.sqrt(nf[tset]))
    rng = np.random.default_rng(3); zn = []
    vg = grp[ok].copy(); vf = f[ok]
    for _ in range(5):
        pg = rng.permutation(vg); s = np.bincount(pg, weights=vf, minlength=G); n = np.bincount(pg, minlength=G)
        zn.append(np.var(((s / np.maximum(n, 1))[tset] - f0) / (sdf / np.sqrt(n[tset]))))
    res["trbc"]["overdispersion"] = dict(n_targets=int(len(tset)), var_z=float(np.var(z)), var_z_permuted=float(np.mean(zn)))
    for nm, (ma, mb) in (("random_half", (r1 == 0, r1 == 1)), ("guide_half", (gp == 0, gp == 1))):
        (fa, na), (fb, nb) = group_means(ok & ma), group_means(ok & mb)
        sel = np.array([i for i in range(G - 1) if na[i] >= MINF and nb[i] >= MINF])
        da, db = fa[sel] - f0, fb[sel] - f0
        tn = [groups[i] for i in sel]; ine = np.array([t in texp and texp[t] < 0.01 for t in tn])
        res["trbc"]["reproducibility_" + nm] = dict(n=int(len(sel)), pearson=float(np.corrcoef(da, db)[0, 1]),
                                                     pearson_inert=float(np.corrcoef(da[ine], db[ine])[0, 1]) if ine.sum() > 10 else None)
    # link: TRBC skew vs expression clone signal on inert targets (random halves)
    A, B = eff["random_half"]; Ai, Bi = A - A[inert].mean(0), B - B[inert].mean(0); ecross = (Ai * Bi).sum(1)
    li = [k for k, i in enumerate(tids) if inert[k] and nf[i] >= 2 * MINF]
    skew = np.array([abs(mf[tids[k]] - f0) / (sdf / np.sqrt(nf[tids[k]])) for k in li])
    (fa, na), (fb, nb) = group_means(ok & (r1 == 0)), group_means(ok & (r1 == 1))
    tcross = np.array([(fa[tids[k]] - f0) * (fb[tids[k]] - f0) for k in li])
    q = np.quantile(skew, [0.25, 0.75])
    res["trbc"]["link_inert"] = dict(n=len(li), spearman_skew_vs_expr_cross=float(spearmanr(skew, ecross[li])[0]),
                                     spearman_trbc_cross_vs_expr_cross=float(spearmanr(tcross, ecross[li])[0]),
                                     expr_cross_top_quartile_skew=float(np.mean(ecross[li][skew >= q[1]])),
                                     expr_cross_bottom_quartile_skew=float(np.mean(ecross[li][skew <= q[0]])))
    print(label, "trbc", json.dumps(res["trbc"]), flush=True)
    return res


def main():
    args = sys.argv[1:]; out = args[-1]; os.makedirs(out, exist_ok=True)
    res = {"rules": __doc__}
    for i in range(0, len(args) - 1, 2):
        res[args[i + 1]] = run_state(args[i], args[i + 1])
    json.dump(res, open(out + "/metrics_clone_tests_cd4.json", "w"), indent=1)
    print("saved", flush=True)


if __name__ == "__main__":
    main()
