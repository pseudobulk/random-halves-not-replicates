"""Confidence intervals for the headline numbers (preprint writing step; A100). Readout predictions are computed once
(5-fold CV, verify_headline.cv_readout) and held fixed; the bootstrap resamples the units (targets for CD4, conditions
for Tahoe, 1,000 draws) and recomputes ratio-of-sums statistics. Split-to-split variation is shown separately.
cd4   donor-1 per-cell aggregates after the audit fixes (cells_D1_*_v2; code = gp * 4 + r1 * 2 + r2). Designs:
        random_r1, random_r2, random_xor   three different random halvings of the same cells
        guide                              gp = 0 vs 1 (independent transductions)
      Targets present in both states with >= 10 cells on every side of every design, target gene measured; genes with
      control mean >= 0.1 in both states; inert = control mean < 0.01 in both states; expressed >= 0.5 in both.
      Per design: inert frac_pos per state (inert baseline), gamma share (template removed), readout share k = 20 and
      a shuffled control, mean per-target Pearson (rest) and its Spearman-Brown value, and the share of targets that
      Wang et al. 2026 would call reliable (rho = SB of the per-target Pearson, median over the three random halvings for
      the random design; rho >= 0.5), for all / inert / expressed targets.
tahoe fixed-target readout (k = 10, 5-fold) on the plate 6 / 14 arrays of tahoe_plate_recompute.py, beta from the same
      well vs the other plate, and the across-plate gamma share; bootstrap over conditions.
Usage: ci_bootstrap.py cd4 <cells_Rest> <cells_Stim> <out_dir> | tahoe <D_arrays.npz> <out_dir>. Exploratory."""
import os, sys, json
import numpy as np
from verify_headline import cv_readout

NB = 1000
DES = {"random_r1": ([0, 1, 4, 5], [2, 3, 6, 7]), "random_r2": ([0, 2, 4, 6], [1, 3, 5, 7]),
       "random_xor": ([0, 3, 4, 7], [1, 2, 5, 6]), "guide": ([0, 1, 2, 3], [4, 5, 6, 7])}


def sb(r):
    r = np.asarray(r, float); return np.where(r > 0, 2 * r / (1 + r), np.nan)


def ci(fun, n, rng):
    v = [fun(rng.integers(0, n, n)) for _ in range(NB)]
    return [float(np.nanpercentile(v, 2.5)), float(np.nanpercentile(v, 97.5))]


def prow(A, B):
    a = A - A.mean(1, keepdims=True); b = B - B.mean(1, keepdims=True)
    return (a * b).sum(1) / np.sqrt((a ** 2).sum(1) * (b ** 2).sum(1))


def cd4(dr, ds, out):
    def load(d):
        m = json.load(open(d + "/meta.json")); return m["genes"], m["groups"], np.load(d + "/sum8.npy", mmap_mode="r"), np.load(d + "/n8.npy")
    st = {"rest": load(dr), "stim": load(ds)}
    genes = sorted(set(st["rest"][0]) & set(st["stim"][0])); gpos = {g: j for j, g in enumerate(genes)}
    ok, cms = None, {}
    for k, (gl, groups, S, n8) in st.items():
        gi = {g: j for j, g in enumerate(gl)}; col = np.array([gi[g] for g in genes]); ci_ = len(groups) - 1
        cms[k] = np.asarray(S[ci_].sum(0))[col] / n8[ci_].sum()
        s_ = {t for i, t in enumerate(groups[:-1]) if all(n8[i, c].sum() >= 10 for a, b in DES.values() for c in (a, b))}
        ok = s_ if ok is None else ok & s_
    targets = sorted(t for t in ok if t in gpos)
    keep = np.where(np.minimum(cms["rest"], cms["stim"]) >= 0.1)[0]; kpos = {j: q for q, j in enumerate(keep)}
    e = np.array([max(cms["rest"][gpos[t]], cms["stim"][gpos[t]]) for t in targets]); inert = e < 0.01
    expd = np.array([min(cms["rest"][gpos[t]], cms["stim"][gpos[t]]) for t in targets]) >= 0.5
    zc = np.array([kpos.get(gpos[t], -1) for t in targets])
    data = {}
    for k, (gl, groups, S, n8) in st.items():
        gi = {g: j for j, g in enumerate(gl)}; col = np.array([gi[g] for g in genes])[keep]; ci_ = len(groups) - 1
        ti = {t: i for i, t in enumerate(groups)}; idx = np.array([ti[t] for t in targets])
        data[k] = (np.asarray(S[idx])[:, :, col].astype(np.float64), n8[idx].astype(float), np.asarray(S[ci_])[:, col].astype(np.float64), n8[ci_].astype(float))
    res = {"rules": __doc__, "n_targets": len(targets), "n_inert": int(inert.sum()), "n_expressed": int(expd.sum()), "n_genes": int(len(keep))}
    rng = np.random.default_rng(13); P = len(targets); rows = np.arange(P)
    per_target_r = {}
    for dn, (a, b) in DES.items():
        H = {}
        for k in st:
            Sg, nt, Sc, nc = data[k]; sides = []
            for codes in (a, b):
                X = Sg[:, codes].sum(1) / nt[:, codes].sum(1)[:, None] - Sc[codes].sum(0) / nc[codes].sum()
                m = zc >= 0; X[rows[m], zc[m]] = 0; sides.append(X)
            H[k] = sides
        d = {}
        for k in st:
            A, B = H[k]; Ai, Bi = A - A[inert].mean(0), B - B[inert].mean(0); x = (Ai * Bi).sum(1)
            xi = x[inert]; d["inert_frac_pos_" + k] = [float((xi > 0).mean())] + ci(lambda ix: (xi[ix] > 0).mean(), len(xi), rng)
        r = prow(*H["rest"]); per_target_r[dn] = r
        for sname, msk in (("all", np.ones(P, bool)), ("inert", inert), ("expressed", expd)):
            rr = r[msk]
            d["pearson_rest_" + sname] = [float(np.nanmean(rr))] + ci(lambda ix: np.nanmean(rr[ix]), len(rr), rng)
            d["pearson_rest_sb_" + sname] = float(sb(np.nanmean(rr)))
        c = lambda M: M - M.mean(0)
        XA, XB, YA, YB = c(H["rest"][0]), c(H["rest"][1]), c(H["stim"][0]), c(H["stim"][1])
        bA, gA, gB, bB = (XA + YA) / 2, (XA - YA) / 2, (XB - YB) / 2, (XB + YB) / 2
        gs, bs = (gA * gB).sum(1), (bA * bB).sum(1)
        d["gamma_share"] = [float(gs.sum() / (gs.sum() + bs.sum()))] + ci(lambda ix: gs[ix].sum() / (gs[ix].sum() + bs[ix].sum()), P, rng)
        F = cv_readout(bA, gB, 20, 7); Fs = cv_readout(bA, gB, 20, 7, shuffle=True)
        ex = (gB ** 2).sum(1) - ((gB - F) ** 2).sum(1); exs = (gB ** 2).sum(1) - ((gB - Fs) ** 2).sum(1)
        d["readout_k20"] = [float(ex.sum() / gs.sum())] + ci(lambda ix: ex[ix].sum() / gs[ix].sum(), P, rng)
        d["readout_k20_shuffled"] = float(exs.sum() / gs.sum())
        d["gamma_cross_per_target"] = float(gs.mean()); d["beta_cross_per_target"] = float(bs.mean())
        res[dn] = d
        print(dn, json.dumps(d), flush=True)
    rel = {"random": sb(np.median(np.stack([per_target_r[x] for x in ("random_r1", "random_r2", "random_xor")]), 0)), "guide": sb(per_target_r["guide"])}
    res["wang_reliable_share_rest"] = {dn: {s: float(np.nanmean(np.nan_to_num(rel[dn][m], nan=0.0) >= 0.5)) for s, m in (("all", np.ones(P, bool)), ("inert", inert), ("expressed", expd))} for dn in rel}
    print("wang", json.dumps(res["wang_reliable_share_rest"]), flush=True)
    return res


def tahoe(npz, out):
    z = np.load(npz); X = {k: z[k].astype(np.float64) for k in ("plate6_0", "plate6_1", "plate6_pool", "plate14_0", "plate14_1", "plate14_pool")}
    def twoway(M):
        m = M.mean((0, 1)); return M - M.mean(1, keepdims=True) - M.mean(0, keepdims=True) + m, M.mean(0) - m
    C, Pn, G = X["plate6_0"].shape
    res = {"rules": __doc__, "n_lines": C, "n_conditions": Pn}; rng = np.random.default_rng(17)
    for Pl, O in (("plate14", "plate6"), ("plate6", "plate14")):
        g1, _ = twoway(X[Pl + "_1"]); g0, b0 = twoway(X[Pl + "_0"]); _, bo = twoway(X[O + "_pool"])
        den = (g0 * g1).sum((0, 2))                          # per condition, summed over lines
        d = {}
        for nm, Fe in (("same_well", b0), ("other_plate", bo)):
            ex = np.zeros(Pn)
            for c in range(C):
                Pc = cv_readout(Fe, g1[c], 10, 11 + c); ex += (g1[c] ** 2).sum(1) - ((g1[c] - Pc) ** 2).sum(1)
            d[nm] = [float(ex.sum() / den.sum())] + ci(lambda ix: ex[ix].sum() / den[ix].sum(), Pn, rng)
        res["readout_target_" + Pl] = d
        print(Pl, json.dumps(d), flush=True)
    gA, bA = twoway(X["plate6_pool"]); gB, bB = twoway(X["plate14_pool"])
    gs, bs = (gA * gB).sum((0, 2)), (bA * bB).sum(1) * C
    res["gamma_share_across_plates"] = [float(gs.sum() / (gs.sum() + bs.sum()))] + ci(lambda ix: gs[ix].sum() / (gs[ix].sum() + bs[ix].sum()), Pn, rng)
    for p in ("plate6", "plate14"):
        g0, b0 = twoway(X[p + "_0"]); g1, b1 = twoway(X[p + "_1"])
        gs, bs = (g0 * g1).sum((0, 2)), (b0 * b1).sum(1) * C
        res["gamma_share_within_well_" + p] = [float(gs.sum() / (gs.sum() + bs.sum()))] + ci(lambda ix: gs[ix].sum() / (gs[ix].sum() + bs[ix].sum()), Pn, rng)
    print("gamma", json.dumps({k: v for k, v in res.items() if k.startswith("gamma")}), flush=True)
    return res


def main():
    mode = sys.argv[1]
    if mode == "cd4":
        out = sys.argv[4]; os.makedirs(out, exist_ok=True); res = cd4(*sys.argv[2:5])
    else:
        out = sys.argv[3]; os.makedirs(out, exist_ok=True); res = tahoe(*sys.argv[2:4])
    json.dump(res, open(out + "/metrics_ci_bootstrap.json", "w"), indent=1)
    print("saved", flush=True)


if __name__ == "__main__":
    main()
