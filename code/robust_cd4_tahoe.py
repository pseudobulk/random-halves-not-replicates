"""Robustness of Z1 (am) and (aj) (A100).
cd4   donor-1 single cells, Rest vs Stim48hr, random cell halves vs guide halves, now on the SAME targets (the
      (am) table used 9,623 vs 6,673 targets). Grid: k in (5, 10, 20, 50); gene set all expressed (control >= 0.1) vs
      top 2,000 by control mean; target subsets all common vs guide-balanced (the two guide halves' cell counts within
      a factor 1.5 of each other in both states). Per cell: gamma share, readout share (state = top-k PCs of half-A
      beta, ridge with lambda from inner 2-fold CV, 2 folds over targets) and its shuffled control, inert targets'
      frac_pos (referenced to the inert baseline, Z1 an). Template removed for the readout.
tahoe target fixed (two-way gamma of plate P half 1), beta from the same well vs the other plate, as
      tahoe_readout_source.py, now with k in (5, 10, 20), lambda by inner CV instead of fixed 1, and gene set all vs
      top 2,000 by mean absolute effect on the OTHER plate (so the gene selection does not reuse the scored plate).
Usage: robust_cd4_tahoe.py cd4 <out_Rest> <out_Stim> <out_Rest_guide> <out_Stim_guide> <out_dir>
       robust_cd4_tahoe.py tahoe <D_arrays.npz> <out_dir>. Exploratory."""
import os, sys, json
import numpy as np

LAMS = (1e-3, 1e-2, 1e-1, 1, 10, 100, 1000)


def ridge(Z, Y, lam):
    k = Z.shape[1]
    return np.linalg.solve(Z.T @ Z + lam * np.trace(Z.T @ Z) / k * np.eye(k), Z.T @ Y)


def pick_lam(Z, Y, rng):
    n = len(Z); perm = rng.permutation(n); f = (perm[: n // 2], perm[n // 2:]); best, bl = -np.inf, 1.0
    for lam in LAMS:
        e = sum(((Y[f[1 - i]] ** 2).sum() - ((Y[f[1 - i]] - Z[f[1 - i]] @ ridge(Z[f[i]], Y[f[i]], lam)) ** 2).sum()) for i in (0, 1))
        if e > best:
            best, bl = e, lam
    return bl


def predict(bA, gB, k, seed=11, shuffle=False):
    P = len(gB); rng = np.random.default_rng(seed); perm = rng.permutation(P); folds = (perm[: P // 2], perm[P // 2:])
    F = np.zeros_like(gB); src = bA[rng.permutation(P)] if shuffle else bA
    for i in (0, 1):
        tr, te = folds[i], folds[1 - i]
        V = np.linalg.svd(src[tr], full_matrices=False)[2][:k].T; Z = src @ V
        F[te] = Z[te] @ ridge(Z[tr], gB[tr], pick_lam(Z[tr], gB[tr], rng))
    return F


def expl(F, g):
    return float((g ** 2).sum() - ((g - F) ** 2).sum())


def run_cd4(dr, ds, drg, dsg, out):
    from cd4_readout import load, MIN_HALF, MIN_EXPR
    D = {"cell": (load(dr), load(ds)), "guide": (load(drg), load(dsg))}
    gst = {"rest": list(D["cell"][0][1]), "stim": list(D["cell"][1][1])}
    for nm, (r, s_) in D.items():
        assert list(r[1]) == gst["rest"] and list(s_[1]) == gst["stim"]
    genes = sorted(set(gst["rest"]) & set(gst["stim"]))
    col = {st: np.array([{g: j for j, g in enumerate(gst[st])}[g] for g in genes]) for st in gst}
    ctrl_mean = {}
    for nm, (r, s_) in D.items():
        for st, (groups, _, S, n) in (("rest", r), ("stim", s_)):
            ci = groups.index("__ctrl__"); ctrl_mean[(nm, st)] = np.stack([np.asarray(S[ci, h])[col[st]] / n[ci, h] for h in range(2)])
    cm = (ctrl_mean[("cell", "rest")].mean(0) + ctrl_mean[("cell", "stim")].mean(0)) / 2
    ok = None
    for nm, (r, s_) in D.items():
        for groups, _, S, n in (r, s_):
            st_ = set(t for i, t in enumerate(groups) if t != "__ctrl__" and n[i].min() >= MIN_HALF)
            ok = st_ if ok is None else ok & st_
    targets = sorted(ok)
    gall = {g: j for j, g in enumerate(genes)}
    texp = np.array([max(ctrl_mean[("cell", "rest")].mean(0)[gall[t]], ctrl_mean[("cell", "stim")].mean(0)[gall[t]]) if t in gall else 0.0 for t in targets])   # both states (10-18 fix)
    inert = (texp < 0.01) & np.array([t in gall for t in targets])   # unmeasured targets are not "inert" (10-18 fix)
    gsets = {"all_expressed": np.where(np.minimum(ctrl_mean[("cell", "rest")].mean(0), ctrl_mean[("cell", "stim")].mean(0)) >= MIN_EXPR)[0]}
    gsets["top2000"] = gsets["all_expressed"][np.argsort(-cm[gsets["all_expressed"]])[:2000]]
    H, ncell = {}, {}
    for nm, (r, s_) in D.items():
        arr = []
        for st, (groups, _, S, n) in (("rest", r), ("stim", s_)):
            gi = {g: i for i, g in enumerate(groups)}; idx = np.array([gi[t] for t in targets])
            X = np.stack([np.asarray(S[idx, h])[:, col[st]] / n[idx, h][:, None] - ctrl_mean[(nm, st)][h] for h in range(2)])
            for k_, t in enumerate(targets):
                if t in gall:
                    X[:, k_, gall[t]] = 0
            arr.append(X); ncell[(nm, st)] = n[idx]
        H[nm] = arr
    bal = np.ones(len(targets), bool)
    for st in ("rest", "stim"):
        c = ncell[("guide", st)]; bal &= (np.minimum(c[:, 0], c[:, 1]) / np.maximum(c[:, 0], c[:, 1])) >= 1 / 1.5
    res = {"rules": __doc__, "n_targets": len(targets), "n_inert": int(inert.sum()), "n_balanced": int(bal.sum()),
           "median_cells_per_half": {"%s_%s" % k: np.median(v, 0).tolist() for k, v in ncell.items()}}
    for subset, m in (("all_common", np.ones(len(targets), bool)), ("guide_balanced", bal)):
        for gs_name, gk in gsets.items():
            for nm in H:
                R_, S_ = H[nm]
                XA, XB, YA, YB = R_[0][m][:, gk], R_[1][m][:, gk], S_[0][m][:, gk], S_[1][m][:, gk]
                inert_m = inert[m]
                d = {}
                for lab, A, B in (("rest", XA, XB), ("stim", YA, YB)):
                    a, b = A - A[inert_m].mean(0), B - B[inert_m].mean(0)
                    d["inert_frac_pos_" + lab] = float(((a * b).sum(1)[inert_m] > 0).mean())
                c0 = lambda M: M - M.mean(0)
                XA, XB, YA, YB = c0(XA), c0(XB), c0(YA), c0(YB)
                bA, gA, gB = (XA + YA) / 2, (XA - YA) / 2, (XB - YB) / 2
                den = float((gA * gB).sum()); d["gamma_share"] = den / (den + float((bA * (XB + YB) / 2).sum()))
                for k in (5, 10, 20, 50):
                    d["readout_k%d" % k] = expl(predict(bA, gB, k), gB) / den
                d["readout_k20_shuffled"] = expl(predict(bA, gB, 20, shuffle=True), gB) / den
                key = "%s|%s|%s" % (subset, gs_name, nm)
                res[key] = d
                print(key, json.dumps({k: round(v, 4) for k, v in d.items()}), flush=True)
    return res


def run_tahoe(npz, out):
    z = np.load(npz); X = {k: z[k].astype(np.float64) for k in ("plate6_0", "plate6_1", "plate6_pool", "plate14_0", "plate14_1", "plate14_pool")}
    def two_way(M):
        m = M.mean((0, 1)); a = M.mean(1) - m; b = M.mean(0) - m
        return M - a[:, None] - b - m, b
    res = {"rules": __doc__}
    for P_, O in (("plate14", "plate6"), ("plate6", "plate14")):
        sel = np.abs(X[O + "_pool"]).mean((0, 1)); top = np.argsort(-sel)[:2000]
        for gs_name, gk in (("all", np.arange(X[P_ + "_1"].shape[2])), ("top2000_other_plate", top)):
            g1, _ = two_way(X[P_ + "_1"][..., gk]); g0, b0 = two_way(X[P_ + "_0"][..., gk])
            den = float((g0 * g1).sum())
            srcs = {"same_well": b0, "other_plate": two_way(X[O + "_pool"][..., gk])[1]}
            d = {}
            for k in (5, 10, 20):
                for sname, b in srcs.items():
                    e = 0.0; es = 0.0
                    for seed in (3, 4):
                        rng = np.random.default_rng(seed); Pn = b.shape[0]; perm = rng.permutation(Pn); folds = (perm[: Pn // 2], perm[Pn // 2:])
                        F = np.zeros_like(g1); Fs = np.zeros_like(g1); bs_ = b[rng.permutation(Pn)]
                        for i in (0, 1):
                            tr, te = folds[i], folds[1 - i]
                            for src, FF in ((b, F), (bs_, Fs)):
                                V = np.linalg.svd(src[tr], full_matrices=False)[2][:k].T; Zm = src @ V
                                for c in range(g1.shape[0]):
                                    FF[c, te] = Zm[te] @ ridge(Zm[tr], g1[c, tr], pick_lam(Zm[tr], g1[c, tr], rng))
                        e += expl(F, g1) / 2; es += expl(Fs, g1) / 2
                    d["%s_k%d" % (sname, k)] = e / den; d["%s_k%d_shuffled" % (sname, k)] = es / den
            key = "target_%s|%s" % (P_, gs_name)
            res[key] = d
            print(key, json.dumps({k: round(v, 4) for k, v in d.items()}), flush=True)
    return res


def main():
    mode = sys.argv[1]
    if mode == "cd4":
        out = sys.argv[6]; os.makedirs(out, exist_ok=True); res = run_cd4(*sys.argv[2:7])
    else:
        out = sys.argv[3]; os.makedirs(out, exist_ok=True); res = run_tahoe(*sys.argv[2:4])
    json.dump(res, open(out + "/metrics_robust.json", "w"), indent=1)
    print("saved", flush=True)


if __name__ == "__main__":
    main()
