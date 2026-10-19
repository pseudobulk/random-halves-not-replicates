"""Does the readout structure of Z1 (r) survive independent screens? (A100; Replogle 2022 raw single cells.)
Z1 (aj): in Tahoe the same target gamma is explained ~26% by beta of the same well but 3-7% by beta of a replicate
plate, i.e. most of the within-well readout reads fluctuations shared by beta and gamma of the same well. In (r) the
two halves are random cells of the same screen, so a K562 clone's deviation enters half-A beta (= (K + R) / 2) and
half-B gamma (= (K - R) / 2) alike. K562 has two independent screens (day 6 essential, day 8 genome-wide: different
transductions, library, time), RPE1 one (day 7 essential).
Data: K6 = K562_essential_raw_singlecell_01, K8 = K562_gwps_raw_singlecell_01 (part 1), R = rpe1_raw_singlecell_01.
Perturbations with >= 20 cells in all three, minus p53 / proliferation genes and TP53 (F26 rule); genes common to the
three files with control mean log1p CP10K >= 0.1 in all three; per-lane control correction and CP10K as explore_V
(Line). Responses per line: random cell halves A / B and all cells; target gene zeroed; centred over perturbations
(template removal). gamma = (K - R) / 2, beta = (K + R) / 2.
Target fixed: half-B gamma of the pair (K_T, R), T = K8 (and, reversed, K6). Score = explained energy of the target
on held-out perturbations / within-design gamma cross product (half A x half B), so all feature sources are in the
same units. Readout: state = top-k PCs of the feature over training perturbations, ridge (lambda from inner 2-fold CV),
2 folds over perturbations (seed 11); shuffled control permutes the training states. Feature sources:
  pair_same     half-A beta of the same pair (as r)             shares K_T clones and R clones with the target
  pair_indepK   (other K562 screen, all cells + R half A) / 2   shares R clones only
  K_same        K_T half A                                      shares K_T clones
  K_indep       other K562 screen, all cells                    shares nothing with the target
  K_indep_half  other K562 screen, half A (cell numbers ~ K_same)
  R_same        R half A                                        shares R clones
Also signal correlations between the screens: cross product of their halves / geometric mean of each screen's
own half x half cross product.
Expectation written before running: if the readout is a stable property of the lines, K_indep keeps roughly rho^2 or
more of K_same, rho = K6-K8 signal correlation (the other screen is an imperfect proxy only through day / library
differences). A Tahoe-like collapse (K_indep well below rho^2 x K_same) would mean that the same-screen readout
mostly reads clone-level coupling. The R side cannot be tested this way (one RPE1 screen).
Usage: replogle_indep_readout.py <K6.h5ad> <K8.h5ad> <R.h5ad> <out_dir>. Exploratory."""
import os, sys, json, time
import numpy as np
from explore_V_day6_day8 import obs_info, read_rows, Line, CTRL_NAMES, PRO, P53, load, MIN_CELLS

KS = (10, 20)
LAMS = (1e-3, 1e-2, 1e-1, 1, 10, 100)
MIN_EXPR = 0.1


def ridge(Z, Y, lam):
    k = Z.shape[1]
    return np.linalg.solve(Z.T @ Z + lam * np.trace(Z.T @ Z) / k * np.eye(k), Z.T @ Y)


def expl(F, g):
    return float((g ** 2).sum() - ((g - F) ** 2).sum())


def pick_lam(Z, Y, rng):
    P = len(Z); perm = rng.permutation(P); f = (perm[: P // 2], perm[P // 2:])
    best, bl = -np.inf, 1.0
    for lam in LAMS:
        e = sum(expl(Z[f[1 - i]] @ ridge(Z[f[i]], Y[f[i]], lam), Y[f[1 - i]]) for i in (0, 1))
        if e > best:
            best, bl = e, lam
    return bl


def readout(feat, g, k, seed=11, shuffle=False):
    P = len(g); rng = np.random.default_rng(seed); perm = rng.permutation(P); folds = (perm[: P // 2], perm[P // 2:])
    F = np.zeros_like(g)
    for i in (0, 1):
        tr, te = folds[i], folds[1 - i]
        V = np.linalg.svd(feat[tr], full_matrices=False)[2][:k].T; Z = feat @ V
        Ztr = Z[tr][rng.permutation(len(tr))] if shuffle else Z[tr]
        F[te] = Z[te] @ ridge(Ztr, g[tr], pick_lam(Ztr, g[tr], rng))
    return expl(F, g)


def effects(L, perts, seed):
    rng = np.random.default_rng(seed); mu = L.Y[L.ctrl].mean(0); G = L.Y.shape[1]
    H = np.zeros((2, len(perts), G)); Fu = np.zeros((len(perts), G))
    for k, p in enumerate(perts):
        rows = rng.permutation(np.where(L.codes == L.cats.index(p))[0]); h = len(rows) // 2
        H[0, k] = L.Y[rows[:h]].mean(0) - mu; H[1, k] = L.Y[rows[h:2 * h]].mean(0) - mu; Fu[k] = L.Y[rows].mean(0) - mu
        if p in L.g2j:
            H[:, k, L.g2j[p]] = 0; Fu[k, L.g2j[p]] = 0
    return H, Fu, mu


def cross(a, b):
    return float((a * b).sum())


def main():
    p6, p8, pr, out = sys.argv[1:5]; os.makedirs(out, exist_ok=True); t0 = time.time()
    paths = {"K6": p6, "K8": p8, "R": pr}
    meta = {}
    for nm, p in paths.items():
        codes, cats, gem, genes = obs_info(p); n = np.bincount(codes, minlength=len(cats))
        meta[nm] = dict(ok={c for c, k in zip(cats, n) if k >= MIN_CELLS}, genes=list(genes))
        print(nm, "cells", len(codes), "perts", len(cats), "genes", len(genes), flush=True)
    excl = set(P53) | set(PRO) | {"TP53"}
    perts = sorted(set.intersection(*[m["ok"] for m in meta.values()]) - set(CTRL_NAMES) - excl)
    common = sorted(set.intersection(*[set(m["genes"]) for m in meta.values()]))
    print("perturbations", len(perts), "common genes", len(common), flush=True)
    E = {}
    for seed, (nm, p) in enumerate(paths.items()):
        d = load(p, set(perts)); gidx = {g: j for j, g in enumerate(d["genes"])}
        X = read_rows(p, d["rows"], np.array([gidx[g] for g in common]))
        L = Line(X, d["codes"][d["rows"]], d["cats"], common, d["gem"][d["rows"]], d["ctrl"]); del X
        E[nm] = effects(L, perts, seed); del L
        print(nm, "effects done", round(time.time() - t0), "s", flush=True)
    gk = np.where(np.min([E[nm][2] for nm in E], 0) >= MIN_EXPR)[0]
    c = lambda X: X - X.mean(-2, keepdims=True)
    H = {nm: c(E[nm][0][..., gk]) for nm in E}; Fu = {nm: c(E[nm][1][:, gk]) for nm in E}
    np.savez(out + "/responses.npz", perts=np.array(perts), genes=np.array([common[j] for j in gk]),
             **{nm + "_H": H[nm].astype(np.float32) for nm in H}, **{nm + "_full": Fu[nm].astype(np.float32) for nm in Fu})
    res = {"rules": __doc__, "n_perturbations": len(perts), "n_genes": int(len(gk))}
    print("genes kept", len(gk), flush=True)
    sc = {}
    for a, b in (("K6", "K8"), ("K6", "R"), ("K8", "R")):
        wa, wb = cross(H[a][0], H[a][1]), cross(H[b][0], H[b][1])
        x = (cross(H[a][0], H[b][1]) + cross(H[a][1], H[b][0])) / 2
        sc[a + "_" + b] = dict(within_a=wa, within_b=wb, cross=x, signal_corr=float(x / np.sqrt(wa * wb)))
    res["signal_corr"] = sc; print("signal_corr", json.dumps(sc), flush=True)
    for T, O in (("K8", "K6"), ("K6", "K8")):
        gA, g = (H[T][0] - H["R"][0]) / 2, (H[T][1] - H["R"][1]) / 2
        bA, bB = (H[T][0] + H["R"][0]) / 2, (H[T][1] + H["R"][1]) / 2
        den, bs = cross(gA, g), cross(bA, bB)
        r = dict(gamma_share=den / (den + bs), gamma_cross=den, beta_cross=bs)
        feats = {"pair_same": bA, "pair_indepK": (Fu[O] + H["R"][0]) / 2, "K_same": H[T][0], "K_indep": Fu[O],
                 "K_indep_half": H[O][0], "R_same": H["R"][0]}
        for k in KS:
            for nm, f in feats.items():
                r["%s_k%d" % (nm, k)] = readout(f, g, k) / den
                r["%s_k%d_shuffled" % (nm, k)] = readout(f, g, k, shuffle=True) / den
        res["target_" + T] = r
        print("target", T, json.dumps(r), round(time.time() - t0), "s", flush=True)
    json.dump(res, open(out + "/metrics_replogle_indep_readout.json", "w"), indent=1)
    print("saved", flush=True)


if __name__ == "__main__":
    main()
