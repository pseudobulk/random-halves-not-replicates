"""Minimal model of non-independent replicate units in pooled perturbation screens (Z1 ap). Simulation, no data.
Cell i of perturbation p in context c, from unit k (clone or well):
    y_i = mu_pc (1 + g_k) + a_k + e_i
  additive a_k ~ N(0, SIG_A^2 I): the unit's own state, present whatever the perturbation does (inert targets too)
  gain     g_k ~ N(0, SIG_G^2):   the unit scales the perturbation's response (e.g. drug exposure of a well,
                                  responsiveness / activation of a clone); zero for inert targets
  e_i ~ N(0, SIG_E^2 I). K units per perturbation and transduction, Dirichlet(ALPHA) weights (jackpots).
  truth: beta_p = s_p (z_p W); gamma_p = s_p (0.5 (z_p Rm) W2 + 0.3 d_p); X = beta + gamma, Y = beta - gamma;
         20% inert (s_p = 0). The true readout share (state from beta -> gamma) is about 0.7.
  shared units (one well holding all lines; one clone in two states): the same units enter X and Y, the additive
  part scaled by RHO in Y, the gain identical. separate (lines screened apart): X and Y have their own units.
Designs with the same cells per half: random = one transduction, cells split at random (both halves share its units);
independent = each half from its own transduction.
Expectations (one context, gene-summed cross products, w = unit weights, K_eff = 1 / sum w^2):
  random:      E[x_A . x_B] = |mu|^2 (1 + SIG_G^2 / K_eff) + G SIG_A^2 / K_eff
  independent: E[x_A . x_B] = |mu|^2
  -> the additive term makes inert targets reproducible and inflates ceilings; neither unit term shrinks with more
     cells, so the Spearman-Brown doubling is wrong for them; the gain inflates beta and gamma by the same factor, so
     the gamma share is unchanged; with shared units the gain scales beta_A and gamma_B together, so a readout from
     beta to gamma keeps working within a transduction and loses accuracy across transductions.
Attenuation block (10-18 audit; 1.5-2.5 added 10-19): true readout ~0.7, shared additive units with SIG_A 0.2 / 0.6 / 1.2 / 1.5 / 2.0 / 2.5 (iso or along
W); reports inert vs expressed cross energy per target, so the unit-to-signal ratio can be matched to CD4 (inert
random-half cross about 4x the expressed guide-half cross), and the readout seen by random vs independent designs.
Fabrication block: truth without any readout (gamma unrelated to the state), weak effects (scale 0.1), additive units
along the shared response directions W, shared vs separate: can non-independent halves create a readout from nothing?
Outputs per scenario x design (seed-averaged): inert frac_pos, gamma share, readout share (k = 10, 2-fold, ridge),
pearson_delta ceiling (Spearman-Brown) for all and inert perturbations in context X.
Usage: clone_model_sim.py <out_json>. Exploratory."""
import sys, json
import numpy as np

G, P, R, KREAD, NCELL, K, ALPHA, SIG_E, RHO = 1000, 800, 10, 10, 120, 6, 0.5, 1.0, 0.5
SIG_A, SIG_G = 0.2, 0.6
SEEDS = (1, 2, 3)


def truth(rng, scale=0.3, read_w=0.5):
    W, W2 = rng.normal(size=(R, G)), rng.normal(size=(R, G))
    Rm = rng.normal(size=(R, R)) / np.sqrt(R)
    z = rng.normal(size=(P, R)) / np.sqrt(R)
    s = np.exp(rng.normal(np.log(scale), 0.5, P)); s[: P // 5] = 0
    beta = s[:, None] * (z @ W)
    gamma = s[:, None] * (read_w * (z @ Rm) @ W2 + 0.3 * rng.normal(size=(P, G)))
    return beta, gamma, s == 0, W


def transduction(rng, add, gain, W=None, sig_a=None):
    """add: False, True (isotropic) or "low" (along the shared response directions W, same per-gene variance)"""
    sa = SIG_A if sig_a is None else sig_a
    if add == "low":
        A = sa * (rng.normal(size=(P, K, R)) / np.sqrt(R)) @ W
    else:
        A = sa * rng.normal(size=(P, K, G)) if add else np.zeros((P, K, G))
    return dict(w=rng.dirichlet(np.full(K, ALPHA), size=P), A=A,
                g=SIG_G * rng.normal(size=(P, K)) if gain else np.zeros((P, K)))


def half(rng, mu, T, n, a_scale=1.0):
    f = np.stack([rng.multinomial(n, wi) for wi in T["w"]]) / n
    return mu * (1 + (f * T["g"]).sum(1))[:, None] + a_scale * np.einsum("pk,pkg->pg", f, T["A"]) + rng.normal(size=mu.shape) * SIG_E / np.sqrt(n)


def ridge(Z, Y, lam=0.1):
    k = Z.shape[1]
    return np.linalg.solve(Z.T @ Z + lam * np.trace(Z.T @ Z) / k * np.eye(k), Z.T @ Y)


LAMS = (1e-2, 1e-1, 1, 10, 100, 1000)


def pick_lam(Z, Y, rng):
    n = len(Z); perm = rng.permutation(n); f = (perm[: n // 2], perm[n // 2:]); best, bl = -np.inf, 1.0
    for lam in LAMS:
        e = sum(((Y[f[1 - i]] ** 2).sum() - ((Y[f[1 - i]] - Z[f[1 - i]] @ ridge(Z[f[i]], Y[f[i]], lam)) ** 2).sum()) for i in (0, 1))
        if e > best:
            best, bl = e, lam
    return bl


def readout(bA, gA, gB, rng):
    """as the data analyses: top-k PCs of beta_A on training perturbations, ridge with lambda from inner 2-fold CV"""
    perm = rng.permutation(P); folds = (perm[: P // 2], perm[P // 2:]); F = np.zeros_like(gB)
    for i in (0, 1):
        tr, te = folds[i], folds[1 - i]
        V = np.linalg.svd(bA[tr], full_matrices=False)[2][:KREAD].T; Z = bA @ V
        F[te] = Z[te] @ ridge(Z[tr], gB[tr], pick_lam(Z[tr], gB[tr], rng))
    return float(((gB ** 2).sum() - ((gB - F) ** 2).sum()) / (gA * gB).sum())


def metrics(XA, XB, YA, YB, inert, rng):
    c = lambda M: M - M.mean(0)
    XA, XB, YA, YB = c(XA), c(XB), c(YA), c(YB)
    bA, bB, gA, gB = (XA + YA) / 2, (XB + YB) / 2, (XA - YA) / 2, (XB - YB) / 2
    gs, bs = (gA * gB).sum(), (bA * bB).sum()
    ref = lambda M: M - M[inert].mean(0)
    xa, xb = ref(XA), ref(XB)
    sig = (xa * xb).sum(1)
    a, b = xa - xa.mean(1, keepdims=True), xb - xb.mean(1, keepdims=True)
    pr = (a * b).sum(1) / np.sqrt((a ** 2).sum(1) * (b ** 2).sum(1))
    sb = lambda r: 2 * r / (1 + r) if r > 0 else float("nan")
    return dict(inert_frac_pos=float((sig[inert] > 0).mean()), gamma_share=float(gs / (gs + bs)),
                readout_share=readout(bA, gA, gB, rng), ceiling_pearson_all=sb(pr.mean()),
                ceiling_pearson_inert=sb(pr[inert].mean()))


def main():
    out = sys.argv[1]
    res = {"rules": __doc__, "params": dict(G=G, P=P, R=R, KREAD=KREAD, NCELL=NCELL, K=K, ALPHA=ALPHA, SIG_E=SIG_E,
                                            RHO=RHO, SIG_A=SIG_A, SIG_G=SIG_G, seeds=SEEDS)}
    scen = {"none": (False, False, False), "additive_separate": (True, False, False), "additive_shared": (True, False, True),
            "gain_separate": (False, True, False), "gain_shared": (False, True, True), "both_shared": (True, True, True)}
    acc = {}
    n = NCELL // 2
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        beta, gamma, inert, W = truth(rng)
        X, Y = beta + gamma, beta - gamma
        tr = metrics(X, X, Y, Y, inert, rng)
        acc.setdefault("truth", []).append(dict(gamma_share=tr["gamma_share"], readout_share=tr["readout_share"]))
        for nm, (add, gain, shared) in scen.items():
            for design in ("random", "independent"):
                if design == "random":
                    TX = transduction(rng, add, gain); TY = TX if shared else transduction(rng, add, gain)
                    TXs, TYs = (TX, TX), (TY, TY)
                else:
                    TXs, TYs = [], []
                    for _ in range(2):
                        TX = transduction(rng, add, gain); TXs.append(TX); TYs.append(TX if shared else transduction(rng, add, gain))
                ys = RHO if shared else 1.0
                XA, XB = half(rng, X, TXs[0], n), half(rng, X, TXs[1], n)
                YA, YB = half(rng, Y, TYs[0], n, ys), half(rng, Y, TYs[1], n, ys)
                acc.setdefault("%s|%s" % (nm, design), []).append(metrics(XA, XB, YA, YB, inert, rng))
        # fabrication block: no true readout (gamma unrelated to the state), weak effects, additive units along the
        # response directions (as an activation-level or exposure axis would be)
        beta, gamma, inert, W = truth(rng, scale=0.1, read_w=0.0)
        X, Y = beta + gamma, beta - gamma
        tr = metrics(X, X, Y, Y, inert, rng)
        acc.setdefault("fab_truth", []).append(dict(gamma_share=tr["gamma_share"], readout_share=tr["readout_share"]))
        for nm, shared in (("fab_low_shared", True), ("fab_low_separate", False), ("fab_none", None)):
            for design in ("random", "independent"):
                add = False if shared is None else "low"
                if design == "random":
                    TX = transduction(rng, add, False, W); TY = TX if shared else transduction(rng, add, False, W)
                    TXs, TYs = (TX, TX), (TY, TY)
                else:
                    TXs, TYs = [], []
                    for _ in range(2):
                        TX = transduction(rng, add, False, W); TXs.append(TX); TYs.append(TX if shared else transduction(rng, add, False, W))
                ys = RHO if shared else 1.0
                XA, XB = half(rng, X, TXs[0], n), half(rng, X, TXs[1], n)
                YA, YB = half(rng, Y, TYs[0], n, ys), half(rng, Y, TYs[1], n, ys)
                acc.setdefault("%s|%s" % (nm, design), []).append(metrics(XA, XB, YA, YB, inert, rng))
        # attenuation block (10-18 audit): a REAL readout (truth 0.7) with unit deviations as large as in CD4
        # (inert cross energy about 4x the independent-unit signal of expressed targets): does the independent design
        # still see the readout, or do noisy features (unit deviations in beta_A) wipe it out?
        beta, gamma, inert, W = truth(rng)
        X, Y = beta + gamma, beta - gamma
        for sa in (0.2, 0.6, 1.2, 1.5, 2.0, 2.5):   # 1.5-2.5 added 10-19 to reach the CD4 unit-to-signal ratio (4.0, 7.5)
            for kind in (True, "low"):
                for design in ("random", "independent"):
                    if design == "random":
                        TX = transduction(rng, kind, False, W, sa); TXs, TYs = (TX, TX), (TX, TX)
                    else:
                        TXs = [transduction(rng, kind, False, W, sa) for _ in range(2)]; TYs = TXs
                    XA, XB = half(rng, X, TXs[0], n), half(rng, X, TXs[1], n)
                    YA, YB = half(rng, Y, TYs[0], n, RHO), half(rng, Y, TYs[1], n, RHO)
                    m = metrics(XA, XB, YA, YB, inert, rng)
                    c = lambda M: M - M.mean(0)
                    ia = c(XA) - c(XA)[inert].mean(0); ib = c(XB) - c(XB)[inert].mean(0)
                    m["inert_cross_per_target"] = float((ia * ib).sum(1)[inert].mean())
                    m["expressed_cross_per_target"] = float((ia * ib).sum(1)[~inert].mean())
                    acc.setdefault("atten_sa%.1f_%s|%s" % (sa, "low" if kind == "low" else "iso", design), []).append(m)
        print("seed", seed, "done", flush=True)
    for key, lst in acc.items():
        res[key] = {m: float(np.nanmean([d[m] for d in lst])) for m in lst[0]}
        print(key, json.dumps({m: round(v, 3) for m, v in res[key].items()}), flush=True)
    json.dump(res, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
