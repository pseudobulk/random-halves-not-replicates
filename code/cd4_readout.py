"""Readout structure on the STATE axis: CD4 T resting vs activated 48 h (Zhu 2026 donor 1; A100). Z1 r / s showed the
cross-line gamma is a line-specific linear readout of a shared low-dim state (CRISPR pairs 50-62%, Tahoe drugs 73%);
Z1 u showed no such structure between iPSC clones. The deepened frame predicts structure here: activation rewires the
alarm -> effector arms (mTORC1 / MYC growth, apoptotic priming, NF-kB / IFN), so the resting-vs-activated gamma should
be a readout of the shared state, loading on those arms.
Data: cd4_sc_extract.py aggregates (out_Rest, out_Stim: per target and random cell half, sums of log1p CP10K; controls
"__ctrl__"). Targets with >= MIN_HALF cells per half in both states; genes with control mean >= MIN_EXPR in both.
Effects x_s = sum / n - ctrl (per half). beta = (x_rest + x_stim) / 2; gamma_rest = x_rest - beta = -gamma_stim. Target
= half-B gamma; state from half-A beta (top k PCs of the training targets; or named sensor scores); ridge readout with
inner-CV lambda; 2 folds over targets (seed 11); share of the denoised gamma signal (explore_P.share). Controls: state
rows shuffled; low-rank ceiling from gamma itself (needs gamma measured); source-only state (resting only; activated
only). Readout content: top genes and Hallmark enrichment; PC alignment of beta with the sensor sets.
Usage: cd4_readout.py <out_Rest> <out_Stim> <hallmark txt> <out_dir>. Exploratory."""
import os, sys, json
import numpy as np
from scipy.stats import hypergeom
from tahoe_arms import ISR, HSF1, NRF2, P53, HALL

MIN_HALF, MIN_EXPR = 10, 0.1
KS = (5, 10, 20, 50)
LAMS = (1e-3, 1e-2, 1e-1, 1, 10, 100)


def load(d):
    m = json.load(open(d + "/meta.json")); S = np.load(d + "/sum_half.npy", mmap_mode="r")
    n = np.load(d + "/ctrl_moments.npz")["n_half"]
    return m["groups"], m["genes"], S, n


def share(g, gA, F):
    sig = (g ** 2).sum() - ((gA - g) ** 2).sum() / 2
    return float(((g ** 2).sum() - ((g - F) ** 2).sum()) / sig)


def ridge(Z, Y, lam):
    k = Z.shape[1]
    return np.linalg.solve(Z.T @ Z + lam * np.trace(Z.T @ Z) / k * np.eye(k), Z.T @ Y)


def fit(Z, Y, rng):
    n = len(Z); perm = rng.permutation(n); a, b = perm[: n // 2], perm[n // 2:]
    best = None
    for lam in LAMS:
        e = sum(((Y[v] - Z[v] @ ridge(Z[t], Y[t], lam)) ** 2).sum() for t, v in ((a, b), (b, a)))
        if best is None or e < best[0]:
            best = (e, lam)
    return ridge(Z, Y, best[1])


def main():
    dr, ds, hpath, out = sys.argv[1:5]; os.makedirs(out, exist_ok=True)
    (gr, genes_r, Sr, nr), (gs, genes_s, Ss, ns) = load(dr), load(ds)
    cr, cs = gr.index("__ctrl__"), gs.index("__ctrl__")
    common = sorted(set(genes_r) & set(genes_s)); pr_, ps_ = {g: j for j, g in enumerate(genes_r)}, {g: j for j, g in enumerate(genes_s)}
    colr, cols = np.array([pr_[g] for g in common]), np.array([ps_[g] for g in common])
    ctr = np.stack([np.asarray(Sr[cr, h])[colr] / nr[cr, h] for h in range(2)]); cts = np.stack([np.asarray(Ss[cs, h])[cols] / ns[cs, h] for h in range(2)])
    gk = np.where((ctr.mean(0) >= MIN_EXPR) & (cts.mean(0) >= MIN_EXPR))[0]; genes = [common[j] for j in gk]
    ir, is_ = {g: i for i, g in enumerate(gr)}, {g: i for i, g in enumerate(gs)}
    targets = sorted(t for t in gr if t != "__ctrl__" and t in is_ and nr[ir[t]].min() >= MIN_HALF and ns[is_[t]].min() >= MIN_HALF)
    P, G = len(targets), len(gk); print("targets", P, "genes", G, flush=True)
    XR = np.zeros((2, P, G)); XS = np.zeros((2, P, G))
    for k, t in enumerate(targets):
        for h in range(2):
            XR[h, k] = np.asarray(Sr[ir[t], h])[colr[gk]] / nr[ir[t], h] - ctr[h, gk]
            XS[h, k] = np.asarray(Ss[is_[t], h])[cols[gk]] / ns[is_[t], h] - cts[h, gk]
    gpos = {g: j for j, g in enumerate(genes)}
    for k, t in enumerate(targets):
        if t in gpos:
            XR[:, k, gpos[t]] = 0; XS[:, k, gpos[t]] = 0
    betaA, betaB = (XR[0] + XS[0]) / 2, (XR[1] + XS[1]) / 2
    gA, g = XR[0] - betaA, XR[1] - betaB                      # gamma of the resting state (activated = -gamma)
    raw_share = float(((gA * g).sum()) / ((gA * g).sum() + (betaA * betaB).sum()))
    res = {"rules": __doc__, "n_targets": P, "n_genes": G, "gamma_share_of_specific_signal": raw_share,
           "signal_fraction_of_gamma_energy": float(((g ** 2).sum() - ((gA - g) ** 2).sum() / 2) / (g ** 2).sum())}
    perm = np.random.default_rng(11).permutation(P); folds = (perm[: P // 2], perm[P // 2:])
    def run_state(src, k):
        Fr = np.zeros_like(g); Fs = np.zeros_like(g); Fl = np.zeros_like(g)
        for f in (0, 1):
            tr, te = folds[f], folds[1 - f]; rng = np.random.default_rng(100 + f)
            V = np.linalg.svd(src[tr], full_matrices=False)[2][:k].T
            Z = src @ V; M = fit(Z[tr], g[tr], rng); Fr[te] = Z[te] @ M
            sh = rng.permutation(len(tr)); Ms = fit(Z[tr][sh], g[tr], rng); Fs[te] = Z[te] @ Ms
            W = np.linalg.svd((gA[tr] + g[tr]) / 2, full_matrices=False)[2][:k].T; Fl[te] = (gA[te] @ W) @ W.T
        return share(g, gA, Fr), share(g, gA, Fs), share(g, gA, Fl), Fr
    for k in KS:
        r, s_, l, _ = run_state(betaA, k)
        res["k%d" % k] = dict(readout=r, shuffled=s_, lowrank_needs_gamma=l,
                              src_rest=run_state(XR[0], k)[0], src_stim=run_state(XS[0], k)[0])
        print("k", k, {a: round(b, 4) for a, b in res["k%d" % k].items()}, flush=True)
    # named sensor state
    H = {}
    for line in open(hpath):
        parts = [x for x in line.rstrip("\n").split("\t") if x]
        if len(parts) > 1:
            H[parts[0]] = set(parts[1:])
    Sets = {h: H[h] for h in HALL if h in H}; Sets.update(ISR_ATF4=set(ISR), HSF1_heat_shock=set(HSF1), NRF2_oxidative=set(NRF2), p53_core=set(P53))
    Sets = {h: np.array(sorted(gpos[x] for x in s if x in gpos)) for h, s in Sets.items()}; Sets = {h: ix for h, ix in Sets.items() if len(ix) >= 10}
    names = list(Sets)
    Zn = np.column_stack([betaA[:, Sets[h]].mean(1) for h in names]); Zn = (Zn - Zn.mean(0)) / np.maximum(Zn.std(0), 1e-12)
    Fr = np.zeros_like(g)
    for f in (0, 1):
        tr, te = folds[f], folds[1 - f]; Fr[te] = Zn[te] @ fit(Zn[tr], g[tr], np.random.default_rng(100 + f))
    res["named_state"] = dict(k=len(names), readout=share(g, gA, Fr))
    M = fit(Zn, g, np.random.default_rng(0)); W = np.array([[M[i, Sets[h2]].mean() for h2 in names] for i in range(len(names))]) / np.abs(g).mean()
    flat = sorted(((abs(W[i, j]), names[i], names[j], round(float(W[i, j]), 3)) for i in range(len(names)) for j in range(len(names))), reverse=True)[:12]
    res["named_wiring_top"] = [(a, b, c) for _, a, b, c in flat]
    # PC alignment of beta
    U, sv, Vt = np.linalg.svd(betaA, full_matrices=False); absb = np.abs(betaA).mean(0)
    dec = np.digitize(absb, np.quantile(absb, np.linspace(0.1, 0.9, 9))); bydec = {q: np.where(dec == q)[0] for q in range(10)}
    rng = np.random.default_rng(3); pcs = []
    for j in range(12):
        v = Vt[j]; zs = {}
        for h in names:
            ix = Sets[h]; m = v[ix].mean()
            nulls = np.array([v[np.concatenate([rng.choice(bydec[q], (dec[ix] == q).sum(), replace=False) for q in range(10) if (dec[ix] == q).sum()])].mean() for _ in range(50)])
            zs[h] = float((m - nulls.mean()) / max(nulls.std(), 1e-12))
        best = max(zs, key=lambda n_: abs(zs[n_])); pcs.append((j, round(float(sv[j] ** 2 / (sv ** 2).sum()), 3), best, round(zs[best], 1)))
    res["pc_alignment"] = pcs
    # readout content (k = 20, all targets)
    V = np.linalg.svd(betaA, full_matrices=False)[2][:20].T; M = fit(betaA @ V, g, np.random.default_rng(0))
    w = np.sqrt((M ** 2).sum(0)); top = np.argsort(-w)[:150]; sel = {genes[j] for j in top}; enr = []
    for h, hs in H.items():
        hs = hs & set(genes); x = len(hs & sel)
        if x >= 3:
            enr.append((h, x, len(hs), float(hypergeom.sf(x - 1, G, len(hs), len(top)))))
    res["readout_top_genes"] = [genes[j] for j in top[:30]]; res["readout_hallmark"] = sorted(enr, key=lambda t: t[3])[:8]
    print("named", res["named_state"], res["named_wiring_top"][:8], flush=True)
    print("pcs", pcs, flush=True)
    print("top genes", res["readout_top_genes"][:20], res["readout_hallmark"][:5], flush=True)
    json.dump(res, open(out + "/metrics_cd4_readout.json", "w"), indent=1, default=float)
    print("saved", flush=True)


if __name__ == "__main__":
    main()
