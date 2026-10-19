"""CD4 state axis, inert-target control done properly (A100). cd4_inert.py (and cd4_readout.py) did not remove the
template -- the shift shared by all perturbed groups relative to controls in each state and half -- so every target,
inert or not, carried it. Here, for the same data (Zhu 2026 donor 1, resting vs activated 48 h):
  variant A  original (no template removal), as cd4_readout.py
  variant B  template removed: each state x half, subtract the mean over all targets
  variant C  B plus projecting out the top-3 directions of inert targets' responses (fit on half A, per state)
  variant D  B, evaluated on expressed targets only (readout fit and scored within them)
For each variant: inert vs expressed per-target cross-half signal (median, fraction > 0); gamma signal share carried by
inert targets; gamma share of the specific signal; readout share at k = 20 (state = top-20 PCs of half-A beta, ridge,
2 folds, as cd4_readout.py) with its shuffled control. Inert = target control mean < 0.01 in both states.
Usage: cd4_inert2.py <out_Rest> <out_Stim> <out_dir>. Exploratory."""
import os, sys, json
import numpy as np
from cd4_readout import load, share, fit, MIN_HALF, MIN_EXPR

K = 20


def readout(betaA, gA, g, seed=11):
    P = len(g); perm = np.random.default_rng(seed).permutation(P); folds = (perm[: P // 2], perm[P // 2:])
    Fr = np.zeros_like(g); Fs = np.zeros_like(g)
    for f in (0, 1):
        tr, te = folds[f], folds[1 - f]; rng = np.random.default_rng(100 + f)
        V = np.linalg.svd(betaA[tr], full_matrices=False)[2][:K].T
        Z = betaA @ V; Fr[te] = Z[te] @ fit(Z[tr], g[tr], rng)
        sh = rng.permutation(len(tr)); Fs[te] = Z[te] @ fit(Z[tr][sh], g[tr], rng)
    return share(g, gA, Fr), share(g, gA, Fs)


def stats(XR, XS, inert, expd):
    betaA, betaB = (XR[0] + XS[0]) / 2, (XR[1] + XS[1]) / 2
    gA, gB = XR[0] - betaA, XR[1] - betaB
    gsig = (gA * gB).sum(1); bsig = (betaA * betaB).sum(1)
    d = {"gamma_share_of_specific": float(gsig.sum() / (gsig.sum() + bsig.sum())),
         "inert_gamma_signal_share": float(gsig[inert].sum() / gsig.sum()), "inert_target_share": float(inert.mean())}
    for name, m in (("inert", inert), ("expressed", expd)):
        sr, ss = (XR[0][m] * XR[1][m]).sum(1), (XS[0][m] * XS[1][m]).sum(1)
        d[name] = dict(rest_signal_median=float(np.median(sr)), rest_frac_pos=float((sr > 0).mean()),
                       stim_signal_median=float(np.median(ss)), stim_frac_pos=float((ss > 0).mean()))
    r, s_ = readout(betaA, gA, gB); d["readout_k20"] = r; d["readout_k20_shuffled"] = s_
    return d


def main():
    dr, ds, out = sys.argv[1:4]; os.makedirs(out, exist_ok=True)
    (gr, genes_r, Sr, nr), (gs, genes_s, Ss, ns) = load(dr), load(ds)
    cr, cs = gr.index("__ctrl__"), gs.index("__ctrl__")
    common = sorted(set(genes_r) & set(genes_s)); pr_, ps_ = {g: j for j, g in enumerate(genes_r)}, {g: j for j, g in enumerate(genes_s)}
    colr, cols = np.array([pr_[g] for g in common]), np.array([ps_[g] for g in common])
    ctr = np.stack([np.asarray(Sr[cr, h])[colr] / nr[cr, h] for h in range(2)]); cts = np.stack([np.asarray(Ss[cs, h])[cols] / ns[cs, h] for h in range(2)])
    cpos = {g: j for j, g in enumerate(common)}
    gk = np.where((ctr.mean(0) >= MIN_EXPR) & (cts.mean(0) >= MIN_EXPR))[0]; genes = [common[j] for j in gk]
    ir, is_ = {g: i for i, g in enumerate(gr)}, {g: i for i, g in enumerate(gs)}
    targets = sorted(t for t in gr if t != "__ctrl__" and t in is_ and nr[ir[t]].min() >= MIN_HALF and ns[is_[t]].min() >= MIN_HALF and t in cpos)
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
    texp = np.array([max(ctr[:, cpos[t]].mean(), cts[:, cpos[t]].mean()) for t in targets])   # inert = < 0.01 in BOTH states (10-18 fix; was min)
    inert, expd = texp < 0.01, texp >= 0.5
    res = {"rules": __doc__, "n_targets": P, "n_inert": int(inert.sum()), "n_expressed": int(expd.sum())}
    res["A_original"] = stats(XR, XS, inert, expd); print("A", json.dumps(res["A_original"]), flush=True)
    XRb = XR - XR.mean(1, keepdims=True); XSb = XS - XS.mean(1, keepdims=True)
    res["B_template_removed"] = stats(XRb, XSb, inert, expd); print("B", json.dumps(res["B_template_removed"]), flush=True)
    XRc, XSc = XRb.copy(), XSb.copy()
    for X in (XRc, XSc):
        V = np.linalg.svd(X[0][inert], full_matrices=False)[2][:3].T
        for h in range(2):
            X[h] -= (X[h] @ V) @ V.T
    res["C_template_and_inert_dirs_removed"] = stats(XRc, XSc, inert, expd); print("C", json.dumps(res["C_template_and_inert_dirs_removed"]), flush=True)
    m = expd
    res["D_template_removed_expressed_only"] = stats(XRb[:, m] - XRb[:, m].mean(1, keepdims=True), XSb[:, m] - XSb[:, m].mean(1, keepdims=True),
                                                    np.zeros(m.sum(), bool), np.ones(m.sum(), bool))
    print("D", json.dumps({k: v for k, v in res["D_template_removed_expressed_only"].items() if k != "inert"}), flush=True)
    json.dump(res, open(out + "/metrics_cd4_inert2.json", "w"), indent=1)
    print("saved", flush=True)


if __name__ == "__main__":
    main()
