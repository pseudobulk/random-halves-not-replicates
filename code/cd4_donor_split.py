"""CD4 state axis with DONOR halves (A100; Zhu 2026): donors {1,2} vs {3,4} are independent biological replicates and
independent transductions, so clone- and donor-specific states cannot appear reproducible. Guides pooled within each
half (both guide ranks), which keeps more signal than the guide split. Otherwise as cd4_inert_guide.py (whose docstring
follows).
Is the inert-target component (Z1 ag) a clone effect of random cell halves? (A100; CD4, Zhu 2026, donors pooled.)
In cd4_inert2.py the two halves were random cells of the same target, so cells descended from the same few transduced
clones fall in both halves and any clone-specific state looks "reproducible". Here the halves are the target's guides
(even / odd guide rank; independent transductions, so different clones), as in cd4_m3.py. Same comparisons as
cd4_inert2.py variants A (no template) and B (template removed), for Rest vs Stim48hr: per-target cross-half signal
for inert (target control log1p CP10K < 0.01 in both states) vs expressed (>= 0.5) targets; gamma signal carried by
inert targets; gamma share; readout share at k = 20 with shuffled control.
If the inert reproducibility collapses under the guide split, the component is a clone (or guide) effect of the
cell-split design, and guide-split denoising is the appropriate default.
Usage: cd4_inert_guide.py <pseudobulk_merged.h5ad> <out_dir>. Exploratory."""
import os, sys, json
import numpy as np
import h5py
from scipy.sparse import csr_matrix
from cd4_m3 import cat, MIN_CELLS, MIN_EXPR, P53, PRO
from cd4_inert2 import stats, fit, K


def main():
    path, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    with h5py.File(path) as f:
        o = f["obs"]
        donors, dc = cat(o, "donor_id"); conds, cc = cat(o, "culture_condition"); runs, _ = cat(o, "10xrun_id")
        guides, gc = cat(o, "guide_id")
        genes_t, tc = cat(o, "perturbed_gene_name"); gtypes, gtc = cat(o, "guide_type")
        ncell = o["n_cells"][:]
        gname = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var"]["gene_name"][:]])
        X = csr_matrix((f["X"]["data"][:], f["X"]["indices"][:], f["X"]["indptr"][:]), shape=tuple(f["X"].attrs["shape"]))
    print("rows", X.shape, "donors", list(donors), "conditions", list(conds), "runs", list(runs), "guide types", list(gtypes), flush=True)
    ctrl_type = [i for i, g in enumerate(gtypes) if "non" in g.lower()][0]
    isc = gtc == ctrl_type
    D, S, R = len(donors), len(conds), 2
    T = len(genes_t)
    tgt = np.where(isc, T, tc)                                   # control rows -> group T
    # half = rank of the guide within its target (controls: within the non-targeting set), even / odd
    up, pinv = np.unique(tgt * (gc.max() + 1) + gc, return_inverse=True)   # sorted by target, then guide
    ut = up // (gc.max() + 1)
    first = np.searchsorted(ut, ut, side="left")
    rc = ((np.arange(len(up)) - first) % 2)[pinv]
    key = ((tgt * D + dc) * S + cc) * R + rc
    uk, inv = np.unique(key, return_inverse=True)
    Gm = csr_matrix((np.ones(len(key)), (inv, np.arange(len(key)))), shape=(len(uk), len(key)))
    agg = (Gm @ X).tocsr()                                       # groups x genes, summed counts
    cells = np.bincount(inv, weights=ncell, minlength=len(uk))
    lib = np.asarray(agg.sum(1)).ravel()
    print("groups", len(uk), flush=True)
    kt = uk // (D * S * R); kd = (uk // (S * R)) % D; ks = (uk // R) % S; kr = uk % R
    pos = {k: i for i, k in enumerate(uk)}
    def row(t, d, s, r):
        return pos.get(((t * D + d) * S + s) * R + r, -1)
    # gene filter on controls
    ctrl_mean = []
    for d in range(D):
        for s in range(S):
            rr = [row(T, d, s, r) for r in range(R) if row(T, d, s, r) >= 0]
            v = np.asarray(agg[rr].sum(0)).ravel(); ctrl_mean.append(np.log1p(v / v.sum() * 1e4))
    gk = np.where(np.min(ctrl_mean, 0) >= MIN_EXPR)[0]
    gpos = {g: k for k, g in enumerate(gname[gk])}
    excl = set(P53) | set(PRO) | {"TP53"}
    print("genes", len(gk), flush=True)
    def effect(t, ctx, r):
        """ctx: list of (donor, state) pooled; returns log1p CP10K diff vs pooled controls of the same contexts, run r."""
        rt = [row(t, d, s, r) for d, s in ctx]; rc_ = [row(T, d, s, r) for d, s in ctx]
        if min(rt) < 0 or min(rc_) < 0 or cells[rt].sum() < MIN_CELLS:
            return None
        a = np.asarray(agg[rt].sum(0)).ravel(); c = np.asarray(agg[rc_].sum(0)).ravel()
        e = (np.log1p(a / a.sum() * 1e4) - np.log1p(c / c.sum() * 1e4))[gk]
        if genes_t[t] in gpos:
            e[gpos[genes_t[t]]] = 0
        return e
    S_rest, S_stim = list(conds).index("Rest"), list(conds).index("Stim48hr")
    def effect2(t, ctx):
        rt = [row(t, d, s, r) for d, s in ctx for r in range(R) if row(t, d, s, r) >= 0]
        rc_ = [row(T, d, s, r) for d, s in ctx for r in range(R) if row(T, d, s, r) >= 0]
        if not rt or not rc_ or cells[rt].sum() < 2 * MIN_CELLS:
            return None
        a = np.asarray(agg[rt].sum(0)).ravel(); c = np.asarray(agg[rc_].sum(0)).ravel()
        e = (np.log1p(a / a.sum() * 1e4) - np.log1p(c / c.sum() * 1e4))[gk]
        if genes_t[t] in gpos:
            e[gpos[genes_t[t]]] = 0
        return e
    ctrl_full = {}
    for s in (S_rest, S_stim):
        rr = [row(T, d, s, r) for d in range(D) for r in range(R) if row(T, d, s, r) >= 0]
        v = np.asarray(agg[rr].sum(0)).ravel(); ctrl_full[s] = np.log1p(v / v.sum() * 1e4)
    gidx = {g: j for j, g in enumerate(gname)}
    keepT, XR, XS = [], [], []
    for t in range(T):
        halves = ([0, 1], [2, 3])
        eR = [effect2(t, [(d, S_rest) for d in hv]) for hv in halves]
        eS = [effect2(t, [(d, S_stim) for d in hv]) for hv in halves]
        if any(x is None for x in eR + eS) or genes_t[t] not in gidx:
            continue
        keepT.append(t); XR.append(eR); XS.append(eS)
    XR = np.array(XR).transpose(1, 0, 2); XS = np.array(XS).transpose(1, 0, 2)
    texp = np.array([max(ctrl_full[S_rest][gidx[genes_t[t]]], ctrl_full[S_stim][gidx[genes_t[t]]]) for t in keepT])   # inert in BOTH states (10-18 audit fix)
    inert, expd = texp < 0.01, texp >= 0.5
    res = {"rules": __doc__, "n_targets": len(keepT), "n_inert": int(inert.sum()), "n_expressed": int(expd.sum()), "n_genes": int(len(gk))}
    print("targets", len(keepT), "inert", int(inert.sum()), "expressed", int(expd.sum()), flush=True)
    res["A_original_donor_split"] = stats(XR, XS, inert, expd); print("A", json.dumps(res["A_original_donor_split"]), flush=True)
    XRb = XR - XR.mean(1, keepdims=True); XSb = XS - XS.mean(1, keepdims=True)
    res["B_template_removed_donor_split"] = stats(XRb, XSb, inert, expd); print("B", json.dumps(res["B_template_removed_donor_split"]), flush=True)
    # C (10-18 Methods audit): referenced to the mean of the inert targets (null baseline), the preprint's main rule
    XRc = XR - XR[:, inert].mean(1, keepdims=True); XSc = XS - XS[:, inert].mean(1, keepdims=True)
    res["C_null_baseline_donor_split"] = stats(XRc, XSc, inert, expd); print("C", json.dumps(res["C_null_baseline_donor_split"]), flush=True)
    # Responsive subsets. Selecting by the cross-half energy that is then evaluated biases the subset toward positive
    # cross-half signal (inert targets inside it were 98.8% positive), so two selection rules that do not reuse it:
    #   kd: knockdown of the target gene itself (pooled over donors and both states; the target gene is zeroed in the
    #       response profiles, so its counting noise is independent of theirs); expressed targets only
    #   eA: energy of donor half A alone. Pure noise in B cannot become reproducible by selecting on A, so inert targets
    #       stay a valid null; with signal present, selection on A inflates the s.n_A term, which makes the share with
    #       the cross-product denominator conservative and the default share optimistic, so both are reported
    rcall = [row(T, d, s, r) for d in range(D) for s in (S_rest, S_stim) for r in range(R) if row(T, d, s, r) >= 0]
    cv = np.asarray(agg[rcall].sum(0)).ravel(); cv = np.log1p(cv / cv.sum() * 1e4)
    def tlev(t):
        rt = [row(t, d, s, r) for d in range(D) for s in (S_rest, S_stim) for r in range(R) if row(t, d, s, r) >= 0]
        a = np.asarray(agg[rt].sum(0)).ravel(); return np.log1p(a[gidx[genes_t[t]]] / a.sum() * 1e4)
    kd = np.array([cv[gidx[genes_t[t]]] - tlev(t) for t in keepT])
    eA = (XRb[0] ** 2).sum(1) + (XSb[0] ** 2).sum(1)
    def readout_x(XR_, XS_, seed=11):
        betaA = (XR_[0] + XS_[0]) / 2; gA = XR_[0] - betaA; g = XR_[1] - (XR_[1] + XS_[1]) / 2
        P = len(g); perm = np.random.default_rng(seed).permutation(P); folds = (perm[: P // 2], perm[P // 2:])
        Fr = np.zeros_like(g); Fs = np.zeros_like(g)
        for f in (0, 1):
            tr, te = folds[f], folds[1 - f]; rng = np.random.default_rng(100 + f)
            V = np.linalg.svd(betaA[tr], full_matrices=False)[2][:K].T
            Z = betaA @ V; Fr[te] = Z[te] @ fit(Z[tr], g[tr], rng)
            sh = rng.permutation(len(tr)); Fs[te] = Z[te] @ fit(Z[tr][sh], g[tr], rng)
        gs = (gA * g).sum(); ex = lambda F: float(((g ** 2).sum() - ((g - F) ** 2).sum()) / gs)
        return ex(Fr), ex(Fs)
    res["kd_summary"] = dict(expressed_kd_median=float(np.median(kd[expd])), expressed_kd_q90=float(np.quantile(kd[expd], 0.9)))
    print("kd", json.dumps(res["kd_summary"]), flush=True)
    for rule, score, pool in (("kd", kd, expd), ("eA", eA, np.ones(len(keepT), bool))):
        for N in (1000, 300):
            cand = np.where(pool)[0]; top = cand[np.argsort(-score[cand])[:N]]
            sub = lambda X: X[:, top] - X[:, top].mean(1, keepdims=True)
            d = stats(sub(XRb), sub(XSb), inert[top], expd[top]); d["n"] = int(len(top)); d["n_inert"] = int(inert[top].sum())
            d["readout_k20_xden"], d["readout_k20_shuffled_xden"] = readout_x(sub(XRb), sub(XSb))
            res["C_%s_top%d" % (rule, N)] = d
            print("C", rule, N, json.dumps(d), flush=True)
    json.dump(res, open(out + "/metrics_cd4_donor_split.json", "w"), indent=1)
    print("saved", flush=True)


if __name__ == "__main__":
    main()
