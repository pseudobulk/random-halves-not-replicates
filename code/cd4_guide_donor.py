"""Same guide vs different guide across independent donors (A100; preprint review: can guide-specific effects explain
what random halves see and guide halves do not?). Released 4-donor pseudobulk of Zhu 2026 (GWCD4i.pseudobulk_merged).
For each target, four pseudobulk effects per state (resting, stimulated 48 h): guide rank g in {0, 1} x donor set {1, 2}
or {3, 4}; both sequencing runs pooled; log1p CP10K of the summed counts minus that of the summed non-targeting controls
of the same donors and state; at least MIN_CELLS cells each. Two designs pair the same four effects:
  same guide   (g0 D12, g0 D34) and (g1 D12, g1 D34): independent transductions and donors, same guide
  diff guide   (g0 D12, g1 D34) and (g1 D12, g0 D34): independent transductions and donors, different guides
Both designs exclude clones and donor-specific responses; they differ only in whether the two halves share the guide.
The difference of their cross-half signals estimates the energy of guide-specific effects (off-target activity, guide
efficiency, guide-level effects at the target locus), which the within-donor guide halves of the single-cell analysis
cannot separate from clones. Reported per design (mean over its two pairings): unexpressed targets' share of positive
cross-half signal and mean signal (null baseline), expressed targets' mean signal (null baseline), context-specific
share and the shared-to-specific predictor (k = 20, template removed; cd4_inert2.stats) with its shuffled control.
Unexpressed = pooled control mean < 0.01 in both states; expressed >= 0.5 in both. Genes: control >= 0.1 in every donor
and state. Usage: cd4_guide_donor.py <pseudobulk_merged.h5ad> <out_dir>. Exploratory."""
import os, sys, json
import numpy as np
import h5py
from scipy.sparse import csr_matrix
from multiprocessing import Pool
from cd4_m3 import cat, MIN_CELLS, MIN_EXPR
from cd4_inert2 import stats

DESIGNS = {"same_guide": ((0, 1), (2, 3)), "diff_guide": ((0, 3), (2, 1))}
NB = 1000   # bootstrap over targets (10-18 review: intervals for the guide-specific ratio and the null excess)
G_ = {}


def pair_stats(args):
    a, b = args
    E, inert, expd = G_["E"], G_["inert"], G_["expd"]
    XR = np.stack([E[:, a, 0], E[:, b, 0]]).astype(np.float64); XS = np.stack([E[:, a, 1], E[:, b, 1]]).astype(np.float64)
    out = {}
    XRc, XSc = XR - XR[:, inert].mean(1, keepdims=True), XS - XS[:, inert].mean(1, keepdims=True)
    for st, X in (("rest", XRc), ("stim", XSc)):
        cr = (X[0] * X[1]).sum(1)
        out["inert_frac_pos_" + st] = float((cr[inert] > 0).mean()); out["inert_mean_" + st] = float(cr[inert].mean())
        out["expressed_mean_" + st] = float(cr[expd].mean()); out["expressed_frac_pos_" + st] = float((cr[expd] > 0).mean())
    XRb, XSb = XR - XR.mean(1, keepdims=True), XS - XS.mean(1, keepdims=True)
    s = stats(XRb, XSb, inert, expd)
    out.update(gamma_share=s["gamma_share_of_specific"], readout_k20=s["readout_k20"], readout_k20_shuffled=s["readout_k20_shuffled"])
    bA, bB = (XRb[0] + XSb[0]) / 2, (XRb[1] + XSb[1]) / 2; gA, gB = XRb[0] - bA, XRb[1] - bB
    out["beta_cross_expressed"] = float((bA[expd] * bB[expd]).sum(1).mean()); out["gamma_cross_expressed"] = float((gA[expd] * gB[expd]).sum(1).mean())
    print("pair", a, b, json.dumps(out), flush=True)
    arrays = {st: (X[0] * X[1]).sum(1) for st, X in (("rest", XRc), ("stim", XSc))}
    return out, arrays


def main():
    path, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    with h5py.File(path) as f:
        o = f["obs"]
        donors, dc = cat(o, "donor_id"); conds, cc = cat(o, "culture_condition")
        guides, gc = cat(o, "guide_id")
        genes_t, tc = cat(o, "perturbed_gene_name"); gtypes, gtc = cat(o, "guide_type")
        ncell = o["n_cells"][:]
        gname = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var"]["gene_name"][:]])
        X = csr_matrix((f["X"]["data"][:], f["X"]["indices"][:], f["X"]["indptr"][:]), shape=tuple(f["X"].attrs["shape"]))
    ctrl_type = [i for i, g in enumerate(gtypes) if "non" in g.lower()][0]
    isc = gtc == ctrl_type
    D, S, R = len(donors), len(conds), 2
    T = len(genes_t)
    tgt = np.where(isc, T, tc)
    up, pinv = np.unique(tgt * (gc.max() + 1) + gc, return_inverse=True)       # guide rank within target, as cd4_donor_split
    ut = up // (gc.max() + 1); first = np.searchsorted(ut, ut, side="left")
    rc = ((np.arange(len(up)) - first) % 2)[pinv]
    key = ((tgt * D + dc) * S + cc) * R + rc
    uk, inv = np.unique(key, return_inverse=True)
    agg = (csr_matrix((np.ones(len(key)), (inv, np.arange(len(key)))), shape=(len(uk), len(key))) @ X).tocsr()
    cells = np.bincount(inv, weights=ncell, minlength=len(uk))
    pos = {k: i for i, k in enumerate(uk)}
    row = lambda t, d, s, r: pos.get(((t * D + d) * S + s) * R + r, -1)
    ctrl_mean = []
    for d in range(D):
        for s in range(S):
            rr = [row(T, d, s, r) for r in range(R) if row(T, d, s, r) >= 0]
            v = np.asarray(agg[rr].sum(0)).ravel(); ctrl_mean.append(np.log1p(v / v.sum() * 1e4))
    gk = np.where(np.min(ctrl_mean, 0) >= MIN_EXPR)[0]; gpos = {g: k for k, g in enumerate(gname[gk])}
    S_rest, S_stim = list(conds).index("Rest"), list(conds).index("Stim48hr")
    cache = {}
    def ctrl(ds, s):
        k = (tuple(ds), s)
        if k not in cache:
            rr = [row(T, d, s, r) for d in ds for r in range(R) if row(T, d, s, r) >= 0]
            v = np.asarray(agg[rr].sum(0)).ravel(); cache[k] = np.log1p(v / v.sum() * 1e4)[gk]
        return cache[k]
    def eff(t, ds, s, r):
        rt = [row(t, d, s, r) for d in ds if row(t, d, s, r) >= 0]
        if not rt or cells[rt].sum() < MIN_CELLS:
            return None
        a = np.asarray(agg[rt].sum(0)).ravel()
        e = np.log1p(a / a.sum() * 1e4)[gk] - ctrl(ds, s)
        if genes_t[t] in gpos:
            e[gpos[genes_t[t]]] = 0
        return e
    full = {}
    for s in (S_rest, S_stim):
        rr = [row(T, d, s, r) for d in range(D) for r in range(R) if row(T, d, s, r) >= 0]
        v = np.asarray(agg[rr].sum(0)).ravel(); full[s] = np.log1p(v / v.sum() * 1e4)
    gidx = {g: j for j, g in enumerate(gname)}
    cellsets = [(0, [0, 1]), (0, [2, 3]), (1, [0, 1]), (1, [2, 3])]
    keepT, E = [], []
    for t in range(T):
        if genes_t[t] not in gidx:
            continue
        es = [[eff(t, ds, s, r) for s in (S_rest, S_stim)] for r, ds in cellsets]
        if any(x is None for rw in es for x in rw):
            continue
        keepT.append(t); E.append(np.array(es, dtype=np.float32))
    E = np.array(E)
    tx = np.array([[full[S_rest][gidx[genes_t[t]]], full[S_stim][gidx[genes_t[t]]]] for t in keepT])
    inert, expd = tx.max(1) < 0.01, tx.min(1) >= 0.5
    print("targets", len(keepT), "inert", int(inert.sum()), "expressed", int(expd.sum()), "genes", len(gk), flush=True)
    G_.update(E=E, inert=inert, expd=expd)
    jobs = [p for ps in DESIGNS.values() for p in ps]
    with Pool(4) as pool:
        both = pool.map(pair_stats, jobs)
    per = dict(zip(jobs, [x[0] for x in both])); arr = dict(zip(jobs, [x[1] for x in both]))
    # per-target cross-half signals (mean over the two pairings of each design) and a bootstrap over targets
    cr = {dn: {st: np.mean([arr[p][st] for p in ps], 0) for st in ("rest", "stim")} for dn, ps in DESIGNS.items()}
    np.savez_compressed(out + "/per_target_cross.npz", targets=np.array([genes_t[t] for t in keepT]), inert=inert, expressed=expd,
                        **{"%s_%s" % (dn, st): cr[dn][st] for dn in cr for st in cr[dn]})
    rng = np.random.default_rng(7); ie, ii = np.where(expd)[0], np.where(inert)[0]; boot = {}
    for st in ("rest", "stim"):
        sg, dg = cr["same_guide"][st], cr["diff_guide"][st]; rat, dpos = [], []
        for _ in range(NB):
            b = rng.choice(ie, len(ie)); rat.append((sg[b].mean() - dg[b].mean()) / dg[b].mean())
            c = rng.choice(ii, len(ii)); dpos.append((sg[c] > 0).mean() - (dg[c] > 0).mean())
        boot["ratio_" + st] = [float(np.percentile(rat, 2.5)), float(np.percentile(rat, 97.5))]
        boot["inert_frac_pos_diff_" + st] = [float(np.mean(dpos)), float(np.percentile(dpos, 2.5)), float(np.percentile(dpos, 97.5))]
    print("bootstrap", json.dumps(boot), flush=True)
    res = {"rules": __doc__, "n_targets": len(keepT), "n_inert": int(inert.sum()), "n_expressed": int(expd.sum()), "n_genes": int(len(gk))}
    for dn, ps in DESIGNS.items():
        res[dn] = {k: float(np.mean([per[p][k] for p in ps])) for k in per[ps[0]]}
        res[dn + "_pairings"] = [per[p] for p in ps]
    g = res["same_guide"]["expressed_mean_rest"] - res["diff_guide"]["expressed_mean_rest"]
    res["guide_specific_over_diff_guide_rest"] = float(g / res["diff_guide"]["expressed_mean_rest"])
    g2 = res["same_guide"]["expressed_mean_stim"] - res["diff_guide"]["expressed_mean_stim"]
    res["guide_specific_over_diff_guide_stim"] = float(g2 / res["diff_guide"]["expressed_mean_stim"])
    res["bootstrap_over_targets"] = boot
    json.dump(res, open(out + "/metrics_cd4_guide_donor.json", "w"), indent=1)
    print("done", json.dumps({k: v for k, v in res.items() if not k.endswith("pairings") and k != "rules"}), flush=True)


if __name__ == "__main__":
    main()
