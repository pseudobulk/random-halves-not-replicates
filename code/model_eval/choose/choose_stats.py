"""CHOOSE: formal tests for the clone-halves result and the formula (cell-type reference, as the main analysis of
choose_clone.py; same clone calls, genes and deviations).
1. Permutation test. Statistic T = mean over perturbations of (random-half signal - clone-half signal) / V. Null: clone
   labels shuffled among the cells of each perturbation (clone sizes kept), which destroys any link between clone
   and expression; the clone halves of shuffled labels are then random groupings. B permutations; one-sided p.
2. Slope = ICC. Bootstrap over perturbations (2,000 draws) of the slope of strength on sum p_k^2 (through the origin)
   and of the median ANOVA ICC, and of their ratio.
3. Out of sample. Perturbations split at random into two halves (500 draws): ICC (median) from one half predicts the
   strength of the other as ICC x sum p_k^2; reported: Spearman and Pearson between predicted and measured, and the
   ratio of summed measured to summed predicted.
Half signals are computed from per-clone sums: with S_k the summed deviation of clone k, a half made of clone set A
has mean sum_{k in A} S_k / n_A, so every split costs a K x K Gram matrix.
Per perturbation, the random-half and clone-half signals average NS splits (argument; 50 in the first run, 1,000 since
10-23, when 50 splits left the sign of one small perturbation, DEAF1 with 64 cells, to the draw); perm_mean is the mean
clone-half signal over the B permutations.
Usage: choose_stats.py <export dir> <cbc dir> <cbc_libs.tsv> <out.json> [B] [NS]"""
import sys, os, json
import numpy as np
from scipy.stats import spearmanr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from choose_clone import load_expr, call_clones, MINC

NS = 50


def cross_from_groups(Sg, ng, side):
    a, b = side, ~side
    na, nb = ng[a].sum(), ng[b].sum()
    if na < 10 or nb < 10:
        return np.nan
    return float((Sg[a].sum(0) / na) @ (Sg[b].sum(0) / nb))


def main():
    exp, cbc, libs, outp = sys.argv[1:5]; B = int(sys.argv[5]) if len(sys.argv) > 5 else 200
    NS = int(sys.argv[6]) if len(sys.argv) > 6 else 50
    rng = np.random.default_rng(1)
    cells, genes, X, meta = load_expr(exp)
    clone, _ = call_clones(cells, cbc, libs)
    C = np.array([clone.get(c) for c in cells], dtype=object); has = np.array([c is not None for c in C])
    lib = np.array([int(c.split("_")[-1]) for c in cells]); P = meta.gRNA.astype(str).to_numpy()
    grp = meta.celltype_cl_coarse2.astype(str).to_numpy()
    gm = np.asarray(X.mean(0)).ravel() >= 0.1; Xd = X[:, gm].toarray()
    key = np.array(["%d|%s" % (a, b) for a, b in zip(lib, grp)]); uk, ki = np.unique(key, return_inverse=True)
    S = np.zeros((len(uk), Xd.shape[1])); np.add.at(S, ki, Xd); N = np.bincount(ki, minlength=len(uk)).astype(float)
    V = float(((Xd - S[ki] / N[ki, None]) ** 2).sum() / len(Xd))
    per = {}
    for p in np.unique(P):
        idx = np.where((P == p) & has)[0]
        if len(idx) < MINC:
            continue
        allp = np.where(P == p)[0]
        Sp = np.zeros_like(S); np.add.at(Sp, ki[allp], Xd[allp]); Np = np.bincount(ki[allp], minlength=len(uk))
        E = Xd[idx] - ((S - Sp) / np.maximum(N - Np, 1)[:, None])[ki[idx]]
        uc, ci = np.unique(C[idx], return_inverse=True); K = len(uc); nk = np.bincount(ci).astype(float)
        Sk = np.zeros((K, E.shape[1])); np.add.at(Sk, ci, E)
        rnd = []
        for _ in range(NS):
            bit = rng.integers(0, 2, len(idx)).astype(bool)
            if bit.sum() >= 10 and (~bit).sum() >= 10:
                rnd.append(float(E[bit].mean(0) @ E[~bit].mean(0)))
        def clone_halves(Sg, ng, ns):
            v = [cross_from_groups(Sg, ng, rng.integers(0, 2, len(ng)).astype(bool)) for _ in range(ns)]
            v = [x for x in v if np.isfinite(x)]
            return float(np.mean(v)) if v else np.nan
        obs = clone_halves(Sk, nk, NS)
        perm = []
        for _ in range(B):
            cp = rng.permutation(ci); Sq = np.zeros_like(Sk); np.add.at(Sq, cp, E)
            perm.append(clone_halves(Sq, nk, 20))
        pk = nk / nk.sum(); m = E.mean(0)
        ssb = float((nk * ((Sk / nk[:, None] - m) ** 2).sum(1)).sum()); ssw = float(((E - (Sk / nk[:, None])[ci]) ** 2).sum())
        Nn = len(idx); n0 = (Nn - (nk ** 2).sum() / Nn) / (K - 1)
        vcl = (ssb / (K - 1) - ssw / (Nn - K)) / n0; icc = vcl / (vcl + ssw / (Nn - K))
        per[p] = dict(n=int(Nn), random=float(np.mean(rnd)), clone=obs, perm=perm, perm_mean=float(np.nanmean(perm)), sum_p2=float((pk ** 2).sum()), ICC=float(icc))
        print(p, Nn, round(per[p]["random"], 3), round(obs, 3), round(float(np.nanmean(perm)), 3), flush=True)
    ps = [p for p in per if np.isfinite(per[p]["clone"])]
    rnd = np.array([per[p]["random"] for p in ps]); cl = np.array([per[p]["clone"] for p in ps])
    PM = np.array([per[p]["perm"] for p in ps])                       # perturbations x B
    T = float((rnd - cl).mean() / V); Tb = (rnd[:, None] - PM).mean(0) / V
    out = {"rules": __doc__, "V": V, "n_perturbations": len(ps), "permutation": dict(
        T=T, T_perm_mean=float(np.nanmean(Tb)), T_perm_sd=float(np.nanstd(Tb)), B=B,
        p_one_sided=float((1 + np.sum(Tb >= T)) / (1 + np.isfinite(Tb).sum())),
        mean_clone_halves=float(cl.mean()), mean_shuffled_clone_halves=float(np.nanmean(PM)), mean_random=float(rnd.mean()))}
    s = (rnd - cl) / V; p2 = np.array([per[p]["sum_p2"] for p in ps]); icc = np.array([per[p]["ICC"] for p in ps])
    slope = lambda i: float((s[i] * p2[i]).sum() / (p2[i] ** 2).sum())
    bs, bi, br = [], [], []
    for _ in range(2000):
        i = rng.integers(0, len(ps), len(ps)); a, b = slope(i), float(np.median(icc[i])); bs.append(a); bi.append(b); br.append(a / b)
    q = lambda v: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
    out["slope_vs_icc"] = dict(slope=slope(np.arange(len(ps))), slope_ci=q(bs), median_icc=float(np.median(icc)), icc_ci=q(bi),
                               ratio_ci=q(br))
    sp, pe, ra = [], [], []
    for _ in range(500):
        perm = rng.permutation(len(ps)); a, b = perm[: len(ps) // 2], perm[len(ps) // 2:]
        pred = np.median(icc[a]) * p2[b]
        sp.append(spearmanr(pred, s[b])[0]); pe.append(np.corrcoef(pred, s[b])[0, 1]); ra.append(s[b].sum() / pred.sum())
    out["out_of_sample"] = dict(spearman_median=float(np.median(sp)), pearson_median=float(np.median(pe)),
                                ratio_measured_to_predicted_median=float(np.median(ra)), ratio_iqr=[float(np.percentile(ra, 25)), float(np.percentile(ra, 75))])
    out["perturbations"] = {p: {k: v for k, v in per[p].items() if k != "perm"} for p in ps}
    print(json.dumps({k: v for k, v in out.items() if k not in ("rules", "perturbations")}), flush=True)
    json.dump(out, open(outp, "w"), indent=1)


if __name__ == "__main__":
    main()
