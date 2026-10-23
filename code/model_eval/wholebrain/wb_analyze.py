"""Whole-brain in vivo Perturb-seq (AAV at P16 into post-mitotic neurons: no clonal expansion): the clone null test of
the preprint (Fig. 2c) and the cross-context A / B comparison of the model evaluation, from wb_aggregate.py sums.
Effects (per unit = guide x class, per random half): mean log1p CP10K minus a non-targeting reference, either
  within  the non-targeting mean of the same mouse, class and half, weighted by the unit's nuclei per mouse (removes
          differences between mice, which random halves of one guide share and two different guides do not), or
  global  the non-targeting mean of the class and half over all mice (as in a screen without a mouse column).
Genes: mean >= 0.1 over all non-targeting nuclei of the class. Units need >= MINC nuclei per half.
Null test (non-targeting guides, per class): cross-half signal sum_g x_A,g x_B,g under random halves; between two
different non-targeting guides (first half of one, second half of the other; 50 random pairings); share positive,
mean, and the mean scaled by the summed per-nucleus variance (comparable to vc_expr of the model evaluation); ratio of
the non-targeting random-half mean to that of targeting perturbations (grouped by target, centred over targets).
Cross-context A / B (Glut -> GABA): per guide, effects with both halves pooled (within-mouse reference), on genes
passing in both classes. Test targets: >= 2 guides with >= MINAB nuclei in both classes (g0, g1 = the two with most
GABA nuclei); split into 2 folds, training on all guides of the other fold's targets plus non-test targets. Models:
context mean, copy, linear (ridge on 50 PCs of the source effects, lambda by 5-fold CV). Truth A = the guide's own GABA
effect, B = the other guide's GABA effect; Pearson over genes without the target gene.
Usage: wb_analyze.py <agg.npz> <out.json> [seed]"""
import sys, json
import numpy as np

import os
MINC, MINAB, NPAIR = int(os.environ.get("WB_MINC", 10)), int(os.environ.get("WB_MINAB", 30)), 50   # env only for smoke tests


def refs(z, c):
    """per unit-half non-targeting reference, within mouse and global; arrays (U, 2, G)"""
    M, m, nm = z["M"], z["m"], z["nm"]; nmouse = len(z["mice"])
    Mc = M.reshape(nmouse, 2, 2, -1)[:, c]; mc = m.reshape(nmouse, 2, 2)[:, c]          # (mouse, half, G), (mouse, half)
    glob = Mc.sum(0) / mc.sum(0)[:, None]                                               # (half, G)
    within_m = np.where(mc[..., None] > 0, Mc / np.maximum(mc, 1)[..., None], glob[None])  # (mouse, half, G)
    return glob, within_m


def effects(z, c, units, mode):
    S, n, nm = z["S"], z["n"], z["nm"]
    glob, wm = refs(z, c)
    E = np.zeros((len(units), 2, S.shape[2]), np.float64)
    for j, u in enumerate(units):
        for h in (0, 1):
            mean = S[u, h] / max(n[u, h], 1)
            if mode == "global":
                ref = glob[h]
            else:
                w = nm[u, h] / max(nm[u, h].sum(), 1); ref = w @ wm[:, h]
            E[j, h] = mean - ref
    return E


def cross(a, b):
    return (a * b).sum(-1)


def null_test(z, c, gmask, mode, rng):
    targets = z["targets"]; ng = len(z["guides"]); n = z["n"]
    U = np.arange(ng) * 2 + c
    nt = [u for u in U if targets[u // 2] == "Non_target" and (n[u] >= MINC).all()]
    tg = [u for u in U if targets[u // 2] != "Non_target" and (n[u] >= MINC).all()]
    En = effects(z, c, nt, mode)[..., gmask]
    rnd = cross(En[:, 0], En[:, 1])
    pairs = []
    for _ in range(NPAIR):
        p = rng.permutation(len(nt)); a, b = p[: len(p) // 2], p[len(p) // 2: 2 * (len(p) // 2)]
        pairs.append(cross(En[a, 0], En[b, 1]))
    diff = np.concatenate(pairs)
    # targeting perturbations grouped by target gene
    S, nm = z["S"], z["nm"]
    tnames = sorted({targets[u // 2] for u in tg})
    Et = []
    glob, wm = refs(z, c)
    for t in tnames:
        us = [u for u in tg if targets[u // 2] == t]
        e = []
        for h in (0, 1):
            ssum = S[us, h].sum(0); nn = n[us, h].sum()
            if mode == "global":
                ref = glob[h]
            else:
                w = nm[us, h].sum(0) / max(nm[us, h].sum(), 1); ref = w @ wm[:, h]
            e.append(ssum / nn - ref)
        Et.append(e)
    if not Et:
        tc = np.array([np.nan])
    else:
        Et = np.array(Et)[..., gmask]; Et = Et - Et.mean(0, keepdims=True)
        tc = cross(Et[:, 0], Et[:, 1])
    return dict(n_nt=len(nt), n_targets=len(tnames), share_pos_random=float((rnd > 0).mean()),
                share_pos_diff_guides=float((diff > 0).mean()), mean_random=float(rnd.mean()), mean_diff_guides=float(diff.mean()),
                se_random=float(rnd.std() / np.sqrt(len(rnd))), mean_target=float(tc.mean()),
                ratio_nt_to_target=float(rnd.mean() / tc.mean()))


def ab_test(z, gmask2, rng):
    n, S, targets = z["n"], z["S"], z["targets"]; ng = len(z["guides"])
    cnt = n.sum(1)                                   # nuclei per unit (both halves)
    E = {c: effects(z, c, list(np.arange(ng) * 2 + c), "within") for c in (0, 1)}
    W = {c: ((E[c][:, 0] * n[np.arange(ng) * 2 + c, 0, None] + E[c][:, 1] * n[np.arange(ng) * 2 + c, 1, None])
             / np.maximum(cnt[np.arange(ng) * 2 + c], 1)[:, None])[:, gmask2] for c in (0, 1)}
    ok = (cnt[np.arange(ng) * 2] >= MINAB) & (cnt[np.arange(ng) * 2 + 1] >= MINAB)
    genes = z["genes"][gmask2]; gi = {g: j for j, g in enumerate(genes)}
    by_t = {}
    for g in np.where(ok)[0]:
        t = targets[g]
        if t != "Non_target":
            by_t.setdefault(t, []).append(g)
    test = {t: sorted(gs, key=lambda g: -cnt[g * 2 + 1])[:2] for t, gs in by_t.items() if len(gs) >= 2}
    tnames = sorted(test); fold = dict(zip(tnames, rng.integers(0, 2, len(tnames))))
    usable = np.where(ok & (targets != "Non_target"))[0]
    res = {m: {"A": [], "B": []} for m in ("mean", "copy", "linear")}; ref = []
    for f in (0, 1):
        held = {t for t in tnames if fold[t] == f}
        tr = np.array([g for g in usable if targets[g] not in held])
        X, Y = W[0][tr], W[1][tr]; mx, my = X.mean(0), Y.mean(0)
        U_, s_, Vt = np.linalg.svd(X - mx, full_matrices=False); V = Vt[:50].T; Z = (X - mx) @ V
        folds = rng.integers(0, 5, len(tr)); best = None
        for lam in (1, 10, 100, 1000, 10000):
            err = 0.0
            for k in range(5):
                a, b = folds != k, folds == k
                Wt = np.linalg.solve(Z[a].T @ Z[a] + lam * np.eye(50), Z[a].T @ (Y[a] - my))
                err += ((Z[b] @ Wt + my - Y[b]) ** 2).sum()
            best = (err, lam) if best is None or err < best[0] else best
        lam = best[1]; Wt = np.linalg.solve(Z.T @ Z + lam * np.eye(50), Z.T @ (Y - my))
        for t in held:
            for gA, gB in (test[t], test[t][::-1]):
                keep = np.ones(len(genes), bool)
                if t in gi:
                    keep[gi[t]] = False
                pred = {"mean": my, "copy": W[0][gA], "linear": ((W[0][gA] - mx) @ V) @ Wt + my}
                yA, yB = W[1][gA][keep], W[1][gB][keep]
                for m, p in pred.items():
                    res[m]["A"].append(np.corrcoef(p[keep], yA)[0, 1]); res[m]["B"].append(np.corrcoef(p[keep], yB)[0, 1])
                ref.append(np.corrcoef(yA, yB)[0, 1])
    out = {"n_test_targets": len(tnames), "truth_A_vs_B": float(np.mean(ref))}
    for m, v in res.items():
        a, b = np.array(v["A"]), np.array(v["B"]); d = (a - b).reshape(-1, 2).mean(1)   # per target (two directions)
        out[m] = dict(A=float(a.mean()), B=float(b.mean()), A_minus_B=float(d.mean()), se=float(d.std() / np.sqrt(len(d))))
    return out


def main():
    aggp, outp = sys.argv[1:3]; seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    z = dict(np.load(aggp, allow_pickle=False)); rng = np.random.default_rng(seed)
    out = {"rules": __doc__, "n_nuclei": int(z["n"].sum()), "n_mice": int(len(z["mice"]))}
    nmouse = len(z["mice"]); gm = {}
    for c, name in ((0, "Glut"), (1, "GABA")):
        Mc = z["M"].reshape(nmouse, 2, 2, -1)[:, c].sum((0, 1)); mc = z["m"].reshape(nmouse, 2, 2)[:, c].sum()
        mean = Mc / mc; gmask = mean >= 0.1; gm[c] = gmask
        var = z["Q"][c] / mc - mean ** 2
        r = {"n_genes": int(gmask.sum()), "n_nt_nuclei": int(mc)}
        for mode in ("within", "global"):
            t = null_test(z, c, gmask, mode, rng)
            t["mean_random_scaled"] = t["mean_random"] / float(var[gmask].sum())
            r[mode] = t
        out[name] = r
        print(name, json.dumps(r), flush=True)
    out["ab_Glut_to_GABA"] = ab_test(z, gm[0] & gm[1], rng)
    print("AB", json.dumps(out["ab_Glut_to_GABA"]), flush=True)
    json.dump(out, open(outp, "w"), indent=1)


if __name__ == "__main__":
    main()
