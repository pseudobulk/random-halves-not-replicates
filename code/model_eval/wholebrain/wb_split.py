"""Whole-brain: random halves against mouse halves (wb_mouse_halves.py), the independent replicate of this design.
Null test and targeting signal: wb_analyze.null_test (within-mouse reference) on the random-half and the mouse-half sums.
Cross-context Glut -> GABA with mouse halves: models trained on group-0 mice (Glut effect -> GABA effect of the training
guides); for a held-out guide g0 (other guide g1 of the same target) the prediction from its group-0 Glut effect is
scored against
  A    GABA effect of g0 in group-0 mice (same guide, same mice: what a single-mouse-set evaluation uses)
  A'   GABA effect of g0 in group-1 mice (same guide, other mice)
  B'   GABA effect of g1 in group-1 mice (other guide, other mice)
A - A' is the part shared only through the mice; A' - B' the guide-specific part. Units need >= MINH nuclei in each
class and half; otherwise as wb_analyze.ab_test.
Usage: wb_split.py <agg_random.npz> <agg_mouse.npz> <out.json> [seed]"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wb_analyze import null_test, effects

MINH = 15


def gmasks(z):
    nmouse = len(z["mice"]); out = {}
    for c in (0, 1):
        Mc = z["M"].reshape(nmouse, 2, 2, -1)[:, c].sum((0, 1)); mc = z["m"].reshape(nmouse, 2, 2)[:, c].sum()
        out[c] = (Mc / mc >= 0.1, z["Q"][c] / mc - (Mc / mc) ** 2)
    return out


def ab_split(z, gmask2, rng):
    n, targets = z["n"], z["targets"]; ng = len(z["guides"]); G = np.arange(ng)
    E = {c: effects(z, c, list(G * 2 + c), "within")[..., gmask2] for c in (0, 1)}     # (guide, half, gene)
    ok = (n[G * 2] >= MINH).all(1) & (n[G * 2 + 1] >= MINH).all(1)
    genes = z["genes"][gmask2]; gi = {g: j for j, g in enumerate(genes)}
    by_t = {}
    for g in np.where(ok)[0]:
        if targets[g] != "Non_target":
            by_t.setdefault(targets[g], []).append(g)
    test = {t: sorted(gs, key=lambda g: -n[g * 2 + 1].sum())[:2] for t, gs in by_t.items() if len(gs) >= 2}
    tnames = sorted(test); fold = dict(zip(tnames, rng.integers(0, 2, len(tnames))))
    usable = np.where(ok & (targets != "Non_target"))[0]
    res = {m: {"A": [], "A1": [], "B1": []} for m in ("mean", "copy", "linear")}; ref = {"A_A1": [], "A1_B1": []}
    for f in (0, 1):
        held = {t for t in tnames if fold[t] == f}
        tr = np.array([g for g in usable if targets[g] not in held])
        X, Y = E[0][tr, 0], E[1][tr, 0]; mx, my = X.mean(0), Y.mean(0)
        _, _, Vt = np.linalg.svd(X - mx, full_matrices=False); V = Vt[:50].T; Z = (X - mx) @ V
        folds = rng.integers(0, 5, len(tr)); best = None
        for lam in (1, 10, 100, 1000, 10000):
            err = 0.0
            for k in range(5):
                a, b = folds != k, folds == k
                Wt = np.linalg.solve(Z[a].T @ Z[a] + lam * np.eye(50), Z[a].T @ (Y[a] - my))
                err += ((Z[b] @ Wt + my - Y[b]) ** 2).sum()
            best = (err, lam) if best is None or err < best[0] else best
        Wt = np.linalg.solve(Z.T @ Z + best[1] * np.eye(50), Z.T @ (Y - my))
        for t in held:
            for gA, gB in (test[t], test[t][::-1]):
                keep = np.ones(len(genes), bool)
                if t in gi:
                    keep[gi[t]] = False
                pred = {"mean": my, "copy": E[0][gA, 0], "linear": ((E[0][gA, 0] - mx) @ V) @ Wt + my}
                tA, tA1, tB1 = E[1][gA, 0][keep], E[1][gA, 1][keep], E[1][gB, 1][keep]
                for m, p in pred.items():
                    for k, y in (("A", tA), ("A1", tA1), ("B1", tB1)):
                        res[m][k].append(np.corrcoef(p[keep], y)[0, 1])
                ref["A_A1"].append(np.corrcoef(tA, tA1)[0, 1]); ref["A1_B1"].append(np.corrcoef(tA1, tB1)[0, 1])
    out = {"n_test_targets": len(tnames), "truth_A_vs_A1": float(np.mean(ref["A_A1"])), "truth_A1_vs_B1": float(np.mean(ref["A1_B1"]))}
    for m, v in res.items():
        a, a1, b1 = (np.array(v[k]).reshape(-1, 2).mean(1) for k in ("A", "A1", "B1"))   # per target, two directions
        se = lambda d: float(d.std() / np.sqrt(len(d)))
        out[m] = dict(A=float(a.mean()), A1=float(a1.mean()), B1=float(b1.mean()), A_minus_A1=float((a - a1).mean()),
                      A_minus_A1_se=se(a - a1), A1_minus_B1=float((a1 - b1).mean()), A1_minus_B1_se=se(a1 - b1))
    return out


def main():
    pr, pm, outp = sys.argv[1:4]; rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 0)
    zr = dict(np.load(pr, allow_pickle=False)); zm = dict(np.load(pm, allow_pickle=False))
    out = {"rules": __doc__}; gm = gmasks(zr)
    for c, name in ((0, "Glut"), (1, "GABA")):
        mask, var = gm[c]; out[name] = {}
        for lab, z in (("random", zr), ("mouse", zm)):
            t = null_test(z, c, mask, "within", rng); t["mean_random_scaled"] = t["mean_random"] / float(var[mask].sum())
            t["mean_target_scaled"] = t["mean_target"] / float(var[mask].sum()); out[name][lab] = t
        r, m = out[name]["random"], out[name]["mouse"]
        out[name]["nt_mouse_over_random"] = m["mean_random"] / r["mean_random"]
        out[name]["target_mouse_over_random"] = m["mean_target"] / r["mean_target"]
        print(name, json.dumps(out[name]), flush=True)
    out["ab_mouse_halves"] = ab_split(zm, gm[0][0] & gm[1][0], rng)
    print("AB", json.dumps(out["ab_mouse_halves"]), flush=True)
    json.dump(out, open(outp, "w"), indent=1)


if __name__ == "__main__":
    main()
