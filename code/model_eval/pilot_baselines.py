"""Model pilot: baselines that predict a guide's stimulated effect from its resting effect (pilot_common.py).
  mean      mean stimulated effect over training guides (no perturbation information; same for both guides)
  copy      the guide's own resting effect
  additive  resting effect plus the mean stimulated-minus-resting shift over training guides
  linear    ridge map from resting to stimulated effects, fit on training guides (resting effects reduced to 50
            principal components; ridge strength by 5-fold cross-validation on the training guides)
Training pairs use the same guide's resting and stimulated cells, as any model trained on these data does.
Also reports the reference agreement between the two truths (A vs B). Usage: pilot_baselines.py <data_dir> <out_json>"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "8")
import sys, json
import numpy as np
from pilot_common import STATES, effects, score, test_pairs


def ridge_fit(X, Y, lam):
    p = X.shape[1]
    return np.linalg.solve(X.T @ X + lam * np.eye(p), X.T @ Y)


def main():
    d, outp = sys.argv[1:3]
    des, res = effects(d)
    R, S = res[STATES[0]], res[STATES[1]]; genes = R["genes"]
    ri, si = {l: i for i, l in enumerate(R["labels"])}, {l: i for i, l in enumerate(S["labels"])}
    train = [l for l in S["labels"] if l.split("|")[0] in set(des["train"]) and l in ri]
    Xtr = np.stack([R["eff"][ri[l]] for l in train]); Ytr = np.stack([S["eff"][si[l]] for l in train])
    tests = [g for _, _, g, _ in test_pairs(des) if g in ri]
    preds = {}
    preds["mean"] = {g: Ytr.mean(0) for g in tests}
    preds["copy"] = {g: R["eff"][ri[g]] for g in tests}
    shift = (Ytr - Xtr).mean(0)
    preds["additive"] = {g: R["eff"][ri[g]] + shift for g in tests}
    # linear: PCA of resting effects, ridge to stimulated effects
    mu_x, mu_y = Xtr.mean(0), Ytr.mean(0)
    U, s, Vt = np.linalg.svd(Xtr - mu_x, full_matrices=False); V = Vt[:50].T
    Z = (Xtr - mu_x) @ V
    rng = np.random.default_rng(0); fold = rng.integers(0, 5, len(train)); cv = {}
    for lam in (1e-2, 1e-1, 1, 10, 100, 1000):
        err = 0.0
        for k in range(5):
            tr, te = fold != k, fold == k
            W = ridge_fit(Z[tr], Ytr[tr] - mu_y, lam)
            err += ((Z[te] @ W + mu_y - Ytr[te]) ** 2).sum()
        cv[lam] = err
    lam = min(cv, key=cv.get); W = ridge_fit(Z, Ytr - mu_y, lam)
    preds["linear"] = {g: ((R["eff"][ri[g]] - mu_x) @ V) @ W + mu_y for g in tests}
    out = {"n_train_guides": len(train), "n_test_guides": len(tests), "ridge_lambda": lam, "models": {}}
    for m, p in preds.items():
        out["models"][m] = score(p, des, res, genes)
    # reference: agreement of the two truths themselves (stimulated guide g vs the other guide)
    out["truth_A_vs_B"] = score({g: S["eff"][si[g]] for g in tests if g in si}, des, res, genes)
    json.dump(out, open(outp, "w"), indent=1)
    for m, v in out["models"].items():
        print("%-9s" % m, "  ".join("%s A %.3f B %.3f leak %.3f" % (k, v[k]["pearson_A"], v[k]["pearson_B"], v[k]["leak"]) for k in ("inert", "strong")))
    np.savez_compressed(os.path.join(os.path.dirname(outp) or ".", "baseline_preds.npz"),
                        **{"%s__%s" % (m, g): v for m, p in preds.items() for g, v in p.items()})


if __name__ == "__main__":
    main()
