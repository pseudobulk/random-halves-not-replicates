"""Model evaluation step 1: per-donor scores with standard errors, from saved predictions, by the pilot's scoring rule
(pilot_common.score: Pearson over genes of predicted vs observed stimulated effect, target gene excluded; truth A =
the guide's own stimulated cells, truth B = the other guide of the same target).
Standard errors are clustered by target: a target's two guides share the target and swap truths (A of one guide is
B of the other), so the guide-level values are first averaged within target, then the SE is taken over targets.
Models: the four baselines (baseline_preds.npz of pilot_baselines.py) and STATE Transition (adata_pred.h5ad, predicted
effect relative to the predicted non-targeting stimulated cells, as pilot_state_score.py). Also the reference
agreement of the two truths. Groups: inert, strong, and all test targets together.
Usage (env PILOT_DONOR or PILOT_STATES): me_score_detail.py <data_dir> <baseline_preds.npz> <adata_pred.h5ad | -> <out_json>"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "4")
import sys, json
import numpy as np
from pilot_common import STATES, DONOR, effects, score, test_pairs, kinds

BASELINES = ("mean", "copy", "additive", "linear")


def per_guide(pred, des, res, genes):
    """guide-level Pearson against A and B, same rule as pilot_common.score"""
    S = res[STATES[1]]; si = {l: i for i, l in enumerate(S["labels"])}; gi = {g: j for j, g in enumerate(genes)}
    rows = []
    for t, kind, g, o in test_pairs(des):
        if g not in pred or g not in si or o not in si:
            continue
        keep = np.ones(len(genes), bool)
        if t in gi:
            keep[gi[t]] = False
        p = pred[g][keep]
        r = [np.corrcoef(p, S["eff"][si[lab]][keep])[0, 1] if p.std() > 0 else np.nan for lab in (g, o)]
        rows.append((t, kind, r[0], r[1]))
    return rows


def summarise(rows, groups):
    out = {}
    for grp in tuple(groups) + ("all",):
        rr = [r for r in rows if grp == "all" or r[1] == grp]
        if not rr:
            continue
        tg = sorted({r[0] for r in rr})
        M = {k: np.array([np.nanmean([r[i] for r in rr if r[0] == t]) for t in tg]) for k, i in (("A", 2), ("B", 3))}
        M["AB"] = M["A"] - M["B"]
        a, b = np.array([r[2] for r in rr]), np.array([r[3] for r in rr])
        o = dict(n_guides=len(rr), n_targets=len(tg), A=float(np.nanmean(a)), B=float(np.nanmean(b)), AB=float(np.nanmean(a - b)))
        for k in ("A", "B", "AB"):   # means over guides (as the pilot); SEs over target means
            o[k + "_se"] = float(np.nanstd(M[k], ddof=1) / np.sqrt(np.isfinite(M[k]).sum()))
        out[grp] = o
    return out


def main():
    d, bp, sp, outp = sys.argv[1:5]
    des, res = effects(d); genes = res[STATES[1]]["genes"]
    preds = {m: {} for m in BASELINES}
    z = np.load(bp)
    for k in z.files:
        m, g = k.split("__", 1); preds[m][g] = z[k]
    if sp != "-":   # '-' = baselines only
        import anndata as ad
        p = ad.read_h5ad(sp)
        X = p.X.toarray() if hasattr(p.X, "toarray") else np.asarray(p.X)
        assert X.shape[1] == len(genes)
        lab = p.obs["gene"].astype(str).to_numpy()
        stim = (p.obs["cell_type"].astype(str).to_numpy() == "stim") if "cell_type" in p.obs else np.ones(len(lab), bool)
        ctrl = stim & (lab == "non-targeting"); assert ctrl.sum() > 0
        base = X[ctrl].mean(0)
        preds["state_transition"] = {g: X[stim & (lab == g)].mean(0) - base for g in np.unique(lab[stim]) if g != "non-targeting"}
    S = res[STATES[1]]; si = {l: i for i, l in enumerate(S["labels"])}
    preds["truth_A_vs_B"] = {g: S["eff"][si[g]] for _, _, g, _ in test_pairs(des) if g in si}
    out = dict(donor=DONOR, n_test_targets={k: len(des[k]) for k in kinds(des)},
               n_dropped_test=len(des["dropped_test"]), n_genes=len(genes), models={}, rules=__doc__)
    for m, pr in preds.items():
        rows = per_guide(pr, des, res, genes)
        out["models"][m] = summarise(rows, kinds(des))
        out["models"][m]["per_guide"] = [[t, k, float(a), float(b)] for t, k, a, b in rows]   # target, group, A, B
        ref = score(pr, des, res, genes)   # identical means by construction; checked
        for kind in kinds(des):
            if ref[kind]["n"]:
                assert abs(ref[kind]["leak"] - out["models"][m][kind]["AB"]) < 1e-9
        v = out["models"][m]["all"]
        print("%s %-17s A %.3f (%.3f)  B %.3f (%.3f)  A-B %.3f (%.3f)" % (DONOR, m, v["A"], v["A_se"], v["B"], v["B_se"], v["AB"], v["AB_se"]), flush=True)
    json.dump(out, open(outp, "w"), indent=1)


if __name__ == "__main__":
    main()
