"""Model evaluation: fit check for STATE Transition. Does the trained model reproduce the per-guide stimulated effects
of TRAINING guides it was fit on? A model that cannot (as CPA and biolord in the pilot) can show neither skill nor
clone gain, so its A - B says nothing about leakage.
How: `state tx predict --toml <fit.toml>` with the trained checkpoint, where fit.toml lists NFIT training guides (seed 2,
not the validation guides) as the "test" set, plus NCHK of the real test guides as a check that the override reproduces
the saved test predictions (Pearson between the two predicted effects should be ~1). The training data are rebuilt
with pilot_state_prep.py (same rules and seeds as the run).
Predicted effect = mean over a guide's predicted stimulated cells minus mean over predicted non-targeting stimulated
cells; observed = pilot_common.effects; Pearson over genes (target gene excluded), mean and SE over guides. Off-diagonal
baselines: the same correlations with OTHER guides' observed effects / saved predictions (what a shared response alone,
or a label mix-up, would give). Re-prediction is not exact: STATE conditions on sampled control cells.
Usage (env PILOT_DONOR or PILOT_STATES):
  me_fitcheck.py toml <state_work_dir> <data_dir>           -> <state_work_dir>/fit.toml
  me_fitcheck.py score <data_dir> <adata_pred.h5ad> <preds.npz> <out_json>"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "4")
import sys, json
import numpy as np

NFIT, NCHK = 200, 20


def main():
    mode = sys.argv[1]
    if mode == "toml":
        wd, dd = sys.argv[2], sys.argv[3]
        from pilot_common import design, load, test_targets, STATES
        des = design(dd); prep = json.load(open(wd + "/prep.json")); val = sorted(prep["val"])
        present = set(load(dd, STATES[1])["lab"])   # training guide labels with stimulated cells (Replogle: g0 only)
        trains = sorted({"%s|g%d" % (t, p) for t in des["train"] for p in (0, 1)} & present - set(val))
        rng = np.random.default_rng(2)
        fit = sorted(rng.choice(trains, NFIT, replace=False).tolist())
        held = sorted({"%s|g%d" % (t, p) for t in test_targets(des) for p in (0, 1)})
        chk = held[:NCHK]
        q = lambda xs: "[" + ", ".join('"%s"' % x for x in xs) + "]"
        open(wd + "/fit.toml", "w").write('[datasets]\ncd4 = "%s/data/"\n\n[training]\ncd4 = "train"\n\n[fewshot]\n[fewshot."cd4.stim"]\nval = %s\ntest = %s\n' % (
            wd, q(val), q(fit + chk)))
        json.dump(dict(fit=fit, check=chk), open(wd + "/fit_guides.json", "w"), indent=1)
        print("fit guides", len(fit), "check guides", len(chk))
        return
    from pilot_common import STATES, DONOR, effects
    d, pf, pp, outp = sys.argv[2:6]
    des, res = effects(d); S = res[STATES[1]]; genes = S["genes"]; si = {l: i for i, l in enumerate(S["labels"])}
    gpos = {g: j for j, g in enumerate(genes)}; train = set(des["train"])
    import anndata as ad
    p = ad.read_h5ad(pf)
    X = p.X.toarray() if hasattr(p.X, "toarray") else np.asarray(p.X)
    assert X.shape[1] == len(genes)
    lab = p.obs["gene"].astype(str).to_numpy(); stim = p.obs["cell_type"].astype(str).to_numpy() == "stim"
    ctrl = stim & (lab == "non-targeting"); base = X[ctrl].mean(0)
    saved = np.load(pp)
    P, keeps, fitg, chkg = {}, {}, [], []
    for g in np.unique(lab[stim]):
        if g == "non-targeting":
            continue
        t = g.split("|")[0]; keep = np.ones(len(genes), bool)
        if t in gpos:
            keep[gpos[t]] = False
        P[g] = X[stim & (lab == g)].mean(0) - base; keeps[g] = keep
        if t in train and g in si:
            fitg.append(g)
        elif "state_transition__" + g in saved.files:
            chkg.append(g)
    # diagonal: a guide's prediction vs its own observed effect (fit) / its saved prediction (check);
    # off-diagonal baseline: vs the other guides' observed effects / saved predictions (shared response, wrong labels)
    fit = np.array([np.corrcoef(P[g][keeps[g]], S["eff"][si[g]][keeps[g]])[0, 1] for g in fitg])
    fit_off = np.array([np.mean([np.corrcoef(P[g][keeps[g]], S["eff"][si[h]][keeps[g]])[0, 1] for h in fitg[:50] if h != g]) for g in fitg])
    chk = np.array([np.corrcoef(P[g], saved["state_transition__" + g])[0, 1] for g in chkg])
    chk_off = np.array([np.mean([np.corrcoef(P[g], saved["state_transition__" + h])[0, 1] for h in chkg if h != g]) for g in chkg])
    out = dict(donor=DONOR, model="state_transition", n_fit_guides=int(len(fit)), fit_pearson=float(np.nanmean(fit)),
               fit_se=float(np.nanstd(fit, ddof=1) / np.sqrt(np.isfinite(fit).sum())), n_check=int(len(chk)),
               check_min_r=float(np.min(chk)) if len(chk) else None, check_mean_r=float(np.mean(chk)) if len(chk) else None,
               check_r=chk.tolist(), fit_offdiag_pearson=float(np.nanmean(fit_off)), check_offdiag_mean_r=float(np.mean(chk_off)) if len(chk) else None,
               rules=__doc__)
    json.dump(out, open(outp, "w"), indent=1)
    print(DONOR, "STATE fit on training guides r = %.3f (%.3f) n = %d; reproduction of saved test predictions r:" % (
        out["fit_pearson"], out["fit_se"], len(fit)), np.round(np.sort(chk), 3),
          "| off-diagonal: fit %.3f, check %.3f" % (out["fit_offdiag_pearson"], out["check_offdiag_mean_r"] or np.nan))


if __name__ == "__main__":
    main()
