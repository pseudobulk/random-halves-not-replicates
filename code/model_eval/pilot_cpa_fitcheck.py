"""Model pilot: fit check for the trained CPA (as pilot_biolord's): the test guides' resting cells were in training;
decode resting control cells under each test guide's label (state rest) and correlate the predicted effect with the
guide's observed resting effect. A model that cannot reproduce these has not learned per-guide effects.
Usage: CUDA_VISIBLE_DEVICES=<gpu> pilot_cpa_fitcheck.py <data_dir> <cpa_out_dir>"""
import os, sys, json
import numpy as np
import torch
_load = torch.load
torch.load = lambda *a, **k: _load(*a, **{**k, "weights_only": False})   # our own checkpoint; torch >= 2.6 default refuses it
import cpa
from pilot_cpa import build
from pilot_common import effects


def main():
    d, out = sys.argv[1:3]
    des, adata, genes = build(d)
    cpa.CPA.setup_anndata(adata, perturbation_key="condition", control_group="ctrl", dosage_key=None,
                          categorical_covariate_keys=["state"], is_count_data=False, max_comb_len=1)
    model = cpa.CPA.load(os.path.join(out, "model"), adata=adata, use_gpu=True)
    obs = adata.obs; lab = obs["condition"].to_numpy(); rest = (obs["state"] == "rest").to_numpy()
    tests = sorted(set(lab[(obs["split"] == "ood").to_numpy()]))
    rows = rest & np.isin(lab, tests)
    rctrl = rest & (lab == "ctrl")
    rng = np.random.default_rng(44)
    a = adata.copy(); Xn = a.X.copy(); Xn[rows] = adata.X[rng.choice(np.where(rctrl)[0], rows.sum(), replace=True)]; a.X = Xn
    model.predict(a, batch_size=2048)
    P = np.asarray(a.obsm["CPA_pred"]); base = P[rctrl].mean(0)
    _, res = effects(d); R = res["D1_Rest"]; ri = {l: i for i, l in enumerate(R["labels"])}
    fit = [np.corrcoef(P[rows & (lab == g)].mean(0) - base, R["eff"][ri[g]])[0, 1] for g in tests if g in ri]
    fc = dict(resting_fit_pearson=float(np.nanmean(fit)), n=len(fit))
    s = json.load(open(os.path.join(out, "cpa_scores.json"))); s["fit_check"] = fc
    json.dump(s, open(os.path.join(out, "cpa_scores.json"), "w"), indent=1)
    print("fit check", fc, flush=True)


if __name__ == "__main__":
    main()
