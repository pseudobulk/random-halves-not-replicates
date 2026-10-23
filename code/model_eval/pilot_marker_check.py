"""Model pilot: are the A-vs-B gaps carried by a few clone-marker genes? Rescores saved predictions with PuroR removed
(vector resistance gene; expression may follow the integration, so clone-specific) and with PuroR, TRBC1, TRBC2 and TRAC
removed (TCR constant regions mark clones by construction). Predictions: baselines (baseline_preds.npz) and the models
predicted through arc-state (adata_pred.h5ad, scored as pilot_state_score.py).
Usage: pilot_marker_check.py <data_dir> <out_json> name=<adata_pred.h5ad> ..."""
import sys, json
import numpy as np
from pilot_common import design, effects, score, STATES

SETS = {"all": [], "no_PuroR": ["PuroR"], "no_PuroR_TCR": ["PuroR", "TRBC1", "TRBC2", "TRAC"]}


def state_preds(pf, res):
    import anndata as ad
    p = ad.read_h5ad(pf); X = p.X.toarray() if hasattr(p.X, "toarray") else np.asarray(p.X)
    lab = p.obs["gene"].astype(str).to_numpy(); base = X[lab == "non-targeting"].mean(0)
    return {g: X[lab == g].mean(0) - base for g in np.unique(lab) if g != "non-targeting"}


def main():
    d, outp = sys.argv[1:3]
    des, res = effects(d); genes = res[STATES[1]]["genes"]
    present = {g: g in genes for s in SETS.values() for g in s}
    models = {}
    z = np.load(d.rstrip("/") + "/../out/baseline_preds.npz")
    for k in z.files:
        m, g = k.split("__", 1); models.setdefault(m, {})[g] = z[k]
    for a in sys.argv[3:]:
        name, pf = a.split("=", 1); models[name] = state_preds(pf, res)
    out = {"present": present, "models": {}}
    for m, pred in models.items():
        out["models"][m] = {}
        for sname, drop in SETS.items():
            mask = np.array([g not in drop for g in genes])
            s = score(pred, des, res, genes, gene_mask=mask)
            out["models"][m][sname] = {k: dict(A=round(v["pearson_A"], 4), B=round(v["pearson_B"], 4), leak=round(v["leak"], 4)) for k, v in s.items()}
        print(m, json.dumps(out["models"][m]), flush=True)
    json.dump(out, open(outp, "w"), indent=1)


if __name__ == "__main__":
    main()
