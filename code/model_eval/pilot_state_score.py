"""Model pilot: score STATE's predictions (adata_pred.h5ad from `state tx predict --predict-only`) like the other models.
A test guide's predicted effect is the mean predicted expression over its cells minus the mean over the predicted
non-targeting cells of the stimulated state when the prediction file has them; otherwise minus the observed stimulated
control mean (recorded in the output as baseline). Usage: pilot_state_score.py <data_dir> <adata_pred> <out_json>"""
import sys, json
import numpy as np, anndata as ad
from pilot_common import effects, score, STATES


def main():
    d, pf, outp = sys.argv[1:4]
    des, res = effects(d)
    p = ad.read_h5ad(pf)
    X = p.X.toarray() if hasattr(p.X, "toarray") else np.asarray(p.X)
    genes = res[STATES[1]]["genes"]
    assert X.shape[1] == len(genes), "prediction width differs from the pilot genes"   # X_hvg columns = pilot genes in order
    lab = p.obs["gene"].astype(str).to_numpy()
    ct = p.obs["cell_type"].astype(str).to_numpy() if "cell_type" in p.obs else np.full(len(lab), "stim")
    stim = ct == "stim"
    ctrl = stim & (lab == "non-targeting")
    if ctrl.sum() > 0:
        base, how = X[ctrl].mean(0), "predicted non-targeting cells (n=%d)" % ctrl.sum()
    else:
        base, how = res[STATES[1]]["ctrl"], "observed stimulated control mean"
    pred = {g: X[stim & (lab == g)].mean(0) - base for g in np.unique(lab[stim]) if g != "non-targeting"}
    sc = score(pred, des, res, genes)
    json.dump(dict(model="STATE", baseline=how, n_pred_guides=len(pred), scores=sc, rules=__doc__), open(outp, "w"), indent=1)
    print(how, len(pred), json.dumps(sc))


if __name__ == "__main__":
    main()
