"""Model evaluation: one npz of predicted stimulated effects per donor (the four baselines from baseline_preds.npz and
STATE Transition from adata_pred.h5ad, relative to its predicted non-targeting stimulated cells, as me_score_detail.py),
for the cross-donor control (me_cross_donor.py). Keys '<model>__<guide label>', plus genes.
Usage (env PILOT_DONOR): me_export_preds.py <data_dir> <baseline_preds.npz> <adata_pred.h5ad> <out_npz>"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "4")
import sys
import numpy as np
from pilot_common import STATES, effects


def main():
    d, bp, sp, outp = sys.argv[1:5]
    des, res = effects(d); genes = res[STATES[1]]["genes"]
    out = dict(np.load(bp))
    import anndata as ad
    p = ad.read_h5ad(sp)
    X = p.X.toarray() if hasattr(p.X, "toarray") else np.asarray(p.X)
    lab = p.obs["gene"].astype(str).to_numpy(); stim = p.obs["cell_type"].astype(str).to_numpy() == "stim"
    base = X[stim & (lab == "non-targeting")].mean(0)
    for g in np.unique(lab[stim]):
        if g != "non-targeting":
            out["state_transition__" + g] = X[stim & (lab == g)].mean(0) - base
    np.savez_compressed(outp, genes=np.array(genes), **out)
    print(outp, len(out))


if __name__ == "__main__":
    main()
