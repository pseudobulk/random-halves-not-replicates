"""Model pilot: State Embedding (SE-600M) of every pilot cell, for State Transition on SE embeddings (the full State
recipe). Per state, the all-genes raw counts (pilot_<state>_allgenes.npz, the same cells in the same order as the
pilot subset; checked) are written to a temporary h5ad and embedded with `state emb transform`; the embeddings are kept
as emb_<state>.npy and the temporary h5ad is removed (the 5090 disk is nearly full).
Then <out>/data/cd4_d1_se.h5ad is written with the cells in pilot_common.build_adata order: obsm X_state (SE), obsm
X_hvg (log1p CP10K over the pilot genes, the gene-space targets), X = X_hvg (sparse), obs as pilot_state_prep.py.
Usage: CUDA_VISIBLE_DEVICES=<gpu> pilot_se_embed.py <data_dir> <se_dir> <out_dir> [batch_size]"""
import os, sys, json, subprocess
import numpy as np, pandas as pd, anndata as ad
from scipy.sparse import csr_matrix
from pilot_common import STATES, build_adata

STATE_BIN = os.environ.get("STATE_BIN", "state")   # arc-state command line (STATE_BIN in config.env)


def main():
    d, se, out = sys.argv[1:4]; bs = sys.argv[4] if len(sys.argv) > 4 else "128"
    os.makedirs(out + "/data", exist_ok=True); embs = []
    for st in STATES:
        npy = os.path.join(out, "emb_%s.npy" % st)
        z = np.load(os.path.join(d, "pilot_%s_allgenes.npz" % st), allow_pickle=False)
        p = np.load(os.path.join(d, "pilot_%s.npz" % st), allow_pickle=False)
        for k in ("target", "gp", "kind"):
            assert (z[k] == p[k]).all(), "cell order differs between the all-genes and pilot files (%s, %s)" % (st, k)
        if not os.path.exists(npy):
            X = csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
            tmp = os.path.join(out, "tmp_%s.h5ad" % st)
            ad.AnnData(X=X, obs=pd.DataFrame(index=["c%d" % i for i in range(X.shape[0])]),
                       var=pd.DataFrame(index=[str(g) for g in z["genes"]])).write_h5ad(tmp)
            del X
            from pilot_se_infer import embed   # safetensors weights; `state emb transform` needs the 11.5 GB .ckpt
            np.save(npy, np.asarray(embed(se, tmp, int(bs))))
            os.remove(tmp)
        e = np.load(npy); assert e.shape[0] == len(z["target"]); embs.append(e.astype(np.float32))
        print(st, "embedding", e.shape, flush=True)
    des, a, genes = build_adata(d)
    E = np.vstack(embs); assert E.shape[0] == a.n_obs
    obs = pd.DataFrame({"gene": np.where(a.obs["condition"] == "ctrl", "non-targeting", a.obs["condition"]),
                        "cell_type": a.obs["state"].astype(str), "gem_group": a.obs["state"].astype(str)}, index=a.obs_names).astype("category")
    b = ad.AnnData(X=csr_matrix(a.X), obs=obs, var=a.var)
    b.obsm["X_hvg"] = a.X; b.obsm["X_state"] = E
    b.write_h5ad(os.path.join(out, "data", "cd4_d1_se.h5ad"))
    json.dump(dict(n_cells=int(b.n_obs), se_dim=int(E.shape[1]), n_genes=len(genes)), open(os.path.join(out, "se_prep.json"), "w"))
    print("wrote", b.shape, "X_state", E.shape, flush=True)


if __name__ == "__main__":
    main()
