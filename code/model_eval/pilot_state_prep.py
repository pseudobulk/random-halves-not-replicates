"""Model pilot: input for STATE (State Transition, arc-state) from the same cells and genes as pilot_common.build_adata.
h5ad: X and obsm['X_hvg'] = log1p CP10K over the pilot genes; obs gene = guide label ('non-targeting' for controls),
cell_type = rest / stim, gem_group = state (control cells are matched within a state). TOML: one dataset 'cd4',
fewshot on cd4.stim with test = the test guides (their stimulated cells held out; their resting cells train) and val =
NVAL training guides at random (seed 1), whose stimulated cells are then used for validation only.
Usage: pilot_state_prep.py <data_dir> <state_dir>"""
import os, sys, json
import numpy as np
from pilot_common import build_adata, test_targets

NVAL = 40


def main():
    d, sd = sys.argv[1:3]
    os.makedirs(sd + "/data", exist_ok=True)
    des, a, genes = build_adata(d)
    a.obs["gene"] = np.where(a.obs["condition"] == "ctrl", "non-targeting", a.obs["condition"])
    a.obs["cell_type"] = a.obs["state"].astype(str)
    a.obs["gem_group"] = a.obs["state"].astype(str)
    a.obsm["X_hvg"] = a.X.copy()
    a.obs = a.obs[["gene", "cell_type", "gem_group"]].astype("category")
    a.write_h5ad(sd + "/data/cd4_d1.h5ad")
    held = sorted({g for t in test_targets(des) for g in ("%s|g0" % t, "%s|g1" % t)})
    trains = sorted({"%s|g%d" % (t, p) for t in des["train"] for p in (0, 1)})
    val = sorted(np.random.default_rng(1).choice(trains, NVAL, replace=False).tolist())
    q = lambda xs: "[" + ", ".join('"%s"' % x for x in xs) + "]"
    toml = '[datasets]\ncd4 = "%s/data/"\n\n[training]\ncd4 = "train"\n\n[fewshot]\n[fewshot."cd4.stim"]\nval = %s\ntest = %s\n' % (sd, q(val), q(held))
    open(sd + "/cd4_fewshot.toml", "w").write(toml)
    json.dump(dict(n_cells=int(a.n_obs), n_genes=len(genes), n_test=len(held), val=val), open(sd + "/prep.json", "w"), indent=1)
    print("cells", a.n_obs, "genes", len(genes), "test guides", len(held), "val guides", len(val))


if __name__ == "__main__":
    main()
