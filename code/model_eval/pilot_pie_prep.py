"""Model pilot: inputs for PIE (arc-pie 1.0.0) from the same cells and genes as the other pilot models.
One dataset 'cd4' with contexts rest / stim; perturbations are guide labels '<target>|g<gp>' so that each guide keeps its
own evidence (PIE pools a perturbation's train-row responses per perturbation name and removes only the query's own
row); an aliases file maps each guide label to its target gene for the knowledge sources (used only on a direct miss).
Writes, under <pie_dir>:
  counts/cd4.h5ad       raw counts (pilot genes), obs target_gene (guide label / non-targeting), cell_line
  expression/cd4.h5ad   log1p CP10K with the cell's UMI over all genes (as pilot_common), same obs
  splits/train.json     {"cd4.rest": all guides, "cd4.stim": training guides except val}
  splits/val.json       {"cd4.stim": the 40 validation guides of the STATE run (state/prep.json)}
  splits/test.json      {"cd4.stim": test guides}
  aliases.yaml          {source: {guide label: target}} for the per-perturbation sources
Usage: pilot_pie_prep.py <data_dir> <state_prep_json> <pie_dir>"""
import os, sys, json
import numpy as np, pandas as pd, anndata as ad
from scipy.sparse import vstack
from pilot_common import STATES, load, lognorm, design

PERT_SOURCES = ("esm2", "string_space", "depmap_gene_effect", "perturbation_text")


def main():
    d, sp, out = sys.argv[1:4]
    for sub in ("counts", "expression", "splits"):
        os.makedirs(os.path.join(out, sub), exist_ok=True)
    des = design(d); val = json.load(open(sp))["val"]
    Xc, Xl, obs = [], [], []
    for st, short in zip(STATES, ("rest", "stim")):
        D = load(d, st)
        Xc.append(D["X"].astype(np.float32)); Xl.append(lognorm(D["X"], D["umi"]))
        obs.append(pd.DataFrame({"target_gene": np.where(D["lab"] == "ctrl", "non-targeting", D["lab"]), "cell_line": short}))
        genes = D["genes"]
    o = pd.concat(obs, ignore_index=True); o.index = [str(i) for i in range(len(o))]
    var = pd.DataFrame(index=pd.Index(genes, name="gene_name")); var["gene_name"] = genes
    for X, sub in ((vstack(Xc).tocsr(), "counts"), (vstack(Xl).tocsr(), "expression")):
        ad.AnnData(X=X, obs=o.astype("category"), var=var).write_h5ad(os.path.join(out, sub, "cd4.h5ad"))
    guides = sorted(set(o["target_gene"]) - {"non-targeting"})
    test = sorted({"%s|g%d" % (t, p) for t in des["inert"] + des["strong"] for p in (0, 1)})
    train_stim = sorted(set(guides) - set(test) - set(val))
    json.dump({"cd4.rest": guides, "cd4.stim": train_stim}, open(os.path.join(out, "splits", "train.json"), "w"))
    json.dump({"cd4.stim": sorted(val)}, open(os.path.join(out, "splits", "val.json"), "w"))
    json.dump({"cd4.stim": test}, open(os.path.join(out, "splits", "test.json"), "w"))
    alias = {g: g.split("|")[0] for g in guides}
    with open(os.path.join(out, "aliases.yaml"), "w") as f:
        for s in PERT_SOURCES:
            f.write("%s:\n" % s)
            for g, t in alias.items():
                f.write('  "%s": "%s"\n' % (g, t))
    print("cells", len(o), "genes", len(genes), "guides", len(guides), "test", len(test), "val", len(val), "train stim", len(train_stim))


if __name__ == "__main__":
    main()
