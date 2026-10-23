"""Model pilot: Stack (arc-stack, checkpoint Stack-Large-Aligned) in-context generation, no training on these data.
As the Stack prediction tutorial (notebooks/tutorial-predict.ipynb; prompt 0.25, context 0.4, 5 steps, batch 16):
  base   the resting cells of each test guide (one split per guide label) and resting control cells ('non-targeting')
  test   NTEST stimulated control cells (random, seed 7)
Raw counts over all genes (pilot_<state>_allgenes.npz; Stack applies its own library-size scaling and log1p). For each
split Stack generates the test cells under that condition. Generated cells are put on log1p CP10K (their own totals);
a guide's predicted effect is the mean over its generated cells minus the mean over the generated control cells, on the
pilot genes, scored by pilot_common.score on the pilot genes present in Stack's gene list.
Usage: CUDA_VISIBLE_DEVICES=<gpu> pilot_stack.py <data_dir> <model_dir> <out_dir> [n_test]"""
import os, sys, json, time, logging
import numpy as np, pandas as pd, anndata as ad
from scipy.sparse import csr_matrix, vstack
from pilot_common import design, effects, score, STATES, test_pairs

logging.basicConfig(level=logging.INFO)


def load_all(d, st):
    z = np.load(os.path.join(d, "pilot_%s_allgenes.npz" % st), allow_pickle=False)
    X = csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    lab = np.array(["non-targeting" if k == "control" else "%s|g%d" % (t, p) for t, p, k in zip(z["target"], z["gp"], z["kind"])])
    return X, lab, [str(g) for g in z["genes"]]


def main():
    d, md, out = sys.argv[1:4]; ntest = int(sys.argv[4]) if len(sys.argv) > 4 else 1000
    os.makedirs(out, exist_ok=True); t0 = time.time()
    from stack.cli.generation import generate
    from stack.data.training.datasets import load_gene_list
    des = design(d)
    tests = [g for _, _, g, _ in test_pairs(des)]
    Xr, lr, genes = load_all(d, STATES[0]); Xs, ls, genes_s = load_all(d, STATES[1]); assert genes == genes_s
    base_rows = np.isin(lr, tests + ["non-targeting"])
    base = ad.AnnData(X=Xr[base_rows], obs=pd.DataFrame({"condition": lr[base_rows]}, index=["r%d" % i for i in np.where(base_rows)[0]]),
                      var=pd.DataFrame(index=genes))
    sc = np.where(ls == "non-targeting")[0]; sc = np.random.default_rng(7).choice(sc, min(ntest, len(sc)), replace=False)
    test = ad.AnnData(X=Xs[sc], obs=pd.DataFrame({"condition": "non-targeting"}, index=["s%d" % i for i in sc]), var=pd.DataFrame(index=genes))
    base.write_h5ad(os.path.join(out, "base.h5ad")); test.write_h5ad(os.path.join(out, "test.h5ad"))
    gl = os.path.join(md, "basecount_1000per_15000max.pkl"); ck = [f for f in os.listdir(md) if f.endswith(".ckpt")][0]
    stack_genes = [str(g).upper() for g in load_gene_list(gl)]
    _, res = effects(d); pgenes = res[STATES[1]]["genes"]
    col = {g: j for j, g in enumerate(stack_genes)}
    pidx = np.array([col.get(g.upper(), -1) for g in pgenes]); mask = pidx >= 0
    print("stack genes", len(stack_genes), "pilot genes covered", int(mask.sum()), "of", len(pgenes), "%.0fs" % (time.time() - t0), flush=True)
    means = {}

    def keep(split, a):
        X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
        X = np.log1p(X / np.maximum(X.sum(1, keepdims=True), 1e-8) * 1e4)
        m = np.zeros(len(pgenes), np.float32); m[mask] = X[:, pidx[mask]].mean(0); means[split] = m
        if len(means) % 25 == 0:
            print("generated", len(means), "%.0fs" % (time.time() - t0), flush=True)

    generate(os.path.join(md, ck), os.path.join(out, "base.h5ad"), os.path.join(out, "test.h5ad"), gl,
             split_column="condition", prompt_ratio=0.25, context_ratio=0.4, num_steps=5, batch_size=16,
             random_seed=0, on_generation=keep)
    base_m = means["non-targeting"]
    pred = {g: means[g] - base_m for g in tests if g in means}
    s = score(pred, des, res, pgenes, gene_mask=mask)
    json.dump(dict(model="Stack", checkpoint=ck, n_test_cells=int(len(sc)), pilot_genes_covered=int(mask.sum()), scores=s,
                   seconds=time.time() - t0, rules=__doc__), open(os.path.join(out, "stack_scores.json"), "w"), indent=1)
    print(json.dumps(s), "%.0fs" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
