"""cell-eval itself (0.8.2, PyPI 2026-07-27): random-half ceiling vs guide-half "ceiling" on CD4 donor-1 Rest
(preprint review 2.5; runs in its own environment with cell-eval 0.8.2, see README). Input: cd4_celleval_subset.py output
(150 inert + 150 strong, knockdown-matched targets, <= 400 cells each, 15,000 non-targeting controls; raw counts and
per-cell UMI over all genes). Expression = log1p CP10K.
  random  MetricsEvaluator(adata_pred=None, ...).compute_ceiling(profile="full", seed) for seeds 0, 1, 2: cell-eval's
          own disjoint random halves of every perturbation and of the controls, metrics averaged over perturbations,
          Spearman-Brown on the reliability metrics (its _spearman_brown_correct).
  guide   real = cells of guide parity 0 (controls: non-targeting guides of parity 0), pred = parity 1, and reversed;
          MetricsEvaluator(...).compute(profile="full"); the per-perturbation results are averaged and passed through
          the same _spearman_brown_correct. Guides are independent transductions (no shared clones).
For both, the same averaging + correction is also applied to the inert and the strong targets separately.
10-18 audit fixes: equal numbers of guide-0 and guide-1 cells per perturbation (same cells for both designs); check
that the re-aggregation equals cell-eval's own aggregate; final summary on the perturbations shared by all runs.
Usage: run_celleval_ceiling.py <subset.npz> <out_dir> [threads]. Exploratory."""
import os, sys, json
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import polars as pl
from cell_eval import MetricsEvaluator
from cell_eval._evaluator import _spearman_brown_correct, SB_METRICS


def build(npz):
    z = np.load(npz, allow_pickle=False)
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    X = sp.csr_matrix(X.multiply((1e4 / z["umi"])[:, None])); X.data = np.log1p(X.data); X = X.astype(np.float32)
    obs = pd.DataFrame({"target_gene": z["target"].astype(str), "gp": z["gp"].astype(int), "kind": z["kind"].astype(str)})
    obs.index = obs.index.astype(str)
    A = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=z["genes"].astype(str)))
    kinds = obs.drop_duplicates("target_gene").set_index("target_gene")["kind"].to_dict()
    return A, kinds


def summarise(res, kinds):
    """res: per-perturbation polars frame (column 'perturbation'); mean over perturbations then cell-eval's SB."""
    out = {}
    pert = res["perturbation"].to_list()
    for sname in ("all", "inert", "strong"):
        keep = [p for p in pert if sname == "all" or kinds.get(p) == sname]
        sub = res.filter(pl.col("perturbation").is_in(keep))
        raw = sub.drop("perturbation").mean()
        sbv = _spearman_brown_correct(raw)
        out[sname] = {"n": len(keep), "raw_mean": {c: raw[c][0] for c in raw.columns},
                      "ceiling": {c: sbv[c][0] for c in sbv.columns if c in SB_METRICS}}
    return out


def main():
    npz, out = sys.argv[1:3]; th = int(sys.argv[3]) if len(sys.argv) > 3 else 48
    os.makedirs(out, exist_ok=True)
    A, kinds = build(npz)
    if os.environ.get("SMOKE"):          # quick API / timing check on a small subset
        rng = np.random.default_rng(0)
        keep = [p for k in ("inert", "strong") for p in rng.choice([q for q, v in kinds.items() if v == k], 10, replace=False)]
        ctrl = rng.choice(np.where(A.obs["target_gene"] == "non-targeting")[0], 2000, replace=False)
        m = A.obs["target_gene"].isin(keep).values; m[ctrl] = True
        A = A[m].copy(); kinds = {p: kinds[p] for p in keep}; kinds["non-targeting"] = "control"
    # 10-18 audit fix: guide halves had unequal cell numbers (random halves are equal by construction), which lowers
    # the guide-side scores for reasons unrelated to clones. Downsample every perturbation (and the controls) to equal
    # numbers of guide-0 and guide-1 cells; both designs then use exactly the same cells.
    rng = np.random.default_rng(1); keep = np.zeros(A.n_obs, bool); tg = A.obs["target_gene"].values; gp = A.obs["gp"].values
    for t in np.unique(tg):
        i0, i1 = np.where((tg == t) & (gp == 0))[0], np.where((tg == t) & (gp == 1))[0]
        m = min(len(i0), len(i1))
        if m >= 10:
            keep[rng.choice(i0, m, replace=False)] = True; keep[rng.choice(i1, m, replace=False)] = True
    A = A[keep].copy(); kinds = {t: kinds[t] for t in np.unique(A.obs["target_gene"].values)}
    ncell = A.obs.groupby(["target_gene", "gp"]).size().unstack()
    print("cells", A.n_obs, "genes", A.n_vars, "perturbations", len(kinds), "median cells per guide side", float(ncell.median().median()), flush=True)
    res = {"rules": __doc__, "n_cells": int(A.n_obs), "n_genes": int(A.n_vars)}; raw = {}
    for seed in ((0,) if os.environ.get("SMOKE") else (0, 1, 2)):
        ev = MetricsEvaluator(adata_pred=None, adata_real=A, control_pert="non-targeting", pert_col="target_gene",
                              num_threads=th, outdir="%s/random_seed%d" % (out, seed))
        r, agg = ev.compute_ceiling(profile="full", seed=seed)
        res["random_seed%d" % seed] = summarise(r, kinds)
        # check: our re-aggregation reproduces cell-eval's own aggregate (no double correction)
        mine = res["random_seed%d" % seed]["all"]["ceiling"]
        diffs = [abs(mine[c] - agg[c][0]) for c in mine if c in agg.columns and mine[c] == mine[c] and agg[c][0] is not None and agg[c][0] == agg[c][0]]
        res["random_seed%d" % seed]["check_max_abs_diff_vs_celleval_agg"] = float(max(diffs)) if diffs else None
        raw[("random", seed)] = r
        print("random seed", seed, json.dumps({k: v["ceiling"] for k, v in res["random_seed%d" % seed].items() if isinstance(v, dict)}, default=float),
              "check", res["random_seed%d" % seed]["check_max_abs_diff_vs_celleval_agg"], flush=True)
    for nm, (ra, pa) in (("guide_0real_1pred", (0, 1)), ("guide_1real_0pred", (1, 0))):
        real, pred = A[A.obs["gp"] == ra].copy(), A[A.obs["gp"] == pa].copy()
        ev = MetricsEvaluator(adata_pred=pred, adata_real=real, control_pert="non-targeting", pert_col="target_gene",
                              num_threads=th, outdir="%s/%s" % (out, nm))
        r, agg = ev.compute(profile="full")
        res[nm] = summarise(r, kinds); raw[(nm, 0)] = r
        print(nm, json.dumps({k: v["ceiling"] for k, v in res[nm].items()}, default=float), flush=True)
    # both designs restricted to the perturbations present in every run
    shared = set.intersection(*[set(r["perturbation"].to_list()) for r in raw.values()])
    res["n_shared_perturbations"] = len(shared)
    res["shared"] = {"%s_%d" % k: summarise(r.filter(pl.col("perturbation").is_in(list(shared))), kinds) for k, r in raw.items()}
    json.dump(res, open(out + "/metrics_celleval_ceiling.json", "w"), indent=1, default=float)
    print("saved", flush=True)


if __name__ == "__main__":
    main()
