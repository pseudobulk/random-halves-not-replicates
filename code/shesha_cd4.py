"""Shesha check (5090): split-half reproducibility of CD4 perturbations computed with shesha-geometry 0.2.29
(shesha.bio.split_half_reproducibility: cosine between the shift vectors of two random halves of a perturbation's
cells, each relative to the control centroid; 50 splits, base seed 320, at least 30 cells), and the same cosine with the
two halves formed by the perturbation's two sgRNAs, which share no clone (each sgRNA at least 15 cells; same control
centroid; the formula of shesha's _split_half_cosine applied to the fixed sgRNA halves). Input: shesha_extract_cd4.py.
Output: one JSON with summaries for unexpressed and expressed targets and per-target values.
Usage: shesha_cd4.py <cells.npz> <out.json>"""
import sys, json
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import shesha
from shesha.bio import split_half_reproducibility


def main():
    z = np.load(sys.argv[1], allow_pickle=False)
    X = csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    obs = pd.DataFrame({"perturbation": z["label"], "gp": z["gp"], "kind": z["kind"]})
    A = ad.AnnData(X=X, obs=obs); A.var_names = list(z["genes"])
    rnd = split_half_reproducibility(A, perturbation_key="perturbation", control_label="control", n_splits=50, random_state=320, min_cells=30)
    ctrl = np.asarray(X[(obs["perturbation"] == "control").values].mean(0)).ravel()
    lab = obs["perturbation"].values; gp = obs["gp"].values
    guide = {}
    for t in rnd.index:
        m = lab == t
        a, b = m & (gp == 0), m & (gp == 1)
        if a.sum() < 15 or b.sum() < 15:
            continue
        sa = np.asarray(X[a].mean(0)).ravel() - ctrl; sb = np.asarray(X[b].mean(0)).ravel() - ctrl
        guide[t] = float(sa @ sb / (np.linalg.norm(sa) * np.linalg.norm(sb)))
    kind = obs.drop_duplicates("perturbation").set_index("perturbation")["kind"]
    out = {"rules": __doc__, "shesha_version": shesha.__version__ if hasattr(shesha, "__version__") else "0.2.29", "summary": {}, "targets": {}}
    for k in ("unexpressed", "expressed"):
        ts = [t for t in rnd.index if kind[t] == k and t in guide]
        r = rnd.loc[ts, "split_half_cosine"].values; g = np.array([guide[t] for t in ts])
        out["summary"][k] = dict(n=len(ts), random_median=float(np.median(r)), guide_median=float(np.median(g)),
                                 random_mean=float(r.mean()), guide_mean=float(g.mean()), random_gt_guide=float((r > g).mean()),
                                 random_q=[float(x) for x in np.percentile(r, [25, 75])], guide_q=[float(x) for x in np.percentile(g, [25, 75])],
                                 median_cells=float(np.median(rnd.loc[ts, "n_cells"].values)))
        for t in ts:
            out["targets"][t] = dict(kind=k, random=float(rnd.loc[t, "split_half_cosine"]), guide=guide[t], n_cells=int(rnd.loc[t, "n_cells"]))
    print(json.dumps(out["summary"], indent=1))
    json.dump(out, open(sys.argv[2], "w"), indent=1)


if __name__ == "__main__":
    main()
