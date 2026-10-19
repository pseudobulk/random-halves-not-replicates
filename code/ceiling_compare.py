"""How much does a random-half "ceiling" (cell-eval, Arc Institute) overstate what an independent replicate reaches?
(A100.) cell-eval --ceiling: each perturbation's cells and the control's are split into two disjoint random halves,
one half is scored as "real" and the other as "prediction", per-perturbation metrics are averaged and the
Spearman-Brown correction r' = 2r / (1 + r) maps the mean to full depth (README, _evaluator.py, checked 2026-10-18).
Here the same construction is computed twice on the same perturbations: halves = random cells (as cell-eval) and
halves = independent units, then compared.
  cd4    donor-1 single cells (Zhu / Marson 2025): random cells vs guides (guide rank parity; independent
         transductions). Rest and Stim48hr separately. Aggregates from cd4_sc_extract.py / cd4_sc_extract_guide.py.
  tahoe  48 lines x 91 conditions (plates 6 / 14, tahoe_plate_recompute.py): within-well halves vs replicate plates.
Metrics (cell-eval definitions on pseudobulk means of log1p CP10K; effect = perturbation mean - control mean, each
half against its own control half):
  pearson_delta        Pearson over genes between the two halves' effects, per perturbation
  discrimination_l1    1 - rank / n of the true perturbation when the "prediction" half's effect is matched by L1
                       distance against all "real" effects (target gene excluded); computed on random panels of
                       PANEL perturbations (5 panels, same panels for both designs), since cell-eval's all-pairs
                       loop is too large here; for Tahoe within each line across its 91 conditions
  top_overlap_N        PROXY for cell-eval's DE overlap_at_N (which needs a cell-level Wilcoxon test): overlap of the
                       top-N genes by |effect| between halves, N = 50, 100
Each metric: mean over perturbations, then Spearman-Brown if the mean is > 0 (as cell-eval). Subsets (cd4): all;
expressed targets (control >= 0.5); inert targets (< 0.01; nothing to predict beyond clone / guide structure);
strongest knockdowns (top 1000 expressed targets by the drop of the target gene itself, a criterion that does not
reuse the response genes).
Usage: ceiling_compare.py cd4 <cell_split_dir> <guide_split_dir> <label> <out_dir>
       ceiling_compare.py tahoe <D_arrays.npz> <out_dir>. Exploratory."""
import os, sys, json
import numpy as np

PANEL, NPANEL, NS = 500, 5, (50, 100)


def sb(r):
    return float(2 * r / (1 + r)) if r > 0 else float("nan")


def pearson_rows(A, B):
    a = A - A.mean(1, keepdims=True); b = B - B.mean(1, keepdims=True)
    return (a * b).sum(1) / np.sqrt((a ** 2).sum(1) * (b ** 2).sum(1))


def top_overlap(A, B, N):
    ta = np.argsort(-np.abs(A), 1)[:, :N]; tb = np.argsort(-np.abs(B), 1)[:, :N]
    return np.array([len(np.intersect1d(x, y)) / N for x, y in zip(ta, tb)])


def discrimination(real, pred, tcol=None):
    """real, pred: n x G effects of the same n perturbations. tcol: target gene column per row or -1."""
    D = np.empty((len(pred), len(real)))                             # pred i vs real j, in blocks (memory)
    for s0 in range(0, len(pred), 20):
        D[s0:s0 + 20] = np.abs(pred[s0:s0 + 20, None, :] - real[None, :, :]).sum(2)
    if tcol is not None:
        for i, t in enumerate(tcol):
            if t >= 0:
                D[i] -= np.abs(pred[i, t] - real[:, t])
    rank = (D < D[np.arange(len(D)), np.arange(len(D))][:, None]).sum(1)
    return 1 - rank / len(D)


def summarise(A, B, panels, tcol=None):
    out = {}
    pdl = pearson_rows(A, B)
    out["pearson_delta"] = dict(mean=float(np.nanmean(pdl)), ceiling=sb(np.nanmean(pdl)))
    for N in NS:
        ov = top_overlap(A, B, N); out["top_overlap_%d" % N] = dict(mean=float(ov.mean()), ceiling=sb(ov.mean()))
    if panels is not None:
        ds = np.concatenate([discrimination(A[p], B[p], None if tcol is None else tcol[p]) for p in panels])
        out["discrimination_l1"] = dict(mean=float(ds.mean()), ceiling=sb(ds.mean()))
    return out


def run_cd4(dc, dg, label, out):
    from cd4_readout import load, MIN_HALF, MIN_EXPR
    L = {"cell": load(dc), "guide": load(dg)}
    genes = L["cell"][1]; assert list(genes) == list(L["guide"][1])
    E, ctrlm = {}, None
    for nm, (groups, _, S, n) in L.items():
        ci = groups.index("__ctrl__"); ctrl = np.stack([np.asarray(S[ci, h]) / n[ci, h] for h in range(2)])
        if nm == "cell":
            ctrlm = ctrl.mean(0)
        E[nm] = (groups, S, n, ctrl)
    gk = np.where(ctrlm >= MIN_EXPR)[0]; gpos = {genes[j]: k for k, j in enumerate(gk)}; gall = {g: j for j, g in enumerate(genes)}
    ok = [set(t for i, t in enumerate(E[nm][0]) if t != "__ctrl__" and E[nm][2][i].min() >= MIN_HALF) for nm in E]
    perts = sorted(ok[0] & ok[1])
    H, ncell = {}, {}
    for nm, (groups, S, n, ctrl) in E.items():
        gi = {g: i for i, g in enumerate(groups)}; idx = np.array([gi[t] for t in perts])
        H[nm] = np.stack([np.asarray(S[idx, h])[:, gk] / n[idx, h][:, None] - ctrl[h, gk] for h in range(2)])
        ncell[nm] = np.median(n[idx], 0).tolist()
    texp = np.array([ctrlm[gall[t]] if t in gall else 0.0 for t in perts])
    kd = np.array([(H["cell"][0][k, gpos[t]] + H["cell"][1][k, gpos[t]]) / 2 if t in gpos else 0.0 for k, t in enumerate(perts)])
    tcol = np.array([gpos.get(t, -1) for t in perts])
    measured = np.array([t in gall for t in perts])          # targets outside the probe panel are not "inert" (10-18 fix)
    expd, inert = (texp >= 0.5) & measured, (texp < 0.01) & measured
    strong = np.zeros(len(perts), bool); ce = np.where(expd)[0]; strong[ce[np.argsort(kd[ce])[:1000]]] = True
    res = {"rules": __doc__, "label": label, "n_perturbations": len(perts), "n_genes": int(len(gk)), "median_cells_per_half": ncell,
           "n_expressed": int(expd.sum()), "n_inert": int(inert.sum())}
    rng = np.random.default_rng(1)
    for sname, m in (("all", np.ones(len(perts), bool)), ("expressed", expd), ("inert", inert), ("strong_kd", strong)):
        ix = np.where(m)[0]
        panels = [np.sort(rng.choice(ix, min(PANEL, len(ix)), replace=False)) for _ in range(NPANEL)]
        r = {"n": int(len(ix))}
        for nm in H:
            r[nm] = summarise(H[nm][0][ix], H[nm][1][ix], [np.searchsorted(ix, p) for p in panels], tcol[ix])
        res[sname] = r
        print(label, sname, json.dumps(r), flush=True)
    return res


def run_tahoe(npz, out):
    z = np.load(npz); X = {k: z[k].astype(np.float64) for k in ("plate6_0", "plate6_1", "plate14_0", "plate14_1")}
    C, P, G = X["plate6_0"].shape
    designs = {"within_well_plate6": ("plate6_0", "plate6_1"), "within_well_plate14": ("plate14_0", "plate14_1"),
               "across_plates": ("plate6_0", "plate14_1"), "across_plates_rev": ("plate14_0", "plate6_1")}
    res = {"rules": __doc__, "n_lines": C, "n_conditions": P, "n_genes": G}
    for nm, (a, b) in designs.items():
        A, B = X[a].reshape(C * P, G), X[b].reshape(C * P, G)
        r = summarise(A, B, None)
        ds = np.concatenate([discrimination(X[a][c], X[b][c]) for c in range(C)])
        r["discrimination_l1_within_line"] = dict(mean=float(ds.mean()), ceiling=sb(ds.mean()))
        res[nm] = r
        print(nm, json.dumps(r), flush=True)
    return res


def main():
    mode = sys.argv[1]
    if mode == "cd4":
        out = sys.argv[5]; os.makedirs(out, exist_ok=True); res = run_cd4(*sys.argv[2:6])
    else:
        out = sys.argv[3]; os.makedirs(out, exist_ok=True); res = run_tahoe(*sys.argv[2:4])
    json.dump(res, open(out + "/metrics_ceiling_compare.json", "w"), indent=1)
    print("saved", flush=True)


if __name__ == "__main__":
    main()
