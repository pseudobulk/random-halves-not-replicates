"""Tahoe key numbers recomputed with independent wells (A100; follows Z1 ah). On the 93 conditions shared by plate 6 and
its biological replicate plate 14, the two "halves" are plate 6 (both cell halves pooled) and plate 14 (pooled), i.e.
separate wells, spheroids and DMSO references. For comparison the same quantities with the old design (random cell
halves within the plate-6 wells) on the same conditions and lines.
  gamma share of the specific signal (interaction over lines x conditions vs conditions' shared part), as tahoe_m2
  readout share: half-B gamma predicted from half-A beta by a per-line linear readout (top-K PCs of beta over the
                 training conditions, ridge, 2 folds over conditions)
  diagnostics: raw cross products of gamma and beta for within-well halves (plates 6 and 14) and cross-plate halves
                 (equal cell numbers), so that loss of shared signal can be told apart from added independent noise
  response strength (projection on the condition's shared direction): split reproducibility across lines, and LOO R2
                 from 5 HPA (monolayer) baseline PCs
Usage: tahoe_plate_recompute.py <tahoe_agg dir> <hpa rna_celline.tsv> <out_dir>. Exploratory."""
import os, sys, json, ast, csv
import numpy as np
from scipy.stats import spearmanr
import pyarrow.parquet as pq
import config
from tahoe_arms import kr_loo
from tahoe_hpa_check import norm

MIN_HALF, K = 10, 10


def share(g, gA, F):
    sig = (g ** 2).sum() - ((gA - g) ** 2).sum() / 2
    return float(((g ** 2).sum() - ((g - F) ** 2).sum()) / sig)


def ridge(Z, Y, lam):
    k = Z.shape[1]
    return np.linalg.solve(Z.T @ Z + lam * np.trace(Z.T @ Z) / k * np.eye(k), Z.T @ Y)


def readout_share(XA, XB, seed=3, two_way=False):
    """XA, XB: lines x conds x genes. one-way: gamma = X - mean over lines (keeps the line main effect alpha, as the
    readout target of tahoe_readout.py); two-way: gamma = interaction (alpha removed, as tahoe_m2.py), beta = condition
    main effect. alpha includes the error of the line's DMSO reference, which is constant across conditions; within-well
    halves share their DMSO wells, so that error is reproducible there. Per-line readout from beta."""
    L, P, G = XA.shape
    if two_way:
        mA, mB = XA.mean((0, 1)), XB.mean((0, 1))
        aA, aB = XA.mean(1) - mA, XB.mean(1) - mB
        bA, bB = XA.mean(0) - mA, XB.mean(0) - mB
        gA, gB = XA - aA[:, None] - bA - mA, XB - aB[:, None] - bB - mB
    else:
        bA, bB = XA.mean(0), XB.mean(0); gA, gB = XA - bA, XB - bB
    perm = np.random.default_rng(seed).permutation(P); folds = (perm[: P // 2], perm[P // 2:])
    F = np.zeros_like(gB)
    for f in (0, 1):
        tr, te = folds[f], folds[1 - f]
        V = np.linalg.svd(bA[tr], full_matrices=False)[2][:K].T; Z = bA @ V
        for c in range(L):
            F[c, te] = Z[te] @ ridge(Z[tr], gB[c, tr], 1.0)
    gs = (gA * gB).sum(); bs = (bA * bB).sum() * L
    # raw cross products (unbiased for reproducible energy whatever the independent noise) and own energies; readout
    # share also with the cross product as denominator, which does not assume equal noise in A and B
    expl = (gB ** 2).sum() - ((gB - F) ** 2).sum()
    return dict(gamma_share=float(gs / (gs + bs)), readout_share=share(gB, gA, F), readout_share_xden=float(expl / gs),
                gamma_cross=float(gs), beta_cross=float(bs), gammaA_energy=float((gA ** 2).sum()),
                gammaB_energy=float((gB ** 2).sum()), **({"alpha_cross": float((aA * aB).sum() * P)} if two_way else {}))


def main():
    agg, hpa, out = sys.argv[1:4]; os.makedirs(out, exist_ok=True)
    idx = json.load(open(agg + "/index.json")); lines, samples, genes = idx["lines"], idx["samples"], idx["genes"]
    C, Ps = len(lines), len(samples)
    S = np.load(agg + "/sums64.npy", mmap_mode="r"); cnt = np.load(agg + "/counts.npy"); n = cnt[:, 0]
    sm = pq.read_table(config.tahoe_meta("sample_metadata.parquet")).to_pydict()
    info = {s: (d.strip(), str(ast.literal_eval(dc)[0][1]), p) for s, d, dc, p in zip(sm["sample"], sm["drug"], sm["drugname_drugconc"], sm["plate"])}
    row = lambda c, s, h: (c * Ps + s) * 2 + h
    sidx = {p: [i for i, s in enumerate(samples) if info[s][2] == p] for p in ("plate6", "plate14")}
    cond = {p: {(info[samples[i]][0], info[samples[i]][1]): i for i in sidx[p]} for p in sidx}
    isd = lambda k: k[0].startswith("DMSO")
    shared = sorted(k for k in cond["plate6"] if k in cond["plate14"] and not isd(k))
    dmso = {p: [i for i in sidx[p] if info[samples[i]][0].startswith("DMSO")] for p in sidx}
    rr = [row(c, s, h) for c in range(C) for s in dmso["plate6"] for h in range(2)]
    gm = np.asarray(S[rr].sum(0)) / n[rr].sum(); gk = np.where(gm >= 0.1)[0]
    def mean_rows(r):
        m = n[r].sum(); return np.asarray(S[r][:, gk].sum(0)) / m if m > 0 else None
    # lines with enough cells in every shared condition on both plates and both halves
    avail = np.array([[all(n[row(c, cond[p][k], h)] >= MIN_HALF for p in sidx for h in range(2)) for k in shared] for c in range(C)])
    okl = [c for c in range(C) if avail[c].mean() >= 0.85]
    keepk = [j for j in range(len(shared)) if avail[okl, j].all()]
    print("lines available in >= 85% of conditions", len(okl), "conditions with all of them", len(keepk), flush=True)
    shared = [shared[j] for j in keepk]
    print("conditions", len(shared), "lines", len(okl), "genes", len(gk), flush=True)
    D = {}
    for p in sidx:
        for h in (0, 1, "pool"):
            hs = [0, 1] if h == "pool" else [h]
            ref = {c: mean_rows([row(c, s, hh) for s in dmso[p] for hh in hs]) for c in okl}
            D[(p, h)] = np.array([[mean_rows([row(c, cond[p][k], hh) for hh in hs]) - ref[c] for k in shared] for c in okl])
    res = {"rules": __doc__, "n_conditions": len(shared), "n_lines": len(okl), "lines": [lines[c] for c in okl]}
    np.savez(out + "/D_arrays.npz", **{"%s_%s" % k: v.astype(np.float32) for k, v in D.items()},
             lines=np.array([lines[c] for c in okl]), conditions=np.array(["%s|%s" % k for k in shared]), genes=np.array([genes[j] for j in gk]))
    designs = {"old_within_well": (("plate6", 0), ("plate6", 1)), "within_well_plate14": (("plate14", 0), ("plate14", 1)),
               "cross_plate_halves": (("plate6", 0), ("plate14", 1)), "cross_plate_halves_rev": (("plate14", 1), ("plate6", 0)),
               "new_across_plates": (("plate6", "pool"), ("plate14", "pool")), "new_across_plates_rev": (("plate14", "pool"), ("plate6", "pool"))}
    for kk, (a, b) in designs.items():
        res[kk] = readout_share(D[a], D[b]); res[kk + "_twoway"] = readout_share(D[a], D[b], two_way=True)
        print(kk, json.dumps(res[kk]), flush=True); print(kk + "_twoway", json.dumps(res[kk + "_twoway"]), flush=True)
    # response strength and HPA predictability
    d = pq.read_table(config.tahoe_meta("cell_line_metadata.parquet")).to_pydict(); want = {norm(nm): c for c, nm in zip(d["Cell_ID_Cellosaur"], d["cell_name"]) if nm}
    gsel = [genes[j] for j in gk]; gp = {g: j for j, g in enumerate(gsel)}; H2 = {}
    with open(hpa) as f:
        r = csv.reader(f, delimiter="\t"); hdr = next(r)
        gi, li, vi = hdr.index("Gene name"), hdr.index("Cell line"), hdr.index("nTPM")
        for row_ in r:
            nm = norm(row_[li])
            if nm in want and row_[gi] in gp:
                H2.setdefault(want[nm], {})[row_[gi]] = np.log1p(float(row_[vi]))
    Lnames = [lines[c] for c in okl]; Mh = [i for i, c in enumerate(Lnames) if c in H2]
    hb = np.full((len(Mh), len(gsel)), np.nan)
    for kk, i in enumerate(Mh):
        for g, v in H2[Lnames[i]].items():
            hb[kk, gp[g]] = v
    okg = ~np.isnan(hb).any(0); hb = hb[:, okg] - hb[:, okg].mean(0)
    Fh = np.linalg.svd(hb, full_matrices=False)[0][:, :5]; Fh = Fh / Fh.std(0)
    def strength_stats(A, B, U):
        rep, r2 = [], []
        for j in range(A.shape[1]):
            u = U[:, j].mean(0); nu = np.linalg.norm(u)
            if nu < 1e-9:
                continue
            u /= nu; sa, sb = A[:, j] @ u, B[:, j] @ u
            rep.append(spearmanr(sa, sb)[0]); y = ((sa + sb) / 2)[Mh]; r2.append(kr_loo(Fh, y[:, None]))
        return dict(rep_median=float(np.median(rep)), r2_median=float(np.median(r2)), n=len(rep))
    res["strength_old_within_well"] = strength_stats(D[("plate6", 0)], D[("plate6", 1)], D[("plate14", 1)])
    res["strength_new_across_plates"] = strength_stats(D[("plate6", 0)], D[("plate14", 0)], D[("plate14", 1)])
    print("strength old", res["strength_old_within_well"], "new", res["strength_new_across_plates"], flush=True)
    json.dump(res, open(out + "/metrics_tahoe_plate_recompute.json", "w"), indent=1)
    print("saved", flush=True)


if __name__ == "__main__":
    main()
