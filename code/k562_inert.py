"""Does the clone effect of cell halves (Z1 ah) also sit in Replogle's data, the base of the gamma ledger? (A100.)
K562 genome-wide Perturb-seq (Replogle 2022, raw single cells part 01). Each gene has one dual-guide construct, so a guide
split is not possible; instead compare two split designs on inert targets (target gene absent from the expressed-gene
list of the file, or control mean log1p CP10K < 0.01) vs expressed targets (>= 0.5):
  random: the target's cells split at random (as line_parts / the F26 ledger)
  lane:   the target's cells split by 10x lane (gem group parity); controls split the same way
If inert targets show reproducible responses under the random split, the clone effect is present here too; if the lane
split removes it, lane-split denoising suffices for this design, otherwise neither does.
Response genes: top N_GENES by control mean. CP10K uses the per-cell UMI_count from obs. Template (mean over targets
per half) removed before the per-target cross-half signal sum(x_A * x_B).
10-18: added the inert-referenced check (see the comment in main); the original numbers are unchanged.
Usage: k562_inert.py <h5ad> <out_dir>. Exploratory."""
import os, sys, json
import numpy as np
import h5py
from explore_V_day6_day8 import obs_info, read_rows, CTRL_NAMES

N_GENES, MIN_HALF = 2000, 15


def main():
    path, out = sys.argv[1:3]; os.makedirs(out, exist_ok=True)
    codes, cats, gem, genes = obs_info(path)
    with h5py.File(path) as f:
        o = f["obs"]
        umi = o["UMI_count"][:] if "UMI_count" in o else None
    ctrl = [c for c in CTRL_NAMES if c in cats][0]; ic = cats.index(ctrl)
    gem = np.asarray(gem); gnum = np.array([int(str(g).split("_")[-1]) if str(g).split("_")[-1].isdigit() else hash(str(g)) for g in gem])
    print("cells", len(codes), "targets", len(cats), "genes", len(genes), "ctrl", ctrl, "umi in obs", umi is not None, flush=True)
    crow = np.where(codes == ic)[0]
    Xc = read_rows(path, crow, np.arange(len(genes)))
    if umi is None:
        raise SystemExit("UMI_count missing")
    Yc = np.log1p(Xc / umi[crow][:, None] * 1e4); cm = Yc.mean(0)
    top = np.argsort(-cm)[:N_GENES]; gpos = {g: j for j, g in enumerate(genes)}
    counts = np.bincount(codes, minlength=len(cats))
    targets = [i for i, c in enumerate(cats) if i != ic and counts[i] >= 2 * MIN_HALF]
    texp = np.array([cm[gpos[cats[t]]] if cats[t] in gpos else 0.0 for t in targets])
    inert, expd = texp < 0.01, texp >= 0.5
    print("targets", len(targets), "inert", int(inert.sum()), "expressed", int(expd.sum()), flush=True)
    rows = np.where(np.isin(codes, targets))[0]
    X = read_rows(path, rows, top)
    Y = np.log1p(X / umi[rows][:, None] * 1e4)
    Yc_top = Yc[:, top]
    rng = np.random.default_rng(9)
    res = {"rules": __doc__, "n_targets": len(targets), "n_inert": int(inert.sum()), "n_expressed": int(expd.sum())}
    tpos = {t: k for k, t in enumerate(targets)}; rc = codes[rows]
    for design in ("random", "lane"):
        if design == "random":
            ch = rng.integers(0, 2, len(crow)); hh = rng.integers(0, 2, len(rows))
        else:
            ch = gnum[crow] % 2; hh = gnum[rows] % 2
        cmu = [Yc_top[ch == h].mean(0) for h in range(2)]
        XA = np.full((len(targets), N_GENES), np.nan); XB = np.full_like(XA, np.nan)
        for t in targets:
            m = rc == t
            a, b = m & (hh == 0), m & (hh == 1)
            if a.sum() >= MIN_HALF and b.sum() >= MIN_HALF:
                XA[tpos[t]] = Y[a].mean(0) - cmu[0]; XB[tpos[t]] = Y[b].mean(0) - cmu[1]
        ok = ~np.isnan(XA).any(1)
        # 10-18 check: with the all-target template removed, a null target equals minus the template in both halves,
        # which is "reproducible" without any clone. Also report the template's own cross energy and the inert-referenced
        # version (responses minus the mean of the inert targets = the null-perturbation baseline).
        XAr, XBr = XA.copy(), XB.copy()
        TA, TB = np.nanmean(XAr[ok], 0), np.nanmean(XBr[ok], 0)
        mi = inert & ok
        IA, IB = XAr[mi].mean(0), XBr[mi].mean(0)
        sig_i = np.where(ok, (np.nan_to_num(XAr - IA) * np.nan_to_num(XBr - IB)).sum(1), np.nan)
        ref = {"template_cross": float((TA * TB).sum()), "inert_mean_cross": float((IA * IB).sum())}
        for name, m in (("inert", inert & ok), ("expressed", expd & ok)):
            ref[name] = dict(signal_median=float(np.median(sig_i[m])), frac_pos=float((sig_i[m] > 0).mean()))
        XA, XB = XA - np.nanmean(XA[ok], 0), XB - np.nanmean(XB[ok], 0)
        sig = np.where(ok, (np.nan_to_num(XA) * np.nan_to_num(XB)).sum(1), np.nan)
        d = {"n_ok": int(ok.sum()), "inert_referenced": ref}
        for name, m in (("inert", inert & ok), ("expressed", expd & ok)):
            d[name] = dict(n=int(m.sum()), signal_median=float(np.median(sig[m])), frac_pos=float((sig[m] > 0).mean()))
        A = XA[ok] - XA[ok].mean(0); V = np.linalg.svd(A, full_matrices=False)[2][:3].T
        pa, pb = A @ V, (XB[ok] - XB[ok].mean(0)) @ V
        io = (inert & ok)[ok]
        d["inert_pc_crosshalf_r"] = [float(np.corrcoef(pa[io, j], pb[io, j])[0, 1]) for j in range(3)]
        res[design] = d
        print(design, json.dumps(d), flush=True)
    json.dump(res, open(out + "/metrics_k562_inert.json", "w"), indent=1)
    print("saved", flush=True)


if __name__ == "__main__":
    main()
