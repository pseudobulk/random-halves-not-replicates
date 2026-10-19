"""Clone or guide? The reproducible "response" of K562 targets that are not expressed in K562 (A100; arrays from
replogle_indep_readout.py: K562 day 6 essential screen K6, day 8 genome-wide screen K8, same dual-guide elements for
these genes, independent transductions).
Within one screen, an inert target's response reproducible across random cell halves can come from shared clones or
from the guide itself (off-target binding, neighbouring promoters). Across the two screens only guide-level structure
can reproduce (different transductions, so different clones; also different time point and batch). For inert targets
(K562 control mean log1p CP10K < 0.01 in both screens, or absent from a screen's gene list), responses referenced to
the mean of the inert targets (the null-perturbation baseline, not the all-target template; Z1 an):
  within_K6, within_K8   sum over inert targets of half A . half B (same screen)
  cross                  (K6 half A . K8 half B + K6 half B . K8 half A) / 2
  guide_fraction         cross / sqrt(within_K6 * within_K8)  ~ 0 -> clone-level; ~ 1 -> guide-level
Expressed targets (>= 0.5 in both) as reference (their cross-screen ratio is the K6-K8 signal correlation of real
responses, 0.83 for all targets in Z1 al). Bootstrap over inert targets (1000).
Usage: k562_inert_screens.py <responses.npz> <K6.h5ad> <K8.h5ad> <out_dir>. Exploratory."""
import os, sys, json
import numpy as np
import h5py
from explore_V_day6_day8 import obs_info, read_rows, CTRL_NAMES


def target_expr(path, perts):
    codes, cats, gem, genes = obs_info(path)
    with h5py.File(path) as f:
        umi = f["obs"]["UMI_count"][:]
    ctrl = [c for c in CTRL_NAMES if c in cats][0]
    crow = np.where(codes == list(cats).index(ctrl))[0]
    own = {g: j for j, g in enumerate(genes)}; tg = [t for t in perts if t in own]
    X = read_rows(path, crow, np.array([own[t] for t in tg]))
    e = dict(zip(tg, np.log1p(X / umi[crow][:, None] * 1e4).mean(0)))
    return np.array([e.get(t, 0.0) for t in perts])


def main():
    npz, p6, p8, out = sys.argv[1:5]; os.makedirs(out, exist_ok=True)
    z = np.load(npz); perts = list(z["perts"])
    e6, e8 = target_expr(p6, perts), target_expr(p8, perts)
    inert = (e6 < 0.01) & (e8 < 0.01); expd = (e6 >= 0.5) & (e8 >= 0.5)
    # responses.npz holds template-removed responses; adding back the template and subtracting the inert mean equals
    # referencing to the inert baseline directly (the template is a constant vector per half)
    H6, H8 = z["K6_H"].astype(np.float64), z["K8_H"].astype(np.float64)
    ref = lambda H: H - H[:, inert].mean(1, keepdims=True)
    R6, R8 = ref(H6), ref(H8)
    def parts(m):
        w6 = (R6[0][m] * R6[1][m]).sum(1); w8 = (R8[0][m] * R8[1][m]).sum(1)
        x = ((R6[0][m] * R8[1][m]).sum(1) + (R6[1][m] * R8[0][m]).sum(1)) / 2
        return w6, w8, x
    res = {"rules": __doc__, "n_perturbations": len(perts), "n_inert": int(inert.sum()), "n_expressed": int(expd.sum()),
           "inert_targets": [perts[i] for i in np.where(inert)[0]]}
    rng = np.random.default_rng(3)
    for name, m in (("inert", inert), ("expressed", expd)):
        w6, w8, x = parts(m); n = len(x)
        gf = lambda ix: x[ix].sum() / np.sqrt(w6[ix].sum() * w8[ix].sum())
        boot = [gf(rng.integers(0, n, n)) for _ in range(1000)]
        res[name] = dict(n=n, within_K6=float(w6.sum()), within_K8=float(w8.sum()), cross=float(x.sum()),
                         guide_fraction=float(gf(np.arange(n))), ci95=[float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
                         frac_pos_K6=float((w6 > 0).mean()), frac_pos_K8=float((w8 > 0).mean()), frac_pos_cross=float((x > 0).mean()),
                         energy_per_target_K6=float(w6.mean()), energy_per_target_K8=float(w8.mean()))
        print(name, json.dumps(res[name]), flush=True)
    json.dump(res, open(out + "/metrics_k562_inert_screens.json", "w"), indent=1)
    print("saved", flush=True)


if __name__ == "__main__":
    main()
