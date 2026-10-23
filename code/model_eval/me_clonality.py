"""Model evaluation step 1: per-donor clonality from the TRBC marker (A100; cells.npz of cd4_sc_extract_cells.py).
Per cell f = TRBC1 / (TRBC1 + TRBC2) (>= MINTR counts); a clone carries one TCR beta rearrangement, so cells of one
clone share f, and a group drawn from few clones has a mean f away from the population mean f0.
  var_z     TRBC overdispersion at target level, as clone_tests_cd4.py: var over targets (>= 2 MINF valid cells) of
            z = (mean f - f0) / (sd_f / sqrt(n)); 1 under random assignment. It grows with cells per target as well
            as with clonality (var_z - 1 is about n times the clone variance of a target mean).
  vc        size-free clone variance of a guide mean, in units of the per-cell variance sd_f^2: mean over guides of
            ((mean f - f0)^2 - sd_f^2 / n) / sd_f^2 (guides with >= MINF valid cells). Subsampling cells keeps a
            guide's clone mix, so vc also holds for the pilot's capped subsets; vc_test is over the pilot test guides.
Usage: me_clonality.py <out_json> <donor>:<cells_rest>:<cells_stim>:<pilot_design.json> ..."""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json
import numpy as np

MINTR, MINF = 2, 15   # as clone_tests_cd4.py


def one(cd, test):
    meta = json.load(open(cd + "/meta.json")); groups = meta["groups"]; G = len(groups)
    c = np.load(cd + "/cells.npz"); grp, gp = c["grp"], c["gp"]; t1, t2 = c["trbc1"], c["trbc2"]
    tot = t1 + t2; ok = (grp >= 0) & (c["umi"] > 0) & (tot >= MINTR)
    f = np.where(tot > 0, t1 / np.maximum(tot, 1), np.nan)
    f0, sdf = float(np.nanmean(f[ok])), float(np.nanstd(f[ok]))
    s = np.bincount(grp[ok], weights=f[ok], minlength=G); n = np.bincount(grp[ok], minlength=G)
    tset = np.array([i for i in range(G - 1) if n[i] >= 2 * MINF])
    z = (s[tset] / n[tset] - f0) / (sdf / np.sqrt(n[tset]))
    key = grp[ok] * 2 + gp[ok]   # guide = (target, parity)
    sg = np.bincount(key, weights=f[ok], minlength=2 * G); ng = np.bincount(key, minlength=2 * G)
    gok = np.where((ng >= MINF) & (np.arange(2 * G) // 2 < G - 1))[0]
    ex = ((sg[gok] / ng[gok] - f0) ** 2 - sdf ** 2 / ng[gok]) / sdf ** 2
    gi = {g: i for i, g in enumerate(groups)}
    tk = {gi[t] * 2 + p for t in test if t in gi for p in (0, 1)}
    mt = np.array([k in tk for k in gok])
    vc, vct = float(ex.mean()), float(ex[mt].mean())
    return dict(var_z=float(np.var(z)), n_targets=int(len(tset)), median_cells_per_target=float(np.median(n[tset])),
                vc=vc, vc_se=float(ex.std(ddof=1) / np.sqrt(len(ex))), n_guides=int(len(ex)),
                vc_test=vct, vc_test_se=float(ex[mt].std(ddof=1) / np.sqrt(mt.sum())), n_test_guides=int(mt.sum()),
                f0=f0, sd_f=sdf)


def main():
    outp = sys.argv[1]; res = {"rules": __doc__}
    for arg in sys.argv[2:]:
        d, cr, cs, dj = arg.split(":")
        des = json.load(open(dj)); test = des["inert"] + des["strong"]
        res[d] = {"rest": one(cr, test), "stim48": one(cs, test)}
        print(d, json.dumps(res[d]), flush=True)
    json.dump(res, open(outp, "w"), indent=1)


if __name__ == "__main__":
    main()
