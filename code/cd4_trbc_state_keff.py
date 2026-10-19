"""TRBC clone marker: (cd) do resting and stimulated cells of a perturbation come from the same clones, and (ce) an
upper bound on the effective number of clones per perturbation (clone preprint, review v11 M2a/M2c; A100).
Inputs: per-cell extractions of cd4_sc_extract_cells.py (cells.npz + meta.json) for the two states of one donor.
  TRBC1 fraction f per cell (cells with >= 2 TRBC1 + TRBC2 counts); per perturbation (all guides pooled):
  deviation d = mean f - f0 and z = d / (sd_f / sqrt(n)), with f0, sd_f over all cells of the state.
  (cd) Pearson r of d between resting and stimulated over perturbations with >= 30 such cells in both states
       (matched by perturbation name), and the same between random halves within each state for reference.
  (ce) var(z) - 1 = n ICC / K_eff - ICC for K_eff equal-weight clones and intra-clone correlation ICC <= 1, so
       K_eff <= n / z^2 per perturbation (ICC = 1); reported as the median over perturbations with z^2 > 1, the bound
       from the pooled variance (mean n / var(z)), and the per-perturbation distribution.
Usage: cd4_trbc_state_keff.py <out.json> <label>:<rest_cells_dir>:<stim_cells_dir> [...]"""
import sys, json
import numpy as np


def per_pert(cdir, min_n=30):
    c = np.load(cdir + "/cells.npz"); m = json.load(open(cdir + "/meta.json"))
    names = m["groups"]; NG = len(names)
    ok = (c["grp"] >= 0) & (c["grp"] < NG - 1) & (c["umi"] > 0)
    t1, t2 = c["trbc1"][ok], c["trbc2"][ok]; g = c["grp"][ok]; r1 = c["r1"][ok]
    v = (t1 + t2) >= 2; f = t1[v] / (t1[v] + t2[v]); g, r1 = g[v], r1[v]
    f0, sd = f.mean(), f.std()
    out = {}
    cnt = np.bincount(g, minlength=NG); sm = np.bincount(g, f, minlength=NG)
    halves = [(np.bincount(g[r1 == h], minlength=NG), np.bincount(g[r1 == h], f[r1 == h], minlength=NG)) for h in (0, 1)]
    for k in np.where(cnt >= min_n)[0]:
        n = cnt[k]; d = sm[k] / n - f0
        dh = [hs[k] / hc[k] - f0 if hc[k] >= 15 else np.nan for hc, hs in halves]
        out[names[k]] = (int(n), float(d), float(d / (sd / np.sqrt(n))), dh[0], dh[1])
    return out


def main():
    res = {"rules": __doc__}
    for arg in sys.argv[2:]:
        lab, rest, stim = arg.split(":")
        R, S = per_pert(rest), per_pert(stim)
        common = sorted(set(R) & set(S))
        dr = np.array([R[k][1] for k in common]); ds = np.array([S[k][1] for k in common])
        r = {"n_common": len(common), "pearson_rest_vs_stim": float(np.corrcoef(dr, ds)[0, 1])}
        for nm, P in (("rest", R), ("stim", S)):
            a = np.array([P[k][3] for k in P]); b = np.array([P[k][4] for k in P]); okh = ~np.isnan(a) & ~np.isnan(b)
            r["pearson_random_halves_" + nm] = float(np.corrcoef(a[okh], b[okh])[0, 1])
            n = np.array([P[k][0] for k in P], float); z = np.array([P[k][2] for k in P])
            z2 = z ** 2; sel = z2 > 1
            kb = n[sel] / z2[sel]
            r["keff_" + nm] = dict(n_perts=int(len(n)), n_cells_median=float(np.median(n)), var_z=float(np.var(z)),
                                   frac_z2_gt1=float(sel.mean()), keff_upper_median=float(np.median(kb)),
                                   keff_upper_q25_q75=[float(np.percentile(kb, 25)), float(np.percentile(kb, 75))],
                                   keff_upper_pooled=float(n.mean() / np.var(z)))
        res[lab] = r; print(lab, json.dumps(r), flush=True)
    json.dump(res, open(sys.argv[1], "w"), indent=1)


if __name__ == "__main__":
    main()
