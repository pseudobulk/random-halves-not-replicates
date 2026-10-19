"""Non-targeting guides as null perturbations: random halves vs different guides (clone preprint v9, findings ca).
Inputs per donor-state: <ntc_dir> from cd4_ntc_extract.py and <cells_dir> from cd4_sc_extract_cells.py (TRBC counts).
  effect of guide g in half h = its mean log1p CP10K minus the mean over all kept guides in half h (null baseline)
  genes: mean over all non-targeting cells >= 0.1, TCR genes (TRAC, TRBC, TRAV, TRBV, TRGC, TRDC) removed
  guides: >= 10 cells in each r1 half
  random halves: guide g, r1 = 0 vs r1 = 1; different guides: random pairs (i, j), half r1 = 0 of i vs r1 = 1 of j
  (same depth), 50 random pairings; cross-half signal = sum over genes of the product of the two effects
  TRBC: per guide z of the mean TRBC1 fraction against random cells (cells with >= 2 TRBC counts, guides with >= 30),
  var(z) and after permuting guide labels among non-targeting cells; deviation correlated between r1 halves (>= 15)
Usage: cd4_ntc_null.py <out.json> <label>:<ntc_dir>:<cells_dir> [...]"""
import sys, json
import numpy as np

TCR = ("TRAC", "TRBC", "TRAV", "TRBV", "TRGC", "TRDC")


def share(v):
    return float((v > 0).mean())


def one(ntc, cells, seed=0):
    S = np.load(ntc + "/sum2.npy"); n = np.load(ntc + "/n2.npy"); m = json.load(open(ntc + "/meta.json"))
    genes = m["genes"]
    tot = S.sum((0, 1)) / n.sum()
    gk = (tot >= 0.1) & np.array([not g.startswith(TCR) for g in genes])
    ok = (n >= 10).all(1)
    M = S[ok][:, :, gk] / n[ok][:, :, None]
    X = M - M.mean(0, keepdims=True)                       # null baseline per half
    rnd = (X[:, 0] * X[:, 1]).sum(1)
    rng = np.random.default_rng(seed); K = X.shape[0]
    fp, mp = [], []
    for _ in range(50):
        p = rng.permutation(K); a, b = p[: K // 2], p[K // 2: 2 * (K // 2)]
        c = (X[a, 0] * X[b, 1]).sum(1); fp.append(share(c)); mp.append(float(c.mean()))
    boot = [share(rnd[rng.integers(0, K, K)]) for _ in range(1000)]
    res = dict(n_guides=int(K), n_genes=int(gk.sum()), cells_per_half_median=float(np.median(n[ok])),
               random_frac_pos=share(rnd), random_frac_pos_ci95=[float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
               random_mean_cross=float(rnd.mean()), pairs_frac_pos_mean=float(np.mean(fp)), pairs_frac_pos_range=[float(min(fp)), float(max(fp))],
               pairs_mean_cross=float(np.mean(mp)))
    # TRBC clone marker on the same guides' cells
    c = np.load(cells + "/cells.npz"); cm = json.load(open(cells + "/meta.json"))
    NG = len(cm["groups"]); names = cm["guide_ids"]
    sel = (c["grp"] == NG - 1) & (c["umi"] > 0)
    t1, t2 = c["trbc1"][sel], c["trbc2"][sel]; gid = c["gid"][sel]; r1 = c["r1"][sel]
    v = (t1 + t2) >= 2; f = t1[v] / (t1[v] + t2[v]); gid, r1 = gid[v], r1[v]
    f0, sd = f.mean(), f.std()
    def zvar(lab):
        u, inv, cnt = np.unique(lab, return_inverse=True, return_counts=True)
        mu = np.bincount(inv, f) / cnt; keep = cnt >= 30
        return float(np.var((mu[keep] - f0) / (sd / np.sqrt(cnt[keep])))), int(keep.sum())
    vz, nz = zvar(gid)
    vp = float(np.mean([zvar(rng.permutation(gid))[0] for _ in range(5)]))
    dev = []
    for h in (0, 1):
        s = r1 == h; u, inv, cnt = np.unique(gid[s], return_inverse=True, return_counts=True)
        dev.append(dict(zip(u, zip(np.bincount(inv, f[s]) / cnt - f0, cnt))))
    common = [g for g in dev[0] if g in dev[1] and dev[0][g][1] >= 15 and dev[1][g][1] >= 15]
    a = np.array([dev[0][g][0] for g in common]); b = np.array([dev[1][g][0] for g in common])
    res.update(trbc_var_z=vz, trbc_var_z_permuted=vp, trbc_n_guides=nz, trbc_random_half_pearson=float(np.corrcoef(a, b)[0, 1]),
               trbc_n_guides_halves=len(common))
    return res


def main():
    out = {"rules": __doc__}
    for arg in sys.argv[2:]:
        lab, ntc, cells = arg.split(":")
        out[lab] = one(ntc, cells); print(lab, json.dumps(out[lab]), flush=True)
    json.dump(out, open(sys.argv[1], "w"), indent=1)


if __name__ == "__main__":
    main()
