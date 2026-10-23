"""The same clone ruler for every generality dataset (as for the whole-brain and VCC control tests).
Expression: log1p CP10K; genes with mean >= 0.1 over null cells. Effect of a unit (guide or guide construct) in a half:
mean of its cells minus the mean of the null cells of the same sample (library / mouse / donor), weighted by the unit's
cells per sample; for a null unit its own cells are left out of the reference.
Null units (non-targeting / intergenic / safe-targeting guides), >= MINC cells:
  random     cells split by a random bit; cross-half signal sum_g x_A,g x_B,g averaged over NSPLIT splits
  diff       first half of one null unit against the second half of another (random pairings; noise level)
  replicate  where labelled (mice, donors): replicate groups split at random (NSPLIT draws), same signal
  scaled     mean random signal / summed per-cell variance of null cells around their sample mean (= vc_expr scale)
Targeting units: the same random-half signal (mean and median over units); ratio null / targeting. Where a target has
>= 2 units, the cross-unit signal (one unit against another of the same target, independent transductions as clones
carry one construct) over the random-half signal of the same targets = share reproduced by another guide.
Usage: imported by run_generality.py"""
import numpy as np


def lognorm(X):
    tot = np.asarray(X.sum(1)).ravel().astype(np.float64)
    X = X.multiply(1e4 / np.maximum(tot, 1)[:, None]).tocsr().astype(np.float32); X.data = np.log1p(X.data)
    return X


def ruler(X, sample, unit, target, is_null, rep=None, MINC=30, NSPLIT=20, seed=0, ref_mode="null"):
    """X: csr cells x genes (log1p CP10K); sample, unit, target, rep: per-cell arrays (unit None = unassigned);
    is_null: per-cell bool (cell belongs to a null unit).
    ref_mode 'null': reference = null cells of the sample (leave-one-out for null units). With few null units this
    reference carries the clone states of the other null units, which biases the random-half signal up and the
    different-guide signal down. 'all': reference = all assigned cells of the sample, leaving out the unit's own
    cells; each unit's clones then weigh ~1/(number of units) in it, and a shift by the average targeting effect enters
    the random-half and the different-guide signals alike, so their difference is the clone part.
    The clone part of null units is reported as (random - diff) / V."""
    rng = np.random.default_rng(seed)
    unit = np.asarray(unit, dtype=object); ok = np.array([u is not None for u in unit])
    null = ok & is_null
    gm = np.asarray(X[null].mean(0)).ravel() >= 0.1
    Xg = X[:, gm].tocsr()
    us, si = np.unique(sample, return_inverse=True)
    refcells = null if ref_mode == "null" else ok
    Sn = np.zeros((len(us), gm.sum())); Nn = np.zeros(len(us))
    for s in range(len(us)):
        r = np.where(refcells & (si == s))[0]
        if len(r):
            Sn[s] = np.asarray(Xg[r].sum(0)).ravel(); Nn[s] = len(r)
    mu = Sn / np.maximum(Nn, 1)[:, None]
    Sv = np.zeros((len(us), gm.sum())); Nv = np.zeros(len(us))            # null-cell means for the variance
    for s in range(len(us)):
        r = np.where(null & (si == s))[0]
        if len(r):
            Sv[s] = np.asarray(Xg[r].sum(0)).ravel(); Nv[s] = len(r)
    muv = Sv / np.maximum(Nv, 1)[:, None]
    # per-cell variance of null cells around their sample mean
    V = 0.0; nv = 0
    for s in range(len(us)):
        r = np.where(null & (si == s))[0]
        if len(r) > 1:
            V += float(((Xg[r].toarray() - muv[s]) ** 2).sum()); nv += len(r)
    V /= max(nv, 1)

    def effect(idx, own_null):
        E = Xg[idx].toarray().astype(np.float64)
        if own_null or ref_mode == "all":
            own = np.zeros_like(Sn); np.add.at(own, si[idx], E); on = np.bincount(si[idx], minlength=len(us))
            ref = (Sn - own) / np.maximum(Nn - on, 1)[:, None]
        else:
            ref = mu
        return E - ref[si[idx]]

    units = {}
    for u in np.unique(unit[ok].astype(str)):
        idx = np.where(ok & (unit.astype(str) == u))[0]
        if len(idx) >= MINC:
            units[u] = idx
    res_null, res_tg, halves = {}, {}, {}
    for u, idx in units.items():
        isn = bool(is_null[idx[0]]); E = effect(idx, isn)
        cr = []
        for k in range(NSPLIT):
            bit = rng.integers(0, 2, len(idx)).astype(bool)
            if bit.sum() < 5 or (~bit).sum() < 5:
                continue
            a, b = E[bit].mean(0), E[~bit].mean(0); cr.append(float(a @ b))
            if k == 0:
                halves[u] = (a, b)
        r = dict(n=int(len(idx)), random=float(np.mean(cr)), target=str(target[idx[0]]), mean=E.mean(0))
        if rep is not None:
            ru = np.unique(rep[idx]); rr = []
            if len(ru) >= 2:
                for _ in range(NSPLIT):
                    g = np.isin(rep[idx], rng.permutation(ru)[: len(ru) // 2])
                    if g.sum() >= 5 and (~g).sum() >= 5:
                        rr.append(float(E[g].mean(0) @ E[~g].mean(0)))
            r["replicate"] = float(np.mean(rr)) if rr else np.nan
        (res_null if isn else res_tg)[u] = r
    out = {"ref_mode": ref_mode, "n_genes": int(gm.sum()), "V_cell": V, "n_null_units": len(res_null), "n_target_units": len(res_tg)}
    if res_null:
        nn = list(res_null); rnd = np.array([res_null[u]["random"] for u in nn])
        diff = []
        for _ in range(50):
            p = rng.permutation(len(nn))
            for i in range(len(p) // 2):
                a, b = nn[p[2 * i]], nn[p[2 * i + 1]]
                if a in halves and b in halves:
                    diff.append(float(halves[a][0] @ halves[b][1]))
        diff = np.array(diff)
        out.update(null_share_pos_random=float((rnd > 0).mean()), null_share_pos_diff=float((diff > 0).mean()),
                   null_mean_random=float(rnd.mean()), null_se_random=float(rnd.std() / np.sqrt(len(rnd))),
                   null_mean_diff=float(diff.mean()), null_scaled=float(rnd.mean() / V),
                   null_clone_scaled=float((rnd.mean() - diff.mean()) / V),
                   null_cells_median=float(np.median([res_null[u]["n"] for u in nn])))
        if rep is not None:
            rp = np.array([res_null[u]["replicate"] for u in nn]); m = np.isfinite(rp)
            if m.any():
                out.update(null_mean_replicate=float(rp[m].mean()), null_share_pos_replicate=float((rp[m] > 0).mean()),
                           null_n_replicate=int(m.sum()))
    if res_tg:
        tr = np.array([v["random"] for v in res_tg.values()])
        out.update(target_mean_random=float(tr.mean()), target_median_random=float(np.median(tr)),
                   target_cells_median=float(np.median([v["n"] for v in res_tg.values()])))
        if res_null:
            out["ratio_null_to_target_mean"] = out["null_mean_random"] / out["target_mean_random"]
            out["ratio_null_to_target_median"] = out["null_mean_random"] / out["target_median_random"]
        by_t = {}
        for u, v in res_tg.items():
            by_t.setdefault(v["target"], []).append(u)
        cross, rnd2 = [], []
        for t, us_ in by_t.items():
            if len(us_) >= 2:
                a, b = res_tg[us_[0]], res_tg[us_[1]]
                cross.append(float(a["mean"] @ b["mean"])); rnd2.append((a["random"] + b["random"]) / 2)
        if cross:
            out.update(n_targets_two_units=len(cross), target_cross_unit=float(np.mean(cross)),
                       target_random_same=float(np.mean(rnd2)), share_reproduced_by_other_guide=float(np.mean(cross) / np.mean(rnd2)))
        if rep is not None:
            rp = np.array([v.get("replicate", np.nan) for v in res_tg.values()]); m = np.isfinite(rp)
            if m.any():
                out.update(target_mean_replicate=float(rp[m].mean()), target_replicate_over_random=float(rp[m].mean() / tr[m].mean()))
    return out
