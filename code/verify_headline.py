"""Independent recomputation of the headline numbers (preprint correctness check 1; A100). Written from scratch, no
imports from the earlier analysis modules, and deliberately different implementation choices:
  - CD4 data from the per-cell extraction cd4_sc_extract_cells.py (its random halves r1 are unrelated to the halves
    of the original extraction); guide halves = gp.
  - Readout: state = top-k right singular vectors of half-A beta on the training perturbations; ridge solved as an
    augmented least-squares problem (not the normal equations); lambda (absolute, scaled by the mean squared state)
    chosen by 5-fold CV inside the training set; 5-fold outer CV over perturbations (earlier: 2-fold / 2-fold).
  - Tahoe rebuilt directly from tahoe_agg with a different line filter (line kept if usable in >= 90% of the shared
    conditions; earlier 85%).
Reported: CD4 - inert frac_pos per state (inert baseline), gamma share, readout share k = 20 and shuffled, for random
vs guide halves; Tahoe - fixed-target readout (target = two-way gamma of plate P half 1, scored against the plate's
within-well gamma cross product) with beta from the same well (half 0) vs the other plate (pooled).
Usage: verify_headline.py cd4 <cells_Rest> <cells_Stim> <out_dir> | tahoe <tahoe_agg> <out_dir>. Exploratory."""
import os, sys, json, ast
import numpy as np
import config

K, FOLDS, LAMS = 20, 5, (1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0)


def ridge_lstsq(Z, Y, lam):
    """argmin |Y - Z W|^2 + lam * s |W|^2 with s = mean squared entry of Z, via augmented least squares"""
    s = float((Z ** 2).mean()) * Z.shape[0]
    Za = np.vstack([Z, np.sqrt(lam * s) * np.eye(Z.shape[1])])
    Ya = np.vstack([Y, np.zeros((Z.shape[1], Y.shape[1]))])
    return np.linalg.lstsq(Za, Ya, rcond=None)[0]


def folds(n, k, rng):
    p = rng.permutation(n); return [p[i::k] for i in range(k)]


def fit_predict(F_tr, Y_tr, F_te, k, rng):
    V = np.linalg.svd(F_tr - 0, full_matrices=False)[2][:k].T
    Ztr, Zte = F_tr @ V, F_te @ V
    best, bl = -np.inf, 1.0
    inner = folds(len(Ztr), FOLDS, rng)
    for lam in LAMS:
        e = 0.0
        for f in inner:
            tr = np.setdiff1d(np.arange(len(Ztr)), f)
            P = Ztr[f] @ ridge_lstsq(Ztr[tr], Y_tr[tr], lam)
            e += (Y_tr[f] ** 2).sum() - ((Y_tr[f] - P) ** 2).sum()
        if e > best:
            best, bl = e, lam
    return Zte @ ridge_lstsq(Ztr, Y_tr, bl)


def cv_readout(F, Y, k, seed, shuffle=False):
    rng = np.random.default_rng(seed); n = len(Y)
    Fs = F[rng.permutation(n)] if shuffle else F
    P = np.zeros_like(Y)
    for f in folds(n, FOLDS, rng):
        tr = np.setdiff1d(np.arange(n), f)
        P[f] = fit_predict(Fs[tr], Y[tr], Fs[f], k, rng)
    return P


def cd4(dr, ds, out):
    def load(d):
        m = json.load(open(d + "/meta.json")); return m["genes"], m["groups"], np.load(d + "/sum8.npy", mmap_mode="r"), np.load(d + "/n8.npy")
    st = {"rest": load(dr), "stim": load(ds)}
    genes = sorted(set(st["rest"][0]) & set(st["stim"][0]))
    designs = {"random": ([0, 1, 4, 5], [2, 3, 6, 7]), "guide": ([0, 1, 2, 3], [4, 5, 6, 7])}
    tsets, cms = [], {}
    for k, (gl, groups, S, n8) in st.items():
        gi = {g: j for j, g in enumerate(gl)}; col = np.array([gi[g] for g in genes]); ci = len(groups) - 1
        cms[k] = np.asarray(S[ci].sum(0))[col] / n8[ci].sum()
        ok = set()
        for i, t in enumerate(groups[:-1]):
            if all(n8[i, a].sum() >= 10 and n8[i, b].sum() >= 10 for a, b in designs.values()):
                ok.add(t)
        tsets.append(ok)
    gpos = {g: j for j, g in enumerate(genes)}
    targets = sorted(t for t in tsets[0] & tsets[1] if t in gpos)
    keep = np.where(np.minimum(cms["rest"], cms["stim"]) >= 0.1)[0]; kpos = {j: q for q, j in enumerate(keep)}
    texp = np.array([max(cms["rest"][gpos[t]], cms["stim"][gpos[t]]) for t in targets]); inert = texp < 0.01   # both states
    res = {"n_targets": len(targets), "n_inert": int(inert.sum()), "n_genes": int(len(keep))}
    Sg_by = {}
    for k, (gl, groups, S, n8) in st.items():
        gi = {g: j for j, g in enumerate(gl)}; col = np.array([gi[g] for g in genes])[keep]; ci = len(groups) - 1
        ti = {t: i for i, t in enumerate(groups)}; idx = np.array([ti[t] for t in targets])
        Sg_by[k] = (np.asarray(S[idx])[:, :, col].astype(np.float64), n8[idx].astype(float),
                    np.asarray(S[ci])[:, col].astype(np.float64), n8[ci].astype(float))
    for dn, (a, b) in designs.items():
        H = {}
        for k in st:
            Sg, nt, Sc, nc = Sg_by[k]
            sides = []
            for codes in (a, b):
                X = Sg[:, codes].sum(1) / nt[:, codes].sum(1)[:, None] - Sc[codes].sum(0) / nc[codes].sum()
                for r_, t in enumerate(targets):
                    q = kpos.get(gpos[t])
                    if q is not None:
                        X[r_, q] = 0
                sides.append(X)
            H[k] = sides
        d = {}
        for k in ("rest", "stim"):
            A, B = H[k]; Ai, Bi = A - A[inert].mean(0), B - B[inert].mean(0)
            d["inert_frac_pos_" + k] = float(((Ai * Bi).sum(1)[inert] > 0).mean())
        c = lambda M: M - M.mean(0)
        XA, XB, YA, YB = c(H["rest"][0]), c(H["rest"][1]), c(H["stim"][0]), c(H["stim"][1])
        bA, gA, gB = (XA + YA) / 2, (XA - YA) / 2, (XB - YB) / 2
        den = float((gA * gB).sum()); d["gamma_share"] = den / (den + float((bA * (XB + YB) / 2).sum()))
        for nm, sh in (("readout_k20", False), ("readout_k20_shuffled", True)):
            P = cv_readout(bA, gB, K, 7, sh); d[nm] = float(((gB ** 2).sum() - ((gB - P) ** 2).sum()) / den)
        res[dn] = d
        print(dn, json.dumps(d), flush=True)
    return res


def tahoe(agg, out):
    import pyarrow.parquet as pq
    idx = json.load(open(agg + "/index.json")); lines, samples = idx["lines"], idx["samples"]
    C, Ps = len(lines), len(samples)
    S = np.load(agg + "/sums64.npy", mmap_mode="r"); n = np.load(agg + "/counts.npy")[:, 0]
    sm = pq.read_table(config.tahoe_meta("sample_metadata.parquet")).to_pydict()
    meta = {s: (d.strip(), str(ast.literal_eval(dc)[0][1]), p) for s, d, dc, p in zip(sm["sample"], sm["drug"], sm["drugname_drugconc"], sm["plate"])}
    row = lambda c, s, h: (c * Ps + s) * 2 + h
    byp = {p: {(meta[samples[i]][0], meta[samples[i]][1]): i for i in range(Ps) if meta[samples[i]][2] == p} for p in ("plate6", "plate14")}
    dm = {p: [i for i in range(Ps) if meta[samples[i]][2] == p and meta[samples[i]][0].startswith("DMSO")] for p in byp}
    conds = sorted(k for k in byp["plate6"] if k in byp["plate14"] and not k[0].startswith("DMSO"))
    good = np.array([[all(n[row(c, byp[p][k], h)] >= 10 for p in byp for h in (0, 1)) for k in conds] for c in range(C)])
    lk = [c for c in range(C) if good[c].mean() >= 0.90]
    ck = [j for j in range(len(conds)) if good[lk, j].all()]
    rr = [row(c, s, h) for c in lk for s in dm["plate6"] for h in (0, 1)]
    gm = np.asarray(S[rr].sum(0)) / n[rr].sum(); gk = np.where(gm >= 0.1)[0]
    def eff(p, hs, c, k):
        r = [row(c, byp[p][conds[k]], h) for h in hs]; rd = [row(c, s, h) for s in dm[p] for h in hs]
        return np.asarray(S[r].sum(0))[gk] / n[r].sum() - np.asarray(S[rd].sum(0))[gk] / n[rd].sum()
    D = {(p, hn): np.array([[eff(p, hs, c, k) for k in ck] for c in lk]) for p in byp for hn, hs in (("h0", (0,)), ("h1", (1,)), ("pool", (0, 1)))}
    print("lines", len(lk), "conditions", len(ck), "genes", len(gk), flush=True)
    def twoway(M):
        m = M.mean((0, 1)); return M - M.mean(1, keepdims=True) - M.mean(0, keepdims=True) + m, M.mean(0) - m
    res = {"n_lines": len(lk), "n_conditions": len(ck), "n_genes": int(len(gk))}
    for P, O in (("plate14", "plate6"), ("plate6", "plate14")):
        g1, _ = twoway(D[(P, "h1")]); g0, b0 = twoway(D[(P, "h0")]); _, bo = twoway(D[(O, "pool")])
        den = float((g0 * g1).sum()); d = {}
        for nm, F in (("same_well", b0), ("other_plate", bo)):
            e = 0.0
            for c in range(len(lk)):
                Pc = cv_readout(F, g1[c], 10, 11 + c)
                e += float((g1[c] ** 2).sum() - ((g1[c] - Pc) ** 2).sum())
            d[nm] = e / den
        res["target_" + P] = d
        print(P, json.dumps(d), flush=True)
    return res


def main():
    mode = sys.argv[1]
    if mode == "cd4":
        out = sys.argv[4]; os.makedirs(out, exist_ok=True); res = cd4(*sys.argv[2:5])
    else:
        out = sys.argv[3]; os.makedirs(out, exist_ok=True); res = tahoe(*sys.argv[2:4])
    res["rules"] = __doc__
    json.dump(res, open(out + "/metrics_verify_headline.json", "w"), indent=1)
    print("saved", flush=True)


if __name__ == "__main__":
    main()
