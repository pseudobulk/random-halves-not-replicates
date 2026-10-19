"""Does the readout pass non-biological (clone / guide) structure into its gamma prediction? A within-screen test for
pairs without independent units (Nadig Jurkat-HepG2; the RPE1 side of K562-RPE1). (5090 / A100)
Pair of contexts X, Y (two lines, or CD4 Rest / Stim); responses per random half (template removed); beta = (X + Y) / 2,
gamma = (X - Y) / 2. Readout as (r): state = top-k PCs of half-A beta on training perturbations, ridge (lambda by inner
2-fold CV) to half-B gamma; predictions F for held-out perturbations (2 folds over perturbations, seed 11).
Inert in X = target gene not expressed in X (absent from X's gene list, or control mean log1p CP10K < 0.01). For such
targets X's reproducible response (X_A . X_B > 0) cannot be the knockdown's biology; it is clone / guide structure.
Gamma contains +X / 2, so a readout that reads half-A clone deviations (they enter beta_A) passes them on to its
prediction. Measure, per context and target set S:
  pass_X(S) = 2 * sum_S F . X_B / sum_S X_A . X_B   (sign -1 for Y, since gamma contains -Y / 2)
the share of X's within-screen reproducible response energy that reappears in the gamma prediction. Expressed targets
(>= 0.5) give the same quantity for real responses. Shuffled control: features permuted across perturbations before
fitting and predicting. Bootstrap over the targets in S (500) for a 95% interval.
Reading: pass(inert) near 0 -> the readout does not read clone structure in this screen; pass(inert) comparable to
pass(expressed) -> it does. Calibration on CD4 donor-1 single cells: cell split (coupling present, Z1 am) vs guide
split (coupling absent).
Usage: readout_inert_line.py h5ad <lineX.h5ad> <lineY.h5ad> <labelX> <labelY> <out_dir>
       readout_inert_line.py cd4 <out_Rest> <out_Stim> <out_dir>. Exploratory."""
import os, sys, json, time
import numpy as np
import h5py

KS = (10, 20)
LAMS = (1e-3, 1e-2, 1e-1, 1, 10, 100)
MIN_EXPR, NBOOT = 0.1, 500


def ridge(Z, Y, lam):
    k = Z.shape[1]
    return np.linalg.solve(Z.T @ Z + lam * np.trace(Z.T @ Z) / k * np.eye(k), Z.T @ Y)


def pick_lam(Z, Y, rng):
    P = len(Z); perm = rng.permutation(P); f = (perm[: P // 2], perm[P // 2:])
    best, bl = -np.inf, 1.0
    for lam in LAMS:
        e = 0.0
        for i in (0, 1):
            F = Z[f[1 - i]] @ ridge(Z[f[i]], Y[f[i]], lam); e += ((Y[f[1 - i]] ** 2).sum() - ((Y[f[1 - i]] - F) ** 2).sum())
        if e > best:
            best, bl = e, lam
    return bl


def predict(bA, gB, k, seed=11):
    P = len(gB); rng = np.random.default_rng(seed); perm = rng.permutation(P); folds = (perm[: P // 2], perm[P // 2:])
    F = np.zeros_like(gB)
    for i in (0, 1):
        tr, te = folds[i], folds[1 - i]
        V = np.linalg.svd(bA[tr], full_matrices=False)[2][:k].T; Z = bA @ V
        F[te] = Z[te] @ ridge(Z[tr], gB[tr], pick_lam(Z[tr], gB[tr], rng))
    return F


def analyse(XAr, XBr, YAr, YBr, sets, labels):
    """XAr..YBr: perturbations x genes, responses vs control (no template removal). sets: name -> boolean masks.
    The readout is fitted on template-removed responses as in (r). The pass measure uses responses referenced to the
    mean of that context's inert targets (the null-perturbation baseline), not to the mean over all perturbations:
    with an essential library the all-perturbation template is large, and an inert target minus that template equals
    minus the template in both halves, which looks perfectly "reproducible" without any clone (v1 artefact)."""
    c0 = lambda X: X - X.mean(0)
    XA, XB, YA, YB = c0(XAr), c0(XBr), c0(YAr), c0(YBr)
    ref = {}
    for lab, HAr, HBr in ((labels[0], XAr, XBr), (labels[1], YAr, YBr)):
        m = sets[lab]["inert"]
        if m.sum() < 10:
            ref[lab] = (None, None); continue
        # inert baseline, then remove the template direction (mean over all perturbations) so that expressed and
        # inert targets are scored on the same footing: the template-removed readout cannot predict template energy
        # (10-18 audit fix)
        prj = []
        for H in (HAr, HBr):
            R0 = H - H[m].mean(0); u = H.mean(0) - H[m].mean(0); u = u / np.linalg.norm(u)
            prj.append(R0 - np.outer(R0 @ u, u))
        ref[lab] = tuple(prj)
    bA, gA, gB = (XA + YA) / 2, (XA - YA) / 2, (XB - YB) / 2
    out = {"gamma_share": float((gA * gB).sum() / ((gA * gB).sum() + (bA * (XB + YB) / 2).sum()))}
    rng = np.random.default_rng(5)
    for k in KS:
        F = predict(bA, gB, k); Fs = predict(bA[rng.permutation(len(bA))], gB, k)
        d = {"readout_share": float(((gB ** 2).sum() - ((gB - F) ** 2).sum()) / (gA * gB).sum()),
             "readout_share_shuffled": float(((gB ** 2).sum() - ((gB - Fs) ** 2).sum()) / (gA * gB).sum())}
        for lab, sign in ((labels[0], 1.0), (labels[1], -1.0)):
            HA, HB = ref[lab]
            if HA is None:
                continue
            c = sign * 2 * (F * HB).sum(1); cs = sign * 2 * (Fs * HB).sum(1); w = (HA * HB).sum(1)
            for sname, m in sets[lab].items():
                ix = np.where(m)[0]
                if len(ix) < 10:
                    d["%s_%s" % (lab, sname)] = {"n": int(len(ix))}; continue
                boot = []
                for _ in range(NBOOT):
                    b = rng.choice(ix, len(ix)); boot.append(c[b].sum() / w[b].sum())
                d["%s_%s" % (lab, sname)] = dict(n=int(len(ix)), pass_share=float(c[ix].sum() / w[ix].sum()),
                                                 ci95=[float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
                                                 pass_share_shuffled=float(cs[ix].sum() / w[ix].sum()),
                                                 reproducible_energy=float(w[ix].sum()), frac_pos=float((w[ix] > 0).mean()))
        out["k%d" % k] = d
    return out


def run_h5ad(px, py, lx, ly, out):
    from explore_V_day6_day8 import obs_info, read_rows, Line, CTRL_NAMES, PRO, P53, MIN_CELLS
    t0 = time.time(); info = {}
    for lab, p in ((lx, px), (ly, py)):
        codes, cats, gem, genes = obs_info(p)
        with h5py.File(p) as f:
            umi = f["obs"]["UMI_count"][:]
        info[lab] = dict(codes=codes, cats=list(cats), gem=gem, genes=list(genes), umi=umi, path=p)
    excl = set(P53) | set(PRO) | {"TP53"}
    ok = [{c for c, n in zip(d["cats"], np.bincount(d["codes"], minlength=len(d["cats"]))) if n >= MIN_CELLS} for d in info.values()]
    perts = sorted((ok[0] & ok[1]) - set(CTRL_NAMES) - excl)
    common = sorted(set(info[lx]["genes"]) & set(info[ly]["genes"]))
    print("perturbations", len(perts), "common genes", len(common), flush=True)
    H, texp = {}, {}
    for seed, lab in enumerate((lx, ly)):
        d = info[lab]; ctrl = [c for c in CTRL_NAMES if c in d["cats"]][0]
        # target expression in this line: control mean log1p CP10K over the line's own gene list (absent -> 0)
        crow = np.where(d["codes"] == d["cats"].index(ctrl))[0]
        own = {g: j for j, g in enumerate(d["genes"])}; tg = [t for t in perts if t in own]
        Xc = read_rows(d["path"], crow, np.array([own[t] for t in tg]))
        e = np.log1p(Xc / d["umi"][crow][:, None] * 1e4).mean(0); em = dict(zip(tg, e))
        texp[lab] = np.array([em.get(t, 0.0) for t in perts])
        keep = {d["cats"].index(c) for c in set(perts) | {ctrl}}
        rows = np.where(np.isin(d["codes"], list(keep)))[0]
        X = read_rows(d["path"], rows, np.array([own[g] for g in common]))
        L = Line(X, d["codes"][rows], d["cats"], common, d["gem"][rows], ctrl); del X
        rng = np.random.default_rng(100 + seed); mu = L.Y[L.ctrl].mean(0)
        Hh = np.zeros((2, len(perts), len(common)))
        for k, t in enumerate(perts):
            r = rng.permutation(np.where(L.codes == L.cats.index(t))[0]); h = len(r) // 2
            Hh[0, k], Hh[1, k] = L.Y[r[:h]].mean(0) - mu, L.Y[r[h:2 * h]].mean(0) - mu
            if t in L.g2j:
                Hh[:, k, L.g2j[t]] = 0
        H[lab] = (Hh, mu); del L
        print(lab, "done", round(time.time() - t0), "s", flush=True)
    gk = np.where(np.minimum(H[lx][1], H[ly][1]) >= MIN_EXPR)[0]
    C = {lab: H[lab][0][..., gk] for lab in H}
    sets = {lab: {"inert": texp[lab] < 0.01, "expressed": texp[lab] >= 0.5} for lab in (lx, ly)}
    res = {"rules": __doc__, "labels": [lx, ly], "n_perturbations": len(perts), "n_genes": int(len(gk)),
           "n_inert": {lab: int(sets[lab]["inert"].sum()) for lab in sets},
           "inert_targets": {lab: [perts[i] for i in np.where(sets[lab]["inert"])[0]] for lab in sets}}
    res.update(analyse(C[lx][0], C[lx][1], C[ly][0], C[ly][1], sets, (lx, ly)))
    return res


def run_cd4(dr, ds, out):
    from cd4_readout import load, MIN_HALF, MIN_EXPR as CMIN
    (gr, genes_r, Sr, nr), (gs, genes_s, Ss, ns) = load(dr), load(ds)
    cr, cs = gr.index("__ctrl__"), gs.index("__ctrl__")
    common = sorted(set(genes_r) & set(genes_s)); pr_, ps_ = {g: j for j, g in enumerate(genes_r)}, {g: j for j, g in enumerate(genes_s)}
    colr, cols = np.array([pr_[g] for g in common]), np.array([ps_[g] for g in common])
    ctr = np.stack([np.asarray(Sr[cr, h])[colr] / nr[cr, h] for h in range(2)]); cts = np.stack([np.asarray(Ss[cs, h])[cols] / ns[cs, h] for h in range(2)])
    cpos = {g: j for j, g in enumerate(common)}
    gk = np.where((ctr.mean(0) >= CMIN) & (cts.mean(0) >= CMIN))[0]; genes = [common[j] for j in gk]
    ir, is_ = {g: i for i, g in enumerate(gr)}, {g: i for i, g in enumerate(gs)}
    targets = sorted(t for t in gr if t != "__ctrl__" and t in is_ and nr[ir[t]].min() >= MIN_HALF and ns[is_[t]].min() >= MIN_HALF and t in cpos)
    P, G = len(targets), len(gk)
    XR = np.zeros((2, P, G)); XS = np.zeros((2, P, G))
    for k, t in enumerate(targets):
        for h in range(2):
            XR[h, k] = np.asarray(Sr[ir[t], h])[colr[gk]] / nr[ir[t], h] - ctr[h, gk]
            XS[h, k] = np.asarray(Ss[is_[t], h])[cols[gk]] / ns[is_[t], h] - cts[h, gk]
    gpos = {g: j for j, g in enumerate(genes)}
    for k, t in enumerate(targets):
        if t in gpos:
            XR[:, k, gpos[t]] = 0; XS[:, k, gpos[t]] = 0
    er = np.array([ctr[:, cpos[t]].mean() for t in targets]); es = np.array([cts[:, cpos[t]].mean() for t in targets])
    sets = {"rest": {"inert": np.maximum(er, es) < 0.01, "expressed": er >= 0.5},       # inert in BOTH states (10-18 fix)
            "stim": {"inert": np.maximum(er, es) < 0.01, "expressed": es >= 0.5}}
    res = {"rules": __doc__, "labels": ["rest", "stim"], "n_perturbations": P, "n_genes": G,
           "n_inert": int(sets["rest"]["inert"].sum())}
    res.update(analyse(XR[0], XR[1], XS[0], XS[1], sets, ("rest", "stim")))
    return res


def main():
    mode = sys.argv[1]
    if mode == "h5ad":
        out = sys.argv[6]; res = run_h5ad(*sys.argv[2:7])
    else:
        out = sys.argv[4]; res = run_cd4(*sys.argv[2:5])
    os.makedirs(out, exist_ok=True)
    print(json.dumps({k: v for k, v in res.items() if k not in ("rules", "inert_targets")}), flush=True)
    json.dump(res, open(out + "/metrics_readout_inert_line.json", "w"), indent=1)
    print("saved", flush=True)


if __name__ == "__main__":
    main()
