"""Per-perturbation arrays for the preprint figures (A100). From the per-cell aggregates of cd4_sc_extract_cells.py
(audit-fixed version), for one donor (Rest and Stim48hr):
  per target and state: inert / expressed flags (inert = measured and control mean < 0.01 in both states), cell counts,
  cross-half signal x_A . x_B referenced to the inert mean, and per-target Pearson between halves, for random halves
  (r1) and guide halves (gp); and the TRBC1 fraction deviation (mean f - f0) of each half for both designs, plus the
  overdispersion z of each target and of one label permutation.
Genes: control mean >= 0.1 in both states, TCR genes excluded. Output: <out>.npz
Usage: export_fig_data.py <cells_Rest> <cells_Stim> <out_prefix>"""
import sys, json
import numpy as np

DES = {"random": ([0, 1, 4, 5], [2, 3, 6, 7]), "guide": ([0, 1, 2, 3], [4, 5, 6, 7])}


def main():
    dr, ds, out = sys.argv[1:4]
    st = {}
    for k, d in (("rest", dr), ("stim", ds)):
        m = json.load(open(d + "/meta.json")); st[k] = (m["genes"], m["groups"], np.load(d + "/sum8.npy", mmap_mode="r"), np.load(d + "/n8.npy"), d)
    genes = sorted(set(st["rest"][0]) & set(st["stim"][0])); gpos = {g: j for j, g in enumerate(genes)}
    cms, ok = {}, None
    for k, (gl, groups, S, n8, d) in st.items():
        gi = {g: j for j, g in enumerate(gl)}; col = np.array([gi[g] for g in genes]); ci = len(groups) - 1
        cms[k] = np.asarray(S[ci].sum(0))[col] / n8[ci].sum()
        s_ = {t for i, t in enumerate(groups[:-1]) if all(n8[i, c].sum() >= 10 for a, b in DES.values() for c in (a, b))}
        ok = s_ if ok is None else ok & s_
    targets = sorted(t for t in ok if t in gpos)
    tr = {j for j, g in enumerate(genes) if g.startswith(("TRBC", "TRAC", "TRBV", "TRAV", "TRGC", "TRDC"))}
    keep = np.array([j for j in np.where(np.minimum(cms["rest"], cms["stim"]) >= 0.1)[0] if j not in tr]); kpos = {j: q for q, j in enumerate(keep)}
    e = np.array([max(cms["rest"][gpos[t]], cms["stim"][gpos[t]]) for t in targets]); inert = e < 0.01
    expd = np.array([min(cms["rest"][gpos[t]], cms["stim"][gpos[t]]) for t in targets]) >= 0.5
    zc = np.array([kpos.get(gpos[t], -1) for t in targets]); rows = np.arange(len(targets))
    res = dict(targets=np.array(targets), inert=inert, expressed=expd)
    for k, (gl, groups, S, n8, d) in st.items():
        gi = {g: j for j, g in enumerate(gl)}; col = np.array([gi[g] for g in genes])[keep]; ci = len(groups) - 1
        ti = {t: i for i, t in enumerate(groups)}; idx = np.array([ti[t] for t in targets])
        Sg = np.asarray(S[idx])[:, :, col].astype(np.float64); nt = n8[idx].astype(float); Sc = np.asarray(S[ci])[:, col].astype(np.float64); nc = n8[ci].astype(float)
        res["ncells_%s" % k] = nt.sum(1)
        for dn, (a, b) in DES.items():
            sides = []
            for codes in (a, b):
                X = Sg[:, codes].sum(1) / nt[:, codes].sum(1)[:, None] - Sc[codes].sum(0) / nc[codes].sum()
                m = zc >= 0; X[rows[m], zc[m]] = 0; sides.append(X)
            A, B = sides; Ai, Bi = A - A[inert].mean(0), B - B[inert].mean(0)
            res["cross_%s_%s" % (dn, k)] = (Ai * Bi).sum(1)
            a_, b_ = A - A.mean(1, keepdims=True), B - B.mean(1, keepdims=True)
            res["pearson_%s_%s" % (dn, k)] = (a_ * b_).sum(1) / np.sqrt((a_ ** 2).sum(1) * (b_ ** 2).sum(1))
        # TRBC
        c = np.load(d + "/cells.npz"); grp, gp, r1 = c["grp"], c["gp"], c["r1"]; t1, t2 = c["trbc1"], c["trbc2"]
        okc = (grp >= 0) & (c["umi"] > 0) & (t1 + t2 >= 2); f = np.where(t1 + t2 > 0, t1 / np.maximum(t1 + t2, 1), np.nan)
        f0, sdf = float(np.nanmean(f[okc])), float(np.nanstd(f[okc])); G = len(groups)
        def gm(mask):
            s = np.bincount(grp[mask], weights=f[mask], minlength=G); n = np.bincount(grp[mask], minlength=G); return s / np.maximum(n, 1), n
        mf, nf = gm(okc)
        res["trbc_z_%s" % k] = (mf[idx] - f0) / (sdf / np.sqrt(np.maximum(nf[idx], 1))); res["trbc_n_%s" % k] = nf[idx]
        rng = np.random.default_rng(3); pg = rng.permutation(grp[okc]); s = np.bincount(pg, weights=f[okc], minlength=G); n = np.bincount(pg, minlength=G)
        res["trbc_zperm_%s" % k] = ((s / np.maximum(n, 1)) - f0)[idx] / (sdf / np.sqrt(np.maximum(n[idx], 1)))
        for dn, (ma, mb) in (("random", (r1 == 0, r1 == 1)), ("guide", (gp == 0, gp == 1))):
            (fa, na), (fb, nb) = gm(okc & ma), gm(okc & mb)
            res["trbc_dA_%s_%s" % (dn, k)] = fa[idx] - f0; res["trbc_dB_%s_%s" % (dn, k)] = fb[idx] - f0
            res["trbc_nA_%s_%s" % (dn, k)] = na[idx]; res["trbc_nB_%s_%s" % (dn, k)] = nb[idx]
        print(k, "done", flush=True)
    np.savez_compressed(out + ".npz", **res)
    print("saved", out + ".npz", len(targets), "targets", int(inert.sum()), "inert", flush=True)


if __name__ == "__main__":
    main()
