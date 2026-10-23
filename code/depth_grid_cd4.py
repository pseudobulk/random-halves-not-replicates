"""Depth grid for Fig. 2d: cross-half signal and Pearson agreement at 1/8, 1/4, 3/8 and 1/2 of a target's cells per side
(A100; output of cd4_sc_extract_grid.py; estimator as clone_tests_cd4.py DEPTH).
code16 = gp * 8 + r, r = r1 * 4 + r2 * 2 + r3 in 0..7 (eight equal random parts of each guide's cells).
  random side  both guides, r in a set S of m parts (m / 8 of the cells); A and B use disjoint S
  guide side   one guide, r in a set T of 2m parts (also m / 8 of the cells); A = guide 0, B = guide 1, same T
Several disjoint choices are averaged at each m. Effect = side mean - control mean of the same codes (log1p CP10K; genes
with control mean >= 0.1, TCR genes removed; the target's own gene set to 0); cross-half signal referenced to the mean of
the inert targets. Inert (control mean < 0.01) and expressed (>= 0.5) targets are classified from the full gene list of
the existing extraction (cd4_sc_extract_cells.py). Targets: >= MINS cells on every side of every design.
Usage: depth_grid_cd4.py <cells_dir> <grid_dir> <label> [<cells_dir> <grid_dir> <label>] <out_json>."""
import sys, json
import numpy as np

MINS = 4
RAND = {1: [([0], [1]), ([2], [3]), ([4], [5]), ([6], [7])], 2: [([0, 1], [2, 3]), ([4, 5], [6, 7])],
        3: [([0, 1, 2], [4, 5, 6]), ([1, 2, 3], [5, 6, 7])], 4: [([0, 1, 2, 3], [4, 5, 6, 7]), ([0, 1, 4, 5], [2, 3, 6, 7])]}
GUIDE = {1: [[0, 1], [2, 3], [4, 5], [6, 7]], 2: [[0, 1, 2, 3], [4, 5, 6, 7]], 3: [[0, 1, 2, 3, 4, 5], [2, 3, 4, 5, 6, 7]],
         4: [list(range(8))]}
DES = {"random": {m: [(a + [8 + r for r in a], b + [8 + r for r in b]) for a, b in p] for m, p in RAND.items()},
       "guide": {m: [(list(t), [8 + r for r in t]) for t in ts] for m, ts in GUIDE.items()}}


def pearson_rows(A, B):
    a = A - A.mean(1, keepdims=True); b = B - B.mean(1, keepdims=True)
    return (a * b).sum(1) / np.sqrt((a ** 2).sum(1) * (b ** 2).sum(1))


def run_state(cdir, gdir, label):
    om = json.load(open(cdir + "/meta.json")); S8 = np.load(cdir + "/sum8.npy", mmap_mode="r"); n8 = np.load(cdir + "/n8.npy")
    cm_all = np.asarray(S8[-1].sum(0), dtype=np.float64) / n8[-1].sum()
    texp = {g: cm_all[j] for j, g in enumerate(om["genes"])}
    meta = json.load(open(gdir + "/meta16.json")); genes, groups = meta["genes"], meta["groups"]
    S = np.load(gdir + "/sum16.npy", mmap_mode="r"); n16 = np.load(gdir + "/n16.npy"); ci = len(groups) - 1
    gpos = {g: j for j, g in enumerate(genes)}
    ok = np.ones(ci, bool)
    for des in DES.values():
        for prs in des.values():
            for a, b in prs:
                ok &= (n16[:ci, a].sum(1) >= MINS) & (n16[:ci, b].sum(1) >= MINS)
    tids = [i for i in range(ci) if ok[i] and groups[i] in texp]
    tn = [groups[i] for i in tids]
    inert = np.array([texp[t] < 0.01 for t in tn]); expd = np.array([texp[t] >= 0.5 for t in tn])
    Sg = np.asarray(S[tids], dtype=np.float64); Sc = np.asarray(S[ci], dtype=np.float64)
    nt = n16[tids].astype(float); nc = n16[ci].astype(float)
    own = [(k, gpos[t]) for k, t in enumerate(tn) if t in gpos]

    def side(codes):
        m = Sg[:, codes].sum(1) / nt[:, codes].sum(1)[:, None] - Sc[codes].sum(0) / nc[codes].sum()
        for k, j in own:
            m[k, j] = 0
        return m

    res = {"label": label, "n_targets": len(tids), "n_inert": int(inert.sum()), "n_expressed": int(expd.sum()), "grid": {}}
    for des, byf in DES.items():
        rows = []
        for m, prs in byf.items():
            acc = {"inert": [], "expressed": []}; cells = []
            for a, b in prs:
                A, B = side(a), side(b)
                Ai, Bi = A - A[inert].mean(0), B - B[inert].mean(0)
                cross = (Ai * Bi).sum(1); pr = pearson_rows(A, B)
                cells += [np.median(nt[:, a].sum(1)), np.median(nt[:, b].sum(1))]
                for s, mm in (("inert", inert), ("expressed", expd)):
                    acc[s].append((cross[mm].mean(), (cross[mm] > 0).mean(), np.nanmean(pr[mm])))
            r = {"fraction": m / 8, "n_choices": len(prs), "median_cells_per_side": float(np.mean(cells))}
            for s in acc:
                v = np.array(acc[s]).mean(0)
                r[s] = dict(cross_per_target=float(v[0]), frac_pos=float(v[1]), pearson_mean=float(v[2]))
            rows.append(r); print(label, des, json.dumps(r), flush=True)
        res["grid"][des] = rows
    return res


def main():
    args = sys.argv[1:]; out = args[-1]
    res = {"rules": __doc__}
    for i in range(0, len(args) - 1, 3):
        res[args[i + 2]] = run_state(args[i], args[i + 1], args[i + 2])
    json.dump(res, open(out, "w"), indent=1)
    print("saved", flush=True)


if __name__ == "__main__":
    main()
