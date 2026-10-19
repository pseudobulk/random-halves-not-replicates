"""Per-cell extraction for the clone tests of the preprint review (A100; CD4 single cells, Zhu / Marson 2026).
One streaming pass over a donor-condition h5ad. Cells as cd4_sc_extract.py (low_quality and unassigned cells dropped;
controls = guide_type non-targeting, own group NG - 1; log1p CP10K with row sums over all genes). Each cell gets
three bits:
  gp  parity of its guide's rank within its target (controls: within the non-targeting set); guides are independent
      transductions, so one clone falls in one gp only
  r1, r2  independent random bits (seed 21)
and the aggregates are kept per (group, gp, r1, r2), so that one pass gives
  random halves (r1), random quarters (r1 within r2 = 0), guide halves (gp), guide quarters (gp within r2 = 0).
Per cell it also keeps the raw counts of TRBC1, TRBC2 and TRAC. Every T-cell clone carries one TCR beta rearrangement,
transcribed with either constant region TRBC1 or TRBC2, so the TRBC1 fraction is a 1-bit clone marker (the probe
panel has no V genes).
Outputs: sum8.npy (NG x 8 x G float32; index gp * 4 + r1 * 2 + r2), n8.npy (NG x 8), cells.npz (grp, gid, gp, r1, r2,
umi, trbc1, trbc2, trac), meta.json (genes, groups, guide categories).
Usage: cd4_sc_extract_cells.py <h5ad> <out_dir> [workers]"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json, time
import numpy as np
import h5py
from multiprocessing import Pool
from scipy.sparse import csr_matrix

CHUNK = 4000
G_ = {}


def dec(a):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in a]


def cat(o, k):
    return o[k]["codes"][:], dec(o[k]["categories"][:])


def read_rows(f, r0, r1):
    ip = f["X"]["indptr"][r0:r1 + 1]
    a, b = int(ip[0]), int(ip[-1])
    X = csr_matrix((f["X"]["data"][a:b].astype(np.float32), f["X"]["indices"][a:b].astype(np.int32), ip - a),
                   shape=(r1 - r0, G_["G"]))
    X.sum_duplicates()
    return X


def name_parity(names):
    """guide parity from the guide's own name ("GENE-1" -> 1, "GENE-2" -> 0, "NTC-802" -> 0), identical across files
    (10-18 fix: parity from category-code rank differed between Rest and Stim for 558 targets)"""
    out = np.zeros(len(names), np.int64)
    for i, x in enumerate(names):
        a = x.rsplit("-", 1)
        out[i] = int(a[1]) % 2 if len(a) == 2 and a[1].isdigit() else 0
    return out


def worker(args):
    wid, r0w, r1w = args
    g = G_; Gn, NG = g["G"], g["n_grp"]
    acc = np.zeros((NG * 8, Gn), np.float32)
    n = r1w - r0w
    tr = np.zeros((n, 3), np.float32); um = np.zeros(n, np.float32)
    with h5py.File(g["fn"], "r") as f:
        for r0 in range(r0w, r1w, CHUNK):
            r1 = min(r0 + CHUNK, r1w); o = r0 - r0w
            X = read_rows(f, r0, r1)
            umi = np.asarray(X.sum(1)).ravel()
            grp = g["grp"][r0:r1].copy(); grp[umi <= 0] = -1
            v = np.where(grp >= 0)[0]
            if not len(v):
                continue
            X, umi, grp = X[v], umi[v], grp[v]
            tr[o + v] = X[:, g["trcols"]].toarray(); um[o + v] = umi
            L = X.copy(); L.data = np.log1p(L.data / np.repeat(umi, np.diff(L.indptr)) * 1e4)
            idx = grp * 8 + g["code8"][r0:r1][v]
            R = (csr_matrix((np.ones(len(v), np.float32), (idx, np.arange(len(v)))), shape=(acc.shape[0], len(v))) @ L).tocoo()
            acc[R.row, R.col] += R.data
    np.save("%s/part8_%03d.npy" % (g["out"], wid), acc)
    return wid, r0w, tr, um


def main():
    fn, out = sys.argv[1], sys.argv[2]
    nw = int(sys.argv[3]) if len(sys.argv) > 3 else 16
    os.makedirs(out, exist_ok=True); t0 = time.time()
    with h5py.File(fn, "r") as f:
        o = f["obs"]
        pc, pcats = cat(o, "perturbed_gene_name"); gc, gcats = cat(o, "guide_type"); gid, gidcats = cat(o, "guide_id")
        lq = o["low_quality"][:]
        genes = dec(f["var"]["gene_name"][:]); N, Gn = f["X"].attrs["shape"]
    ctrl = np.isin(gc, [k for k, x in enumerate(gcats) if "non" in x.lower()])
    NG = len(pcats) + 1
    grp = np.where(ctrl, NG - 1, pc).astype(np.int64)
    grp[(pc < 0) & ~ctrl] = -1; grp[lq.astype(bool)] = -1
    gid = gid.astype(np.int64); par = name_parity(gidcats)
    gp = np.where(gid >= 0, par[np.maximum(gid, 0)], 0).astype(np.int64)
    rng = np.random.default_rng(21); r1 = rng.integers(0, 2, int(N)); r2 = rng.integers(0, 2, int(N))
    code8 = gp * 4 + r1 * 2 + r2
    gpos = {x: k for k, x in enumerate(genes)}
    trcols = np.array([gpos[x] for x in ("TRBC1", "TRBC2", "TRAC")])
    G_.update(fn=fn, out=out, G=int(Gn), n_grp=NG, grp=grp, code8=code8, trcols=trcols)
    edges = np.linspace(0, int(N), nw + 1).astype(int)
    tr = np.zeros((int(N), 3), np.float32); um = np.zeros(int(N), np.float32)
    with Pool(nw) as pool:
        for wid, r0, t, u in pool.imap_unordered(worker, [(w, edges[w], edges[w + 1]) for w in range(nw)]):
            tr[r0:r0 + len(u)] = t; um[r0:r0 + len(u)] = u
            print("worker %d done %.0fs" % (wid, time.time() - t0), flush=True)
    acc = None
    for w in range(nw):
        p = np.load("%s/part8_%03d.npy" % (out, w))
        acc = p if acc is None else acc + p
        os.remove("%s/part8_%03d.npy" % (out, w))
    ok = (grp >= 0) & (um > 0)
    n8 = np.bincount(grp[ok] * 8 + code8[ok], minlength=NG * 8).reshape(NG, 8)
    np.save(out + "/sum8.npy", acc.reshape(NG, 8, int(Gn))); np.save(out + "/n8.npy", n8)
    np.savez(out + "/cells.npz", grp=grp, gid=gid, gp=gp, r1=r1, r2=r2, umi=um, trbc1=tr[:, 0], trbc2=tr[:, 1], trac=tr[:, 2])
    json.dump(dict(file=fn, n_cells=int(N), n_kept=int(ok.sum()), genes=genes, groups=pcats + ["__ctrl__"],
                   guide_ids=gidcats, rules=__doc__), open(out + "/meta.json", "w"))
    print("saved %s, kept %d cells, %.0fs" % (out, ok.sum(), time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
