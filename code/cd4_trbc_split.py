"""TRBC halves: a direct clone test in the CD4 screen (A100; Zhu / Marson CD4 single cells). Every T-cell clone carries
one TCR beta rearrangement; cells whose transcripts are all TRBC2 come from clones joined at D-J2 on their expressed
allele(s), cells with a high TRBC1 share from clones with a D-J1 joint. Cells are classed per cell from the TRBC1 and
TRBC2 counts kept by cd4_sc_extract_cells.py (cells.npz, same row order):
  class 0  TRBC2 only, >= MIN0 counts
  class 1  TRBC1 fraction >= F1, >= MIN1 counts
  other cells are not used. A TRBC1 clone can show only TRBC2 counts by chance (for a true fraction 0.5 and 4 counts,
  in 6% of cells), so a few cells of one clone can fall in the other class; this works against the test.
Splitting a perturbation's classed cells into class 0 against class 1 gives two halves with the same guides, donor,
library and state that share (nearly) no clone. The same classed cells split by a random bit give random halves.
Aggregates per (group, class, r1, r2) as in cd4_sc_extract_cells.py (r1, r2 random bits, seed 31).
Usage: cd4_trbc_split.py <h5ad> <cells.npz> <out_dir> [workers]"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json, time
import numpy as np
import h5py
from multiprocessing import Pool
from scipy.sparse import csr_matrix

CHUNK, MIN0, MIN1, F1 = 4000, 4, 2, 0.5
G_ = {}


def dec(a):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in a]


def read_rows(f, r0, r1):
    ip = f["X"]["indptr"][r0:r1 + 1]; a, b = int(ip[0]), int(ip[-1])
    X = csr_matrix((f["X"]["data"][a:b].astype(np.float32), f["X"]["indices"][a:b].astype(np.int32), ip - a),
                   shape=(r1 - r0, G_["G"]))
    X.sum_duplicates()
    return X


def worker(args):
    wid, r0w, r1w = args
    g = G_; NG = g["n_grp"]
    acc = np.zeros((NG * 8, g["G"]), np.float32)
    with h5py.File(g["fn"], "r") as f:
        for r0 in range(r0w, r1w, CHUNK):
            r1 = min(r0 + CHUNK, r1w)
            grp = g["grp"][r0:r1]; v = np.where(grp >= 0)[0]
            if not len(v):
                continue
            X = read_rows(f, r0, r1)[v]; umi = np.asarray(X.sum(1)).ravel()
            keep = umi > 0; X, v, umi = X[keep], v[keep], umi[keep]
            L = X.copy(); L.data = np.log1p(L.data / np.repeat(umi, np.diff(L.indptr)) * 1e4)
            idx = grp[v] * 8 + g["code"][r0:r1][v]
            R = (csr_matrix((np.ones(len(v), np.float32), (idx, np.arange(len(v)))), shape=(acc.shape[0], len(v))) @ L).tocoo()
            acc[R.row, R.col] += R.data
    np.save("%s/part_%03d.npy" % (g["out"], wid), acc)
    return wid


def main():
    fn, cz, out = sys.argv[1:4]; nw = int(sys.argv[4]) if len(sys.argv) > 4 else 16
    os.makedirs(out, exist_ok=True); t0 = time.time()
    z = np.load(cz); grp0 = z["grp"].astype(np.int64); t1, t2 = z["trbc1"], z["trbc2"]
    with h5py.File(fn, "r") as f:
        genes = dec(f["var"]["gene_name"][:]); N, Gn = f["X"].attrs["shape"]
        o = f["obs"]; pcats = dec(o["perturbed_gene_name"]["categories"][:])
    assert len(grp0) == N
    tot = t1 + t2; frac = np.where(tot > 0, t1 / np.maximum(tot, 1), np.nan)
    cls = np.full(int(N), -1, np.int64)
    cls[(t1 == 0) & (tot >= MIN0)] = 0
    cls[(frac >= F1) & (tot >= MIN1)] = 1
    rng = np.random.default_rng(31); r1 = rng.integers(0, 2, int(N)); r2 = rng.integers(0, 2, int(N))
    grp = np.where(cls >= 0, grp0, -1)
    code = np.maximum(cls, 0) * 4 + r1 * 2 + r2
    NG = int(grp0.max()) + 1
    G_.update(fn=fn, out=out, G=int(Gn), n_grp=NG, grp=grp, code=code)
    edges = np.linspace(0, int(N), nw + 1).astype(int)
    with Pool(nw) as pool:
        for wid in pool.imap_unordered(worker, [(w, edges[w], edges[w + 1]) for w in range(nw)]):
            print("worker %d %.0fs" % (wid, time.time() - t0), flush=True)
    acc = None
    for w in range(nw):
        p = np.load("%s/part_%03d.npy" % (out, w)); acc = p if acc is None else acc + p
        os.remove("%s/part_%03d.npy" % (out, w))
    ok = grp >= 0
    n8 = np.bincount(grp[ok] * 8 + code[ok], minlength=NG * 8).reshape(NG, 8)
    np.save(out + "/sum8.npy", acc.reshape(NG, 8, int(Gn))); np.save(out + "/n8.npy", n8)
    json.dump(dict(file=fn, n_cells=int(N), class0=int(((cls == 0) & (grp0 >= 0)).sum()), class1=int(((cls == 1) & (grp0 >= 0)).sum()),
                   genes=genes, groups=pcats + ["__ctrl__"], rules=__doc__), open(out + "/meta.json", "w"))
    print("saved", out, "classed cells", int(ok.sum()), "%.0fs" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
