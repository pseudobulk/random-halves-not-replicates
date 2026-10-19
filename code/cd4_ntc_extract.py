"""Non-targeting guides as null perturbations (clone preprint v9, findings ca; A100; CD4 single cells, Zhu / Marson 2026).
One streaming pass over a donor-condition h5ad, keeping only non-targeting cells (guide_type non-targeting; low_quality
and zero-UMI cells dropped), log1p CP10K with row sums over all genes, as cd4_sc_extract_cells.py. Each cell gets the
same random bits as that extraction (numpy default_rng(21): r1, then r2, over all cells in file order), so random halves
here are the r1 halves used throughout. Aggregates are kept per (non-targeting guide, r1).
Outputs: sum2.npy (n_guides x 2 x G float32), n2.npy (n_guides x 2), meta.json (genes, guide names, rules).
Usage: cd4_ntc_extract.py <h5ad> <out_dir> [workers]"""
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


def worker(args):
    wid, r0w, r1w = args
    g = G_; Gn, K = g["G"], g["K"]
    acc = np.zeros((K * 2, Gn), np.float32); cnt = np.zeros(K * 2, np.int64)
    with h5py.File(g["fn"], "r") as f:
        for r0 in range(r0w, r1w, CHUNK):
            r1 = min(r0 + CHUNK, r1w)
            k = g["k"][r0:r1]
            if not (k >= 0).any():
                continue
            X = read_rows(f, r0, r1)
            umi = np.asarray(X.sum(1)).ravel()
            v = np.where((k >= 0) & (umi > 0))[0]
            if not len(v):
                continue
            X, umi = X[v], umi[v]
            L = X.copy(); L.data = np.log1p(L.data / np.repeat(umi, np.diff(L.indptr)) * 1e4)
            idx = k[v] * 2 + g["r1"][r0:r1][v]
            R = (csr_matrix((np.ones(len(v), np.float32), (idx, np.arange(len(v)))), shape=(K * 2, len(v))) @ L).tocoo()
            acc[R.row, R.col] += R.data
            cnt += np.bincount(idx, minlength=K * 2)
    np.save("%s/part2_%03d.npy" % (g["out"], wid), acc)
    return wid, cnt


def main():
    fn, out = sys.argv[1], sys.argv[2]
    nw = int(sys.argv[3]) if len(sys.argv) > 3 else 16
    os.makedirs(out, exist_ok=True); t0 = time.time()
    with h5py.File(fn, "r") as f:
        o = f["obs"]
        gc, gcats = cat(o, "guide_type"); gid, gidcats = cat(o, "guide_id")
        lq = o["low_quality"][:]
        genes = dec(f["var"]["gene_name"][:]); N, Gn = f["X"].attrs["shape"]
    ctrl = np.isin(gc, [j for j, x in enumerate(gcats) if "non" in x.lower()]) & ~lq.astype(bool) & (gid >= 0)
    rng = np.random.default_rng(21); r1 = rng.integers(0, 2, int(N)); r2 = rng.integers(0, 2, int(N))   # same bits as cd4_sc_extract_cells.py
    used = np.unique(gid[ctrl]); K = len(used)
    pos = -np.ones(len(gidcats), np.int64); pos[used] = np.arange(K)
    k = np.where(ctrl, pos[np.maximum(gid, 0)], -1).astype(np.int64)
    G_.update(fn=fn, out=out, G=int(Gn), K=K, k=k, r1=r1.astype(np.int64))
    edges = np.linspace(0, int(N), nw + 1).astype(int)
    cnt = np.zeros(K * 2, np.int64)
    with Pool(nw) as pool:
        for wid, c in pool.imap_unordered(worker, [(w, edges[w], edges[w + 1]) for w in range(nw)]):
            cnt += c; print("worker %d done %.0fs" % (wid, time.time() - t0), flush=True)
    acc = None
    for w in range(nw):
        p = np.load("%s/part2_%03d.npy" % (out, w))
        acc = p if acc is None else acc + p
        os.remove("%s/part2_%03d.npy" % (out, w))
    np.save(out + "/sum2.npy", acc.reshape(K, 2, int(Gn))); np.save(out + "/n2.npy", cnt.reshape(K, 2))
    json.dump(dict(file=fn, n_cells=int(N), n_ctrl=int(ctrl.sum()), genes=genes, guides=[gidcats[j] for j in used],
                   rules=__doc__), open(out + "/meta.json", "w"))
    print("saved %s, %d guides, %d cells, %.0fs" % (out, K, cnt.sum(), time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
