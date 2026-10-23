"""Whole-brain in vivo Perturb-seq: one streaming pass over the h5ad shards (A100; numpy / h5py / scipy), keeping only
the sums the clone null test needs. Nuclei: wb_meta.npz (QC-passing, one guide, Glut or GABA; random half bit fixed
there). Expression: gene-expression features, log1p CP10K with the nucleus's own gene-expression total (over all
gene-expression features); kept genes = mean >= GMIN over the non-targeting nuclei of the first shard in either class
(the analysis applies its own, stricter threshold on the full sums).
Sums kept (float32):
  S[u, h, :]      per unit u = guide * 2 + class and half h: sum of log1p CP10K; n[u, h] nuclei
  M[k, h, :]      per non-targeting stratum k = mouse * 2 + class and half h: sum over non-targeting nuclei; m[k, h]
  Q[c, :]         per class: sum of squares over non-targeting nuclei (per-nucleus variance for scaling)
  nm[u, h, mouse] nuclei per unit, half and mouse (for references within mouse)
Usage: wb_aggregate.py <wb_meta.npz> <h5ad_dir> <out.npz> [workers]"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, time, glob
import numpy as np, h5py
from multiprocessing import Pool
from scipy.sparse import csr_matrix

CHUNK, GMIN = 4000, 0.05
G_ = {}


def dec(a):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])


def features(f):
    var = f["var"]; names = dec(var[var.attrs.get("_index", "_index")][:])
    ft = var["feature_types"]; cats = list(dec(ft["categories"][:])); codes = ft["codes"][:]
    return names, codes == cats.index("Gene Expression")


def chunks(f, sel, cols_expr):
    """yield (rows in sel, log1p CP10K over kept genes) for contiguous blocks of the shard"""
    N = f["X"].attrs["shape"][0]; nf = f["X"].attrs["shape"][1]
    for r0 in range(0, N, CHUNK):
        r1 = min(r0 + CHUNK, N); rr = sel[(sel >= r0) & (sel < r1)]
        if not len(rr):
            continue
        ip = f["X"]["indptr"][r0:r1 + 1]; a, b = int(ip[0]), int(ip[-1])
        X = csr_matrix((f["X"]["data"][a:b], f["X"]["indices"][a:b].astype(np.int64), ip - a), shape=(r1 - r0, nf))[rr - r0]
        Xe = X[:, G_["expr_cols"]].tocsr()
        tot = np.asarray(Xe.sum(1)).ravel()
        Xk = X[:, cols_expr].tocsr().astype(np.float32)
        Xk.data = np.log1p(Xk.data / np.repeat(np.maximum(tot, 1), np.diff(Xk.indptr)) * 1e4).astype(np.float32)
        yield rr, Xk


def add_rows(A, keys, X):
    uk, inv = np.unique(keys, return_inverse=True)
    R = csr_matrix((np.ones(len(keys), np.float32), (inv, np.arange(len(keys)))), shape=(len(uk), len(keys)))
    A[uk] += (R @ X).toarray()


def shard_pass(path):
    g = G_; t0 = time.time()
    with h5py.File(path, "r") as f:
        names, expr = features(f)
        assert list(names[g["keep_cols"]]) == g["genes"], "gene order differs in %s" % path
        ids = dec(f["obs"][f["obs"].attrs.get("_index", "_index")][:])
        rows = np.array([g["pos"].get(c, -1) for c in ids]); sel = np.where(rows >= 0)[0]
        G = len(g["genes"]); U2 = g["nU"] * 2; K2 = g["nK"] * 2
        S = np.zeros((U2, G), np.float32); M = np.zeros((K2, G), np.float32); Q = np.zeros((2, G), np.float32)
        n = np.zeros(U2, np.int64); m = np.zeros(K2, np.int64); nm = np.zeros((U2, g["nM"]), np.int32)
        for rr, X in chunks(f, sel, g["keep_cols"]):
            mi = rows[rr]; h = g["half"][mi]; c = g["cls"][mi]
            uk = (g["guide"][mi] * 2 + c) * 2 + h
            add_rows(S, uk, X); np.add.at(n, uk, 1); np.add.at(nm, (uk, g["mouse"][mi]), 1)
            nt = g["isnt"][mi]
            if nt.any():
                Xn = X[nt]; kk = (g["mouse"][mi] * 2 + c)[nt] * 2 + h[nt]
                add_rows(M, kk, Xn); np.add.at(m, kk, 1)
                Xq = Xn.multiply(Xn).tocsr(); cc = c[nt]
                for ci in (0, 1):
                    if (cc == ci).any():
                        Q[ci] += np.asarray(Xq[cc == ci].sum(0)).ravel()
    print(os.path.basename(path), "nuclei", len(sel), "%.0fs" % (time.time() - t0), flush=True)
    return S, n, M, m, Q, nm


def main():
    meta, hdir, outp = sys.argv[1:4]; nw = int(sys.argv[4]) if len(sys.argv) > 4 else 7
    z = np.load(meta, allow_pickle=True)   # our own wb_meta.npz (string columns stored as object arrays by pandas)
    z = {k: (z[k].astype(str) if z[k].dtype == object else z[k]) for k in z.files}
    gu, gidx = np.unique(z["guide"], return_inverse=True); mu, midx = np.unique(z["mouse"], return_inverse=True)
    tgt = dict(zip(z["guide"], z["target"]))
    paths = sorted(glob.glob(os.path.join(hdir, "*.h5ad")))
    pos = {c: i for i, c in enumerate(z["cell_id"])}; isnt = z["target"] == "Non_target"
    with h5py.File(paths[0], "r") as f:
        names, expr = features(f)
    expr_cols = np.where(expr)[0]
    G_.update(pos=pos, guide=gidx, cls=z["cls"].astype(np.int64), half=z["half"].astype(np.int64), mouse=midx, isnt=isnt,
              nU=len(gu) * 2, nK=len(mu) * 2, nM=len(mu), expr_cols=expr_cols)
    # gene set from the non-targeting nuclei of the first shard
    with h5py.File(paths[0], "r") as f:
        ids = dec(f["obs"][f["obs"].attrs.get("_index", "_index")][:])
        rows = np.array([pos.get(c, -1) for c in ids]); sel = np.where(rows >= 0)[0]; sel = sel[isnt[rows[sel]]]
        acc = np.zeros((2, len(expr_cols))); cnt = np.zeros(2)
        for rr, X in chunks(f, sel, expr_cols):
            c = G_["cls"][rows[rr]]
            for ci in (0, 1):
                acc[ci] += np.asarray(X[c == ci].sum(0)).ravel(); cnt[ci] += (c == ci).sum()
    mean = acc / np.maximum(cnt, 1)[:, None]
    keep_cols = expr_cols[(mean >= GMIN).any(0)]
    G_.update(keep_cols=keep_cols, genes=list(names[keep_cols]))
    print("genes kept", len(keep_cols), "of", len(expr_cols), "(non-targeting nuclei in first shard:", cnt.tolist(), ")", flush=True)
    t0 = time.time(); tot = None
    with Pool(nw) as pool:
        for res in pool.imap_unordered(shard_pass, paths):
            tot = list(res) if tot is None else [a + r for a, r in zip(tot, res)]
    S, n, M, m, Q, nm = tot
    U, K = len(gu) * 2, len(mu) * 2; G = S.shape[1]
    np.savez(outp, S=S.reshape(U, 2, G), n=n.reshape(U, 2), M=M.reshape(K, 2, G), m=m.reshape(K, 2), Q=Q,
             nm=nm.reshape(U, 2, len(mu)), guides=gu, targets=np.array([tgt[x] for x in gu]), mice=mu,
             genes=np.array(G_["genes"]), shards=np.array([os.path.basename(p) for p in paths]))
    print("saved", outp, "units", U, "genes", G, "nuclei", int(n.sum()), "of", len(z["cell_id"]), "%.0fs" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
