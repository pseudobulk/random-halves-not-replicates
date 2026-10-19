"""Non-targeting constructs as null perturbations in cell-line screens (clone preprint, findings cc, cf; review v11 M1).
Replogle et al. 2022 and Nadig et al. 2025 raw single-cell files (obs: gene, sgID_AB; categorical columns stored with
obs/__categories), and the VCC 2025 H1 hESC training data (obs: target_gene, guide_id). One streaming pass: log1p CP10K with row sums over all genes; a random bit per cell (numpy
default_rng(21) over all cells in file order). Aggregates per (group, bit): non-targeting cells grouped by construct
(sgID_AB), targeting cells grouped by target gene.
  genes: mean >= 0.1 over non-targeting cells; groups: >= 10 cells in each half
  non-targeting: effect = construct mean minus the mean over kept constructs in that half (null baseline);
    random halves vs pairs of different constructs (first half of one, second half of another; 50 random pairings)
  targeting (reference): effect = gene mean minus the mean of all non-targeting cells in that half; random halves
  ratio = mean non-targeting random-half cross-half signal / mean targeting random-half cross-half signal
Usage: cl_ntc_null.py <out.json> <label>:<h5ad> [...]   (workers from env WORKERS, default 16)"""
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


def column(f, k):
    o = f["obs"]
    if isinstance(o[k], h5py.Group):                                   # anndata >= 0.8 categorical
        return o[k]["codes"][:], dec(o[k]["categories"][:])
    codes = o[k][:]
    if "__categories" in o and k in o["__categories"]:                 # anndata 0.7 categorical
        return codes, dec(o["__categories"][k][:])
    vals = dec(codes); cats = sorted(set(vals)); pos = {c: i for i, c in enumerate(cats)}
    return np.array([pos[x] for x in vals]), cats


def read_rows(f, r0, r1):
    X = f["X"]
    if isinstance(X, h5py.Dataset):
        return csr_matrix(X[r0:r1].astype(np.float32))
    ip = X["indptr"][r0:r1 + 1]; a, b = int(ip[0]), int(ip[-1])
    M = csr_matrix((X["data"][a:b].astype(np.float32), X["indices"][a:b].astype(np.int32), ip - a), shape=(r1 - r0, G_["G"]))
    M.sum_duplicates(); return M


def worker(args):
    wid, r0w, r1w = args
    g = G_; acc = np.zeros((g["NG"] * 2, g["G"]), np.float32); cnt = np.zeros(g["NG"] * 2, np.int64)
    with h5py.File(g["fn"], "r") as f:
        for r0 in range(r0w, r1w, CHUNK):
            r1 = min(r0 + CHUNK, r1w); k = g["grp"][r0:r1]
            if not (k >= 0).any():
                continue
            X = read_rows(f, r0, r1); umi = np.asarray(X.sum(1)).ravel()
            v = np.where((k >= 0) & (umi > 0))[0]
            X, umi = X[v], umi[v]
            L = X.copy(); L.data = np.log1p(L.data / np.repeat(umi, np.diff(L.indptr)) * 1e4)
            idx = k[v] * 2 + g["bit"][r0:r1][v]
            R = (csr_matrix((np.ones(len(v), np.float32), (idx, np.arange(len(v)))), shape=(g["NG"] * 2, len(v))) @ L).tocoo()
            acc[R.row, R.col] += R.data; cnt += np.bincount(idx, minlength=g["NG"] * 2)
    return wid, acc, cnt


def extract(fn, nw):
    with h5py.File(fn, "r") as f:
        cols = set(f["obs"].keys())
        gk, sk = ("gene", "sgID_AB") if "sgID_AB" in cols else ("target_gene", "guide_id")    # Replogle/Nadig or VCC 2025 H1
        gene, gcats = column(f, gk); sg, scats = column(f, sk)
        X = f["X"]; N, Gn = (X.shape if isinstance(X, h5py.Dataset) else tuple(X.attrs["shape"]))
    nt = np.array(["non-targeting" in c.lower() for c in gcats])[gene]
    tgt_names = [c for c in gcats if "non-targeting" not in c.lower()]
    tpos = {c: i for i, c in enumerate(tgt_names)}
    gmap = np.array([tpos.get(c, -1) for c in gcats])
    nt_used = np.unique(sg[nt]); spos = -np.ones(len(scats), np.int64); spos[nt_used] = np.arange(len(nt_used))
    T = len(tgt_names); K = len(nt_used)
    grp = np.where(nt, T + spos[sg], gmap[gene]).astype(np.int64)
    bit = np.random.default_rng(21).integers(0, 2, int(N)).astype(np.int64)
    G_.update(fn=fn, G=int(Gn), NG=T + K, grp=grp, bit=bit)
    edges = np.linspace(0, int(N), nw + 1).astype(int); acc = None; cnt = 0
    with Pool(nw) as pool:
        for wid, a, c in pool.imap_unordered(worker, [(w, edges[w], edges[w + 1]) for w in range(nw)]):
            acc = a if acc is None else acc + a; cnt = cnt + c
    S = acc.reshape(T + K, 2, int(Gn)); n = cnt.reshape(T + K, 2)
    return S[:T], n[:T], S[T:], n[T:], int(N), [scats[j] for j in nt_used]


def metrics(St, nt_, Sn, nn, rng):
    tot = Sn.sum((0, 1)) / nn.sum(); gk = tot >= 0.1
    ctrl = Sn.sum(0)[:, gk] / nn.sum(0)[:, None]                         # all non-targeting cells, per half
    okn = (nn >= 10).all(1); Mn = Sn[okn][:, :, gk] / nn[okn][:, :, None]
    Xn = Mn - Mn.mean(0, keepdims=True); rnd = (Xn[:, 0] * Xn[:, 1]).sum(1)
    K = len(Xn); fp, mp = [], []
    for _ in range(50):
        p = rng.permutation(K); a, b = p[: K // 2], p[K // 2: 2 * (K // 2)]
        c = (Xn[a, 0] * Xn[b, 1]).sum(1); fp.append(float((c > 0).mean())); mp.append(float(c.mean()))
    boot = [float((rnd[rng.integers(0, K, K)] > 0).mean()) for _ in range(1000)]
    okt = (nt_ >= 10).all(1); Mt = St[okt][:, :, gk] / nt_[okt][:, :, None]
    Xt = Mt - ctrl[None]; Xt = Xt - Xt.mean(0, keepdims=True)
    tc = (Xt[:, 0] * Xt[:, 1]).sum(1)
    return dict(n_genes=int(gk.sum()), n_constructs=int(K), cells_per_half_median=float(np.median(nn[okn])),
                random_frac_pos=float((rnd > 0).mean()), random_frac_pos_ci95=[float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
                random_mean_cross=float(rnd.mean()), pairs_frac_pos_mean=float(np.mean(fp)), pairs_frac_pos_range=[min(fp), max(fp)],
                pairs_mean_cross=float(np.mean(mp)), n_targets=int(okt.sum()), target_cells_per_half_median=float(np.median(nt_[okt])),
                target_random_mean_cross=float(tc.mean()), target_random_median_cross=float(np.median(tc)),
                ratio_nt_to_target_mean=float(rnd.mean() / tc.mean()))


def main():
    out = {"rules": __doc__}; nw = int(os.environ.get("WORKERS", "16"))
    for arg in sys.argv[2:]:
        lab, fn = arg.split(":", 1); t0 = time.time()
        St, nt_, Sn, nn, N, names = extract(fn, nw)
        r = metrics(St, nt_, Sn, nn, np.random.default_rng(0)); r.update(file=os.path.basename(fn), n_cells=N, seconds=round(time.time() - t0))
        out[lab] = r; print(lab, json.dumps(r), flush=True)
        json.dump(out, open(sys.argv[1], "w"), indent=1)


if __name__ == "__main__":
    main()
