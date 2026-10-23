"""Per-cell aggregates with 16 codes for the depth grid of Fig. 2d (A100; CD4 single cells; same cells, guide parity and
random bits r1, r2 (seed 21) as cd4_sc_extract_cells.py, plus a third random bit r3 (seed 22)).
code16 = gp * 8 + r1 * 4 + r2 * 2 + r3, so the four random bits r1 r2 r3 split each guide's cells into eight equal parts
and both guides enter every random side. Genes: those kept by clone_tests_cd4.py (control mean >= 0.1 in the existing
extraction of the same file, TCR genes removed) to bound memory.
Outputs: sum16.npy (NG x 16 x Gk float32), n16.npy (NG x 16), meta16.json (genes kept, groups).
Usage: cd4_sc_extract_grid.py <h5ad> <existing_cells_dir> <out_dir> [workers]"""
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
    acc = np.zeros((NG * 16, len(g["gk"])), np.float32)
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
            um[o + v] = umi
            L = X.copy(); L.data = np.log1p(L.data / np.repeat(umi, np.diff(L.indptr)) * 1e4); L = L[:, g["gk"]]
            idx = grp * 16 + g["code16"][r0:r1][v]
            R = (csr_matrix((np.ones(len(v), np.float32), (idx, np.arange(len(v)))), shape=(acc.shape[0], len(v))) @ L).tocoo()
            acc[R.row, R.col] += R.data
    np.save("%s/part16_%03d.npy" % (g["out"], wid), acc)
    return wid, r0w, tr, um


def main():
    fn, old, out = sys.argv[1], sys.argv[2], sys.argv[3]
    nw = int(sys.argv[4]) if len(sys.argv) > 4 else 8
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
    r3 = np.random.default_rng(22).integers(0, 2, int(N))
    code16 = gp * 8 + r1 * 4 + r2 * 2 + r3
    om = json.load(open(old + "/meta.json")); assert om["genes"] == genes and om["n_cells"] == int(N)
    S8 = np.load(old + "/sum8.npy", mmap_mode="r"); n8o = np.load(old + "/n8.npy")
    cm = np.asarray(S8[-1].sum(0), dtype=np.float64) / n8o[-1].sum()
    gk = np.array([j for j in np.where(cm >= 0.1)[0] if not genes[j].startswith(("TRBC", "TRAC", "TRBV", "TRAV", "TRGC", "TRDC"))])
    gpos = {x: k for k, x in enumerate(genes)}
    trcols = np.array([gpos[x] for x in ("TRBC1", "TRBC2", "TRAC")])
    G_.update(fn=fn, out=out, G=int(Gn), n_grp=NG, grp=grp, code16=code16, trcols=trcols, gk=gk)
    edges = np.linspace(0, int(N), nw + 1).astype(int)
    tr = np.zeros((int(N), 3), np.float32); um = np.zeros(int(N), np.float32)
    with Pool(nw) as pool:
        for wid, r0, t, u in pool.imap_unordered(worker, [(w, edges[w], edges[w + 1]) for w in range(nw)]):
            tr[r0:r0 + len(u)] = t; um[r0:r0 + len(u)] = u
            print("worker %d done %.0fs" % (wid, time.time() - t0), flush=True)
    acc = None
    for w in range(nw):
        p = np.load("%s/part16_%03d.npy" % (out, w))
        acc = p if acc is None else acc + p
        os.remove("%s/part16_%03d.npy" % (out, w))
    ok = (grp >= 0) & (um > 0)
    n16 = np.bincount(grp[ok] * 16 + code16[ok], minlength=NG * 16).reshape(NG, 16)
    assert (n16.sum(1) == n8o.sum(1)).all(), "cells differ from the existing extraction"   # same cells per group; guide parity may differ (name-based since 10-18)
    np.save(out + "/sum16.npy", acc.reshape(NG, 16, len(gk))); np.save(out + "/n16.npy", n16)
    json.dump(dict(file=fn, n_cells=int(N), n_kept=int(ok.sum()), genes=[genes[j] for j in gk], groups=pcats + ["__ctrl__"],
                   rules=__doc__), open(out + "/meta16.json", "w"))
    print("saved %s, kept %d cells, %.0fs" % (out, ok.sum(), time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
