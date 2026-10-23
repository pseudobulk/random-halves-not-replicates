"""Cells for the Shesha check (A100; CD4 donor 1 resting): log1p CP10K expression of every cell of the unexpressed
targets, of 400 randomly drawn expressed targets (seed 5) and of 20,000 randomly drawn non-targeting cells (the
controls), on the 2,000 genes with the highest control mean among genes with control mean >= 0.1 (TCR genes excluded).
Cells, targets, sgRNA parity and the gene filter are those of cd4_sc_extract_cells.py and clone_tests_cd4.py.
Output: one .npz with the CSR matrix, the label per cell (target name or 'control'), the sgRNA parity, the target
class (unexpressed / expressed / control) and the genes. Read with shesha_cd4.py (5090, with anndata and Shesha).
Usage: shesha_extract_cd4.py <h5ad> <cells_dir of the same file> <out.npz> [workers]"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json
import numpy as np
import h5py
from multiprocessing import Pool
from scipy.sparse import csr_matrix, vstack
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cd4_sc_extract_cells import dec, cat, name_parity

CHUNK = 4000
G_ = {}


def read_rows(f, r0, r1, G):
    ip = f["X"]["indptr"][r0:r1 + 1]
    a, b = int(ip[0]), int(ip[-1])
    X = csr_matrix((f["X"]["data"][a:b].astype(np.float32), f["X"]["indices"][a:b].astype(np.int32), ip - a), shape=(r1 - r0, G))
    X.sum_duplicates()
    return X


def worker(args):
    wid, r0w, r1w = args
    g = G_; keep = g["keep"]
    parts, rows = [], []
    with h5py.File(g["fn"], "r") as f:
        for r0 in range(r0w, r1w, CHUNK):
            r1 = min(r0 + CHUNK, r1w)
            v = np.where(keep[r0:r1])[0]
            if not len(v):
                continue
            X = read_rows(f, r0, r1, g["G"])[v]
            umi = np.asarray(X.sum(1)).ravel(); ok = umi > 0
            X, v, umi = X[ok], v[ok], umi[ok]
            L = X.copy(); L.data = np.log1p(L.data / np.repeat(umi, np.diff(L.indptr)) * 1e4)
            parts.append(L[:, g["gk"]].tocsr()); rows.append(r0 + v)
    return wid, (vstack(parts).tocsr() if parts else None), (np.concatenate(rows) if rows else np.array([], int))


def main():
    fn, cdir, out = sys.argv[1:4]
    nw = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    with h5py.File(fn, "r") as f:
        o = f["obs"]
        pc, pcats = cat(o, "perturbed_gene_name"); gc, gcats = cat(o, "guide_type"); gid, gidcats = cat(o, "guide_id")
        lq = o["low_quality"][:].astype(bool)
        genes = dec(f["var"]["gene_name"][:]); N, Gn = f["X"].attrs["shape"]
    ctrl = np.isin(gc, [k for k, x in enumerate(gcats) if "non" in x.lower()])
    gid = gid.astype(np.int64); par = name_parity(gidcats)
    gp = np.where(gid >= 0, par[np.maximum(gid, 0)], 0).astype(np.int64)
    # gene filter and target classes from the existing extraction of the same file
    om = json.load(open(cdir + "/meta.json")); assert om["genes"] == genes
    S = np.load(cdir + "/sum8.npy", mmap_mode="r"); n8 = np.load(cdir + "/n8.npy")
    cm = np.asarray(S[-1].sum(0), dtype=np.float64) / n8[-1].sum()
    tcr = ("TRBC", "TRAC", "TRBV", "TRAV", "TRGC", "TRDC")
    cand = np.array([j for j in np.where(cm >= 0.1)[0] if not genes[j].startswith(tcr)])
    gk = np.sort(cand[np.argsort(-cm[cand])[:2000]])
    gpos = {x: j for j, x in enumerate(genes)}
    texp = {t: cm[gpos[t]] for t in pcats if t in gpos}
    unexp = {k for k, t in enumerate(pcats) if t in texp and texp[t] < 0.01}
    expd = sorted(k for k, t in enumerate(pcats) if t in texp and texp[t] >= 0.5)
    rng = np.random.default_rng(5); expd = set(rng.choice(expd, 400, replace=False).tolist())
    base = (~lq) & (pc >= 0) & ~ctrl
    sel_t = base & np.isin(pc, list(unexp | expd))
    ci = np.where((~lq) & ctrl)[0]; sel_c = np.zeros(int(N), bool); sel_c[np.random.default_rng(6).choice(ci, 20000, replace=False)] = True
    keep = sel_t | sel_c
    G_.update(fn=fn, G=int(Gn), keep=keep, gk=gk)
    edges = np.linspace(0, int(N), nw + 1).astype(int)
    with Pool(nw) as pool:
        res = sorted(pool.map(worker, [(w, edges[w], edges[w + 1]) for w in range(nw)]))
    X = vstack([r[1] for r in res if r[1] is not None]).tocsr(); rows = np.concatenate([r[2] for r in res])
    lab = np.where(sel_c[rows], "control", np.array(pcats, dtype=object)[np.maximum(pc[rows], 0)]).astype(str)
    kind = np.where(sel_c[rows], "control", np.where(np.isin(pc[rows], list(unexp)), "unexpressed", "expressed")).astype(str)
    np.savez_compressed(out, data=X.data, indices=X.indices, indptr=X.indptr, shape=np.array(X.shape), label=lab, gp=gp[rows],
                        kind=kind, genes=np.array([genes[j] for j in gk]), rules=__doc__)
    print("saved", out, X.shape, {k: int((kind == k).sum()) for k in ("unexpressed", "expressed", "control")}, flush=True)


if __name__ == "__main__":
    main()
