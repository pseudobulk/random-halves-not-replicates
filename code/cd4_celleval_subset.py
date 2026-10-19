"""Cell subset for running cell-eval itself on CD4 donor-1 Rest (preprint review 2.5; A100).
Targets (from cells_D1_Rest aggregates, cd4_sc_extract_cells.py): both guide halves (gp = 0 / 1) with >= MINC cells;
  inert    target measured and control mean log1p CP10K < 0.01: NI chosen at random (seed 4)
  strong   control mean >= 0.5, both guide halves knock the target down (Delta log1p <= KD) and the weaker guide reaches
           >= RATIO of the stronger: the NS strongest by mean knockdown
Cells: up to CAP per target (random), and NCTRL non-targeting control cells (random); each cell keeps its guide
parity gp (controls: parity within the non-targeting set), so the guide split can be built downstream.
Genes: control mean log1p CP10K >= 0.01 (plus all target genes that are measured).
Output npz: raw counts CSR (cells x genes), per-cell UMI over ALL genes (for CP10K), target name, gp, gene names.
Usage: cd4_celleval_subset.py <h5ad> <cells_dir> <out_npz> [workers]. Exploratory."""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json, time
import numpy as np
import h5py
from multiprocessing import Pool
from scipy.sparse import csr_matrix, vstack

NI, NS, MINC, CAP, NCTRL, KD, RATIO, CHUNK = 150, 150, 40, 400, 15000, -0.25, 0.7, 4000
G_ = {}


def worker(args):
    r0w, r1w = args
    g = G_; sel = g["sel"]; cols = g["cols"]
    rows = np.where(sel[r0w:r1w])[0] + r0w
    out, umis, idx = [], [], []
    with h5py.File(g["fn"], "r") as f:
        for r0 in range(r0w, r1w, CHUNK):
            r1 = min(r0 + CHUNK, r1w)
            rr = rows[(rows >= r0) & (rows < r1)]
            if not len(rr):
                continue
            ip = f["X"]["indptr"][r0:r1 + 1]; a, b = int(ip[0]), int(ip[-1])
            X = csr_matrix((f["X"]["data"][a:b].astype(np.float32), f["X"]["indices"][a:b].astype(np.int32), ip - a), shape=(r1 - r0, g["G"]))
            X.sum_duplicates(); X = X[rr - r0]
            umis.append(np.asarray(X.sum(1)).ravel()); out.append(X[:, cols]); idx.append(rr)
    if not out:
        return None
    return np.concatenate(idx), vstack(out).tocsr(), np.concatenate(umis)


def main():
    fn, cd, outp = sys.argv[1:4]; nw = int(sys.argv[4]) if len(sys.argv) > 4 else 16; t0 = time.time()
    meta = json.load(open(cd + "/meta.json")); genes, groups = meta["genes"], meta["groups"]
    S = np.load(cd + "/sum8.npy", mmap_mode="r"); n8 = np.load(cd + "/n8.npy"); ci = len(groups) - 1
    cm = np.asarray(S[ci].sum(0)) / n8[ci].sum(); gpos = {g: j for j, g in enumerate(genes)}
    ng0, ng1 = n8[:, :4].sum(1), n8[:, 4:].sum(1)
    cand = [i for i in range(ci) if ng0[i] >= MINC and ng1[i] >= MINC and groups[i] in gpos]
    inert = [i for i in cand if cm[gpos[groups[i]]] < 0.01]
    kd = {}
    for i in cand:
        j = gpos[groups[i]]
        if cm[j] < 0.5:
            continue
        k0 = S[i, :4, j].sum() / ng0[i] - S[ci, :4, j].sum() / n8[ci, :4].sum()
        k1 = S[i, 4:, j].sum() / ng1[i] - S[ci, 4:, j].sum() / n8[ci, 4:].sum()
        strong, weak = min(k0, k1), max(k0, k1)
        if weak <= KD and weak / strong >= RATIO:
            kd[i] = (k0 + k1) / 2
    rng = np.random.default_rng(4)
    pick_i = sorted(rng.choice(inert, min(NI, len(inert)), replace=False).tolist())
    pick_s = sorted(sorted(kd, key=lambda i: kd[i])[:NS])
    print("candidates", len(cand), "inert", len(inert), "strong matched", len(kd), "picked", len(pick_i), len(pick_s), flush=True)
    c = np.load(cd + "/cells.npz"); grp, gp, umi = c["grp"], c["gp"], c["umi"]
    sel = np.zeros(len(grp), bool)
    for i in pick_i + pick_s:
        r = np.where((grp == i) & (umi > 0))[0]
        sel[rng.choice(r, min(CAP, len(r)), replace=False)] = True
    r = np.where((grp == ci) & (umi > 0))[0]; sel[rng.choice(r, min(NCTRL, len(r)), replace=False)] = True
    gsel = sorted(set(np.where(cm >= 0.01)[0].tolist()) | {gpos[groups[i]] for i in pick_i + pick_s})
    with h5py.File(fn, "r") as f:
        N, Gn = f["X"].attrs["shape"]
    G_.update(fn=fn, sel=sel, cols=np.array(gsel), G=int(Gn))
    edges = np.linspace(0, int(N), nw + 1).astype(int)
    parts = []
    with Pool(nw) as pool:
        for res in pool.imap_unordered(worker, [(edges[w], edges[w + 1]) for w in range(nw)]):
            if res is not None:
                parts.append(res)
    parts.sort(key=lambda x: x[0][0])
    idx = np.concatenate([p[0] for p in parts]); X = vstack([p[1] for p in parts]).tocsr(); U = np.concatenate([p[2] for p in parts])
    tname = np.array(["non-targeting" if grp[k] == ci else groups[grp[k]] for k in idx])
    pis = set(pick_i)
    kind = np.array(["control" if grp[k] == ci else ("inert" if grp[k] in pis else "strong") for k in idx])
    np.savez(outp, data=X.data.astype(np.float32), indices=X.indices, indptr=X.indptr, shape=np.array(X.shape),
             umi=U.astype(np.float32), target=tname, kind=kind, gp=gp[idx], genes=np.array([genes[j] for j in gsel]),
             strong_kd=np.array([kd[i] for i in pick_s]), strong_names=np.array([groups[i] for i in pick_s]),
             inert_names=np.array([groups[i] for i in pick_i]))
    print("saved", outp, "cells", X.shape[0], "genes", X.shape[1], "%.0fs" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
