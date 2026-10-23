"""Model pilot (notes/model_pilot_shortlist.md): cell subset of CD4 donor 1, resting and 48 h stimulated, for
context-transfer models scored against same-clone and independent-replicate truth (A100; numpy / h5py / scipy).
Targets
  test   the 300 targets of the cell-eval comparison (inert_names, strong_names of review_tests/celleval subset),
         less those without MINC cells per guide in both states (recorded as dropped_test)
  train  NTRAIN further targets at random (seed 5) with both guides (gp 0 / 1) >= MINC cells in both states
Cells (per state, from cells_<state>_v2/cells.npz; umi > 0): per (target, gp) up to CAP_TEST (test) or CAP_TRAIN
  (train) at random (seed 6); NCTRL non-targeting control cells, each with its gp.
Genes: NHVG genes with control mean log1p CP10K >= 0.05 in both states, ranked by the sum over states of the
  standardised variance of training-target means (from sum8 aggregates), plus the measured test target genes.
Output per state: <out>/pilot_<state>.npz with raw counts CSR (cells x genes), umi over ALL genes, target, gp, kind
  (control / train / inert / strong), gene names. Usage: cd4_pilot_extract.py <sc_dir> <celleval_npz> <out> [workers]
With env PILOT_ALL_GENES=1 the same cells (same seeds and selection) are written with ALL genes, as
pilot_<state>_allgenes.npz with the h5ad row index (for Stack and State Embedding, which need the full gene axis);
pilot_design.json is then not rewritten. Donor: env PILOT_DONOR (default D1); <sc_dir> must hold cells_<donor>_<state>_v2/
and <donor>_<state>.assigned_guide.h5ad (model_eval step 1 links them per donor)."""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json, time
import numpy as np
import h5py
from multiprocessing import Pool
from scipy.sparse import csr_matrix, vstack

DONOR = os.environ.get("PILOT_DONOR", "D1")
STATES = (DONOR + "_Rest", DONOR + "_Stim48hr")
NTRAIN, MINC, CAP_TEST, CAP_TRAIN, NCTRL, NHVG, CHUNK = 1000, 40, 200, 80, 6000, 2000, 4000
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


def load_state(sc, st):
    d = "%s/cells_%s_v2" % (sc, st)
    meta = json.load(open(d + "/meta.json"))
    S = np.load(d + "/sum8.npy", mmap_mode="r"); n8 = np.load(d + "/n8.npy")
    return dict(genes=meta["genes"], groups=meta["groups"], S=S, n8=n8, cells=np.load(d + "/cells.npz"))


def main():
    sc, cev, out = sys.argv[1:4]; nw = int(sys.argv[4]) if len(sys.argv) > 4 else 24; t0 = time.time()
    os.makedirs(out, exist_ok=True)
    c = np.load(cev, allow_pickle=False)
    inert, strong = [str(x) for x in c["inert_names"]], [str(x) for x in c["strong_names"]]
    test = set(inert) | set(strong)
    D = {st: load_state(sc, st) for st in STATES}
    genes = D[STATES[0]]["genes"]   # reference order; the states hold the same genes in different orders
    gpos = {g: j for j, g in enumerate(genes)}
    for st in STATES:
        assert sorted(D[st]["genes"]) == sorted(genes), "gene sets differ between states"
        sp = {g: j for j, g in enumerate(D[st]["genes"])}
        D[st]["col"] = np.array([sp[g] for g in genes])   # reference -> state column
    # targets with both guides >= MINC in both states
    ok = None
    for st in STATES:
        d = D[st]; ci = len(d["groups"]) - 1
        ng0, ng1 = d["n8"][:, :4].sum(1), d["n8"][:, 4:].sum(1)
        s = {d["groups"][i] for i in range(ci) if ng0[i] >= MINC and ng1[i] >= MINC}
        ok = s if ok is None else ok & s
    dropped = sorted(test - ok)   # too few stimulated cells in one guide; kept out of the pilot and recorded
    inert, strong = [t for t in inert if t in ok], [t for t in strong if t in ok]; test = set(inert) | set(strong)
    rng = np.random.default_rng(5)
    train = sorted(rng.choice(sorted(ok - test), NTRAIN, replace=False).tolist())
    # genes: expressed in both states, ranked by variance of training-target means
    score, expr = np.zeros(len(genes)), np.ones(len(genes), bool)
    for st in STATES:
        d = D[st]; ci = len(d["groups"]) - 1; gi = {g: i for i, g in enumerate(d["groups"])}
        cm = (np.asarray(d["S"][ci].sum(0)) / d["n8"][ci].sum())[d["col"]]; expr &= cm >= 0.05
        M = np.stack([(np.asarray(d["S"][gi[t]].sum(0)) / d["n8"][gi[t]].sum())[d["col"]] for t in train])
        v = M.var(0); v = (v - v[expr].mean()) / v[expr].std(); score += v
    score[~expr] = -np.inf
    hvg = np.argsort(-score)[:NHVG]
    gsel = sorted(set(hvg.tolist()) | {gpos[t] for t in test if t in gpos})
    allg = os.environ.get("PILOT_ALL_GENES") == "1"
    if allg:
        gsel = list(range(len(genes)))
    else:
      json.dump(dict(train=train, inert=inert, strong=strong, dropped_test=dropped, n_genes=len(gsel), states=STATES,
                   caps=dict(test=CAP_TEST, train=CAP_TRAIN, ctrl=NCTRL, minc=MINC)),
              open(out + "/pilot_design.json", "w"), indent=1)
    print("targets ok", len(ok), "train", len(train), "test", len(test), "genes", len(gsel), "%.0fs" % (time.time() - t0), flush=True)
    rng = np.random.default_rng(6)
    for st in STATES:
        d = D[st]; ci = len(d["groups"]) - 1; gi = {g: i for i, g in enumerate(d["groups"])}
        grp, gp, umi = d["cells"]["grp"], d["cells"]["gp"], d["cells"]["umi"]
        sel = np.zeros(len(grp), bool); kind = np.full(len(grp), "", dtype="<U7")
        for names, k, cap in ((train, "train", CAP_TRAIN), (inert, "inert", CAP_TEST), (strong, "strong", CAP_TEST)):
            for t in names:
                for p in (0, 1):
                    r = np.where((grp == gi[t]) & (gp == p) & (umi > 0))[0]
                    pick = rng.choice(r, min(cap, len(r)), replace=False); sel[pick] = True; kind[pick] = k
        r = np.where((grp == ci) & (umi > 0))[0]; pick = rng.choice(r, min(NCTRL, len(r)), replace=False)
        sel[pick] = True; kind[pick] = "control"
        fn = "%s/%s.assigned_guide.h5ad" % (sc, st)
        with h5py.File(fn, "r") as f:
            N, Gn = f["X"].attrs["shape"]
        assert int(N) == len(grp), "cells.npz does not match %s" % fn
        assert int(Gn) == len(genes), "h5ad gene count differs from meta"
        G_.update(fn=fn, sel=sel, cols=d["col"][np.array(gsel)], G=int(Gn))
        edges = np.linspace(0, int(N), nw + 1).astype(int); parts = []
        with Pool(nw) as pool:
            for res in pool.imap_unordered(worker, [(edges[w], edges[w + 1]) for w in range(nw)]):
                if res is not None:
                    parts.append(res)
        parts.sort(key=lambda x: x[0][0])
        idx = np.concatenate([p[0] for p in parts]); X = vstack([p[1] for p in parts]).tocsr(); U = np.concatenate([p[2] for p in parts])
        tname = np.array(["non-targeting" if grp[k] == ci else d["groups"][grp[k]] for k in idx])
        np.savez(out + "/pilot_%s%s.npz" % (st, "_allgenes" if allg else ""), data=X.data.astype(np.float32), indices=X.indices,
                 indptr=X.indptr, shape=np.array(X.shape), umi=U.astype(np.float32), target=tname, gp=gp[idx], kind=kind[idx],
                 genes=np.array([genes[j] for j in gsel]), row=idx)
        print(st, "cells", X.shape[0], "genes", X.shape[1], "%.0fs" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
