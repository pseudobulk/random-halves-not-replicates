"""Model evaluation step 3: Frangieh et al. 2021 (Perturb-CITE-seq, patient-derived melanoma cells, Cas9 knockout;
scPerturb v1.4 FrangiehIzar2021_RNA.h5ad) in the pilot format, for a low-clonality point. One transduced pool was
expanded and plated into three conditions (Control, IFN-gamma, co-culture with TILs), so clones are shared across
conditions, but each sgRNA's cells come from many clones (cell line, ~1,000 cells per sgRNA at transduction).
Cells: single-guide cells (MOI == 1, sgRNA called). UMI = obs ncounts (all genes). Pair: source Control -> target T
(IFNg or Coculture). Units = sgRNAs; a unit qualifies with >= MINC cells in both conditions.
  gene     test targets: genes with >= 2 qualifying sgRNAs; the two with most cells (min over conditions) are g0 / g1
           (by name). Further qualifying sgRNAs of the gene become training units '<gene>~<k>' (k = 3, ...), so
           baselines and STATE see the same training set.
  null     the non-targeting sgRNAs that qualify, shuffled (seed 7); half are paired into pseudo-targets 'NTP<i>' (g0 / g1),
           whose A - B can only come from clones; all other non-targeting cells are the control group.
  train    every other sgRNA with >= MINC_TRAIN cells in both conditions (gp = name rank within its gene).
Cross-fit (about 250 genes only): the gene and null test targets are split into NFOLD folds (seed 8); fold k tests its
own targets and trains on everything else (the other fold's test targets included). Same cells for every fold.
Cells per unit: up to CAP (seed 6); control cells up to NCTRL. Genes: control mean log1p CP10K >= 0.05 in both
conditions; the NHVG with the largest standardised variance of unit means in the SOURCE condition (all units; no
target-condition data of test units is used), plus the measured test target genes.
Output: <out>/<T>/pilot_Control.npz, pilot_<T>.npz (pilot format + unit = sgRNA name), fold<k>/pilot_design.json
(extra_kinds null, gene) with links to the npz. Usage: me_frangieh_extract.py <h5ad> <out_dir>"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "8")
import sys, json, time
import numpy as np
import h5py
from scipy.sparse import csc_matrix

MINC, MINC_TRAIN, CAP, NCTRL, NHVG, NFOLD = 30, 10, 200, 6000, 2000, 2
TARGETS = {"IFNg": "IFNγ", "Coculture": "Co-culture"}


def cat(o, k):
    c = [x.decode() if isinstance(x, bytes) else str(x) for x in o[k]["categories"][:]]
    codes = o[k]["codes"][:]
    return np.array([c[i] if i >= 0 else "nan" for i in codes])


def main():
    fn, out = sys.argv[1:3]; t0 = time.time()
    with h5py.File(fn, "r") as f:
        o = f["obs"]
        sg, pert, cond = cat(o, "sgRNA"), cat(o, "perturbation"), cat(o, "perturbation_2")
        moi, umi = o["MOI"][:], o["ncounts"][:].astype(np.float32)
        g = f["var"]["gene_symbol"]
        genes = [x.decode() if isinstance(x, bytes) else str(x) for x in (g[:] if isinstance(g, h5py.Dataset) else np.array(g["categories"][:])[g["codes"][:]])]
        X = f["X"]; shape = tuple(X.attrs["shape"])
        X = csc_matrix((X["data"][:].astype(np.float32), X["indices"][:], X["indptr"][:]), shape=shape)
    keep = (moi == 1) & (sg != "nan")
    rows = np.where(keep)[0]
    X = X.tocsr()[rows]; sg, pert, cond, umi = sg[rows], pert[rows], cond[rows], umi[rows]
    print("single-guide cells", len(rows), "sgRNAs", len(set(sg)), "%.0fs" % (time.time() - t0), flush=True)
    L = X.copy(); L.data = np.log1p(L.data / np.repeat(umi, np.diff(L.indptr)) * 1e4)
    gpos = {x: j for j, x in enumerate(genes)}
    nt = {s for s, p in zip(sg, pert) if p == "control"}
    for short, T in TARGETS.items():
        od = os.path.join(out, short); os.makedirs(od, exist_ok=True)
        S = ("Control", T)
        n = {c: {} for c in S}
        for s, c in zip(sg, cond):
            if c in n:
                n[c][s] = n[c].get(s, 0) + 1
        cnt = {s: min(n[S[0]].get(s, 0), n[S[1]].get(s, 0)) for s in set(sg)}
        qual = {s for s, v in cnt.items() if v >= MINC}
        by_gene = {}
        for s in sorted(set(sg) - nt):
            by_gene.setdefault(pert[np.where(sg == s)[0][0]], []).append(s)
        unit = {}   # sgRNA -> (target, gp, kind)
        test_genes = []
        for gname, ss in sorted(by_gene.items()):
            q = sorted([s for s in ss if s in qual], key=lambda s: -cnt[s])
            if len(q) >= 2:
                test_genes.append(gname)
                for p, s in enumerate(sorted(q[:2])):
                    unit[s] = (gname, p, "gene")
                for k, s in enumerate(q[2:]):
                    unit[s] = ("%s~%d" % (gname, k + 3), 0, "train")
                rest = [s for s in ss if s not in q and cnt[s] >= MINC_TRAIN]
                for k, s in enumerate(sorted(rest)):
                    unit[s] = ("%s~%d" % (gname, len(q) + k + 1), 0, "train")
            else:
                for p, s in enumerate(sorted(s for s in ss if cnt[s] >= MINC_TRAIN)):
                    unit[s] = (gname, p, "train")
        ntq = sorted(s for s in nt if s in qual)
        rng = np.random.default_rng(7); rng.shuffle(ntq)
        npair = len(ntq) // 4   # half of the qualifying non-targeting sgRNAs, in pairs
        null = []
        for i in range(npair):
            nm = "NTP%02d" % i; null.append(nm)
            for p, s in enumerate(sorted(ntq[2 * i:2 * i + 2])):
                unit[s] = (nm, p, "null")
        ctrl_sg = nt - {s for s in ntq[:2 * npair]}
        rng = np.random.default_rng(8)
        fold = {t: i % NFOLD for i, t in enumerate(rng.permutation(test_genes).tolist())}
        fold.update({t: i % NFOLD for i, t in enumerate(rng.permutation(null).tolist())})
        # cells
        rng = np.random.default_rng(6); files = {}
        for c in S:
            idx, tg, gp, kd, un = [], [], [], [], []
            for s, (t, p, k) in sorted(unit.items()):
                r = np.where((sg == s) & (cond == c))[0]
                pick = np.sort(rng.choice(r, min(CAP, len(r)), replace=False))
                idx += pick.tolist(); tg += [t] * len(pick); gp += [p] * len(pick); kd += [k] * len(pick); un += [s] * len(pick)
            r = np.where(np.isin(sg, sorted(ctrl_sg)) & (cond == c))[0]
            pick = np.sort(rng.choice(r, min(NCTRL, len(r)), replace=False))
            idx += pick.tolist(); tg += ["non-targeting"] * len(pick); gp += [0] * len(pick); kd += ["control"] * len(pick); un += list(sg[pick])
            files[c] = dict(idx=np.array(idx), target=np.array(tg), gp=np.array(gp), kind=np.array(kd), unit=np.array(un))
        # genes
        expr = np.ones(len(genes), bool)
        for c in S:
            m = files[c]["kind"] == "control"
            expr &= np.asarray(L[files[c]["idx"][m]].mean(0)).ravel() >= 0.05
        f0 = files[S[0]]; labs = np.array(["%s|%d" % (t, p) for t, p in zip(f0["target"], f0["gp"])])
        ul = sorted(set(labs[f0["kind"] != "control"]))
        M = np.stack([np.asarray(L[f0["idx"][labs == u]].mean(0)).ravel() for u in ul])
        v = M.var(0); v = (v - v[expr].mean()) / v[expr].std(); v[~expr] = -np.inf
        gsel = sorted(set(np.argsort(-v)[:NHVG].tolist()) | {gpos[t] for t in test_genes if t in gpos})
        for c, short_c in zip(S, ("Control", short)):
            d = files[c]; Xs = X[d["idx"]][:, gsel].tocsr()
            np.savez(os.path.join(od, "pilot_%s.npz" % short_c), data=Xs.data.astype(np.float32), indices=Xs.indices,
                     indptr=Xs.indptr, shape=np.array(Xs.shape), umi=umi[d["idx"]], target=d["target"], gp=d["gp"],
                     kind=d["kind"], genes=np.array([genes[j] for j in gsel]), row=rows[d["idx"]], unit=d["unit"])
            print(short, c, "cells", Xs.shape[0], "genes", Xs.shape[1], "%.0fs" % (time.time() - t0), flush=True)
        alltest = set(test_genes) | set(null)
        train_t = sorted({t for t, p, k in unit.values()} - alltest)
        for k in range(NFOLD):
            fd = os.path.join(od, "fold%d" % k); os.makedirs(fd, exist_ok=True)
            for c in ("Control", short):
                lk = os.path.join(fd, "pilot_%s.npz" % c)
                if not os.path.lexists(lk):
                    os.symlink(os.path.join("..", "pilot_%s.npz" % c), lk)
            gk = sorted(t for t in test_genes if fold[t] == k); nk = sorted(t for t in null if fold[t] == k)
            tr = sorted(train_t + [t for t in alltest if fold[t] != k])
            json.dump(dict(train=tr, inert=[], strong=[], null=nk, gene=gk, extra_kinds=["null", "gene"], dropped_test=[],
                           n_genes=len(gsel), states=["Control", short], fold=k, nfold=NFOLD,
                           units={s: list(map(str, u)) for s, u in unit.items()}, control_sgrnas=sorted(ctrl_sg),
                           caps=dict(test=CAP, train=CAP, ctrl=NCTRL, minc=MINC, minc_train=MINC_TRAIN), rules=__doc__),
                      open(os.path.join(fd, "pilot_design.json"), "w"), indent=1)
        print(short, "test genes", len(test_genes), "null pairs", len(null), "train targets", len(train_t),
              "control sgRNAs", len(ctrl_sg), "genes", len(gsel), "%.0fs" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
