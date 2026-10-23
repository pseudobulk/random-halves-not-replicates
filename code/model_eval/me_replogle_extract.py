"""Model evaluation, negative control (model_eval/plan.md section 7 step 2): the pilot design on two SEPARATE screens,
Replogle 2022 K562 essential (source) -> RPE1 essential (target), where no clone is shared between input and truth.
Writes the pilot format of cd4_pilot_extract.py, so pilot_common / pilot_baselines / pilot_state_prep / STATE run
unchanged (env PILOT_STATES=K562,RPE1).
Perturbation unit: the dual-sgRNA element (sgID_AB); a target's elements are separate transductions. Truth A for an
element = its own RPE1 cells (same element, separate screen); truth B = the RPE1 cells of the target's other element.
Without shared clones, A - B is only what the element itself carries (e.g. element-specific knockdown strength).
Targets
  test   targets with >= 2 elements of >= MINC cells in both screens; the two with most cells (min over screens) are
         kept, gp = their rank by sgID name. Classes as the CD4 pilot's cell-eval rule, required in both screens:
         inert   target gene not measured or control mean log1p CP10K < 0.01
         strong  control mean >= 0.5, both elements knock the target down (Delta log1p <= KD) and the weaker reaches
                 >= RATIO of the stronger (env RATIO, default 0.7)
         mismatched  the other two-element targets. The two elements of a Replogle target mostly sit at different
                 promoters (P1 / P2), and often only one knocks the target down, in both screens alike; here A - B > 0
                 is expected from the element's own effect, without clones. Kept as a separate group (extra_kinds).
         Every two-element target's knockdowns are written to the design (kd: per screen, element 0, element 1,
         control mean). Note: 'inert' here means not in the gene list of either file or control mean < 0.01; the files
         keep only detected genes, so some of these targets may be expressed at low level and have real effects.
  train  NTRAIN further targets at random (seed 5) with an element of >= MINC cells in both screens (the element with
         most cells; gp 0)
Cells: per element up to CAP_TEST (test) / CAP_TRAIN (train) at random (seed 6), NCTRL non-targeting cells per screen.
UMI = obs UMI_count (all genes of the original matrix). Genes: the genes measured in both files; NHVG expressed in
both (control mean >= 0.05), ranked by the summed standardised variance of training-target means, plus the measured
test target genes. Usage: me_replogle_extract.py <k562.h5ad> <rpe1.h5ad> <out_dir>"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "8")
import sys, json, time
import numpy as np
import h5py
from scipy.sparse import csr_matrix, vstack
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))   # explore_V_day6_day8.py in code/
from explore_V_day6_day8 import read_codes_cats, var_gene_names

STATES = ("K562", "RPE1")
NTRAIN, MINC, CAP_TEST, CAP_TRAIN, NCTRL, NHVG, KD, CHUNK = 1000, 40, 200, 80, 6000, 2000, -0.25, 20000
RATIO = float(os.environ.get("RATIO", "0.7"))
CTRL = "non-targeting"


def load(fn, t0):
    with h5py.File(fn, "r") as f:
        gc, gcats = read_codes_cats(f, "gene"); ec, ecats = read_codes_cats(f, "sgID_AB")
        umi = f["obs"]["UMI_count"][:].astype(np.float32); genes = [str(g) for g in var_gene_names(f)]
        M = f["X"]; N = M.shape[0]; parts = []
        for s in range(0, N, CHUNK):
            parts.append(csr_matrix(M[s:min(N, s + CHUNK)].astype(np.float32)))
        X = vstack(parts).tocsr()
    print(fn, X.shape, "%.0fs" % (time.time() - t0), flush=True)
    return dict(target=np.asarray(gcats)[gc], elem=np.asarray(ecats)[ec], umi=umi, genes=genes, X=X)


def main():
    kf, rf, out = sys.argv[1:4]; t0 = time.time(); os.makedirs(out, exist_ok=True)
    D = {st: load(fn, t0) for st, fn in zip(STATES, (kf, rf))}
    common = sorted(set(D["K562"]["genes"]) & set(D["RPE1"]["genes"]))
    for st in STATES:
        d = D[st]; gp_ = {g: j for j, g in enumerate(d["genes"])}
        d["X"] = d["X"][:, [gp_[g] for g in common]].tocsr()
        L = d["X"].copy(); L.data = np.log1p(L.data / np.repeat(d["umi"], np.diff(L.indptr)) * 1e4); d["L"] = L
        d["ctrl_mean"] = np.asarray(L[d["target"] == CTRL].mean(0)).ravel()
        e, n = np.unique(d["elem"], return_counts=True); d["ne"] = dict(zip(e, n))
    gpos = {g: j for j, g in enumerate(common)}
    # elements per target with >= MINC cells in both screens
    t_of = {}
    for st in STATES:
        for t, e in zip(D[st]["target"], D[st]["elem"]):
            t_of[e] = t
    ok = {}
    for e, t in t_of.items():
        if t == CTRL:
            continue
        n = min(D[st]["ne"].get(e, 0) for st in STATES)
        if n >= MINC:
            ok.setdefault(t, []).append((n, e))
    pairs = {t: sorted(e for _, e in sorted(v, reverse=True)[:2]) for t, v in ok.items() if len(v) >= 2}

    def kd(st, e, j):
        d = D[st]; m = d["elem"] == e
        return float(np.asarray(d["L"][m][:, j].mean())) - d["ctrl_mean"][j]
    inert, strong, other, kds = [], [], [], {}
    for t, (e0, e1) in sorted(pairs.items()):
        j = gpos.get(t)
        if j is None or all(D[st]["ctrl_mean"][j] < 0.01 for st in STATES):
            inert.append(t); continue
        good = True
        kds[t] = {st: [float(kd(st, e0, j)), float(kd(st, e1, j)), float(D[st]["ctrl_mean"][j])] for st in STATES}
        for st in STATES:
            if D[st]["ctrl_mean"][j] < 0.5:
                good = False; break
            k0, k1 = kds[t][st][:2]; s, w = min(k0, k1), max(k0, k1)
            good &= w <= KD and w / s >= RATIO
        (strong if good else other).append(t)
    test = set(inert) | set(strong) | set(other)
    rng = np.random.default_rng(5)
    pool = sorted(set(ok) - set(pairs))
    train = sorted(rng.choice(pool, min(NTRAIN, len(pool)), replace=False).tolist())
    elem_of = {t: [max(ok[t])[1]] for t in train}
    elem_of.update({t: list(pairs[t]) for t in test})
    # genes: expressed in both, ranked by variance of training-target means
    score, expr = np.zeros(len(common)), np.ones(len(common), bool)
    for st in STATES:
        d = D[st]; expr &= d["ctrl_mean"] >= 0.05
        M = np.stack([np.asarray(d["L"][d["elem"] == elem_of[t][0]].mean(0)).ravel() for t in train])
        v = M.var(0); v = (v - v[expr].mean()) / v[expr].std(); score += v
    score[~expr] = -np.inf
    gsel = sorted(set(np.argsort(-score)[:NHVG].tolist()) | {gpos[t] for t in test if t in gpos})
    json.dump(dict(train=train, inert=inert, strong=strong, mismatched=other, extra_kinds=["mismatched"], dropped_test=[],
                   n_genes=len(gsel), states=STATES,
                   elements={t: elem_of[t] for t in sorted(test)}, n_pairs=len(pairs), kd=kds, ratio=RATIO,
                   caps=dict(test=CAP_TEST, train=CAP_TRAIN, ctrl=NCTRL, minc=MINC), rules=__doc__),
              open(out + "/pilot_design.json", "w"), indent=1)
    print("targets with 2 elements", len(pairs), "inert", len(inert), "strong", len(strong), "mismatched", len(other),
          "train", len(train), "genes", len(gsel), "%.0fs" % (time.time() - t0), flush=True)
    rng = np.random.default_rng(6)
    for st in STATES:
        d = D[st]; sel, kind, gp = [], [], []
        for names, k, cap in ((train, "train", CAP_TRAIN), (inert, "inert", CAP_TEST), (strong, "strong", CAP_TEST),
                              (other, "mismatched", CAP_TEST)):
            for t in names:
                for p, e in enumerate(elem_of[t]):
                    r = np.where((d["elem"] == e) & (d["umi"] > 0))[0]
                    pick = rng.choice(r, min(cap, len(r)), replace=False); sel += pick.tolist(); kind += [k] * len(pick); gp += [p] * len(pick)
        r = np.where((d["target"] == CTRL) & (d["umi"] > 0))[0]; pick = rng.choice(r, min(NCTRL, len(r)), replace=False)
        sel += pick.tolist(); kind += ["control"] * len(pick); gp += [0] * len(pick)
        sel = np.array(sel); o = np.argsort(sel, kind="stable"); sel = sel[o]
        X = d["X"][sel][:, gsel].tocsr()
        np.savez(out + "/pilot_%s.npz" % st, data=X.data.astype(np.float32), indices=X.indices, indptr=X.indptr,
                 shape=np.array(X.shape), umi=d["umi"][sel], target=d["target"][sel], gp=np.array(gp)[o],
                 kind=np.array(kind)[o], genes=np.array([common[j] for j in gsel]), row=sel)
        print(st, "cells", X.shape[0], "genes", X.shape[1], "%.0fs" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
