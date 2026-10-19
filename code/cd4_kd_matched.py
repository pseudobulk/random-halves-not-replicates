"""Does the readout vanish under the guide split only because the two guides knock down differently? (preprint review
2.7; A100; CD4 donor-1 single cells, Rest vs Stim48hr; aggregates of cd4_sc_extract.py / cd4_sc_extract_guide.py.)
Expressed targets (control >= 0.5 in both states, target measured) present in both designs. Knockdown of each guide
half = target gene mean in that half - control mean of the same half (log1p CP10K), per state. Subsets:
  expressed      all such targets
  kd_matched     both guide halves knock down (<= KD in both states) and the weaker guide reaches >= RATIO of the
                 stronger one in both states, i.e. the two halves are the same perturbation in strength
For each subset and design (random cell halves vs guide halves): gamma share, readout share (k = 10, 20; state = top-k
PCs of half-A beta, ridge with lambda from inner CV, 2 folds; template = subset mean) and shuffled control.
If the readout stays ~0 under the guide split in kd_matched, guide efficiency mismatch does not explain its loss.
Usage: cd4_kd_matched.py <out_Rest> <out_Stim> <out_Rest_guide> <out_Stim_guide> <out_dir>. Exploratory."""
import os, sys, json
import numpy as np
from cd4_readout import load, MIN_HALF
from robust_cd4_tahoe import predict, expl

KD, RATIO = -0.25, 0.7


def main():
    dr, ds, drg, dsg, out = sys.argv[1:6]; os.makedirs(out, exist_ok=True)
    D = {"cell": (load(dr), load(ds)), "guide": (load(drg), load(dsg))}
    gst = {"rest": list(D["cell"][0][1]), "stim": list(D["cell"][1][1])}
    genes = sorted(set(gst["rest"]) & set(gst["stim"]))
    col = {st: np.array([{g: j for j, g in enumerate(gst[st])}[g] for g in genes]) for st in gst}
    gall = {g: j for j, g in enumerate(genes)}
    ctrl = {}
    for nm, (r, s_) in D.items():
        for st, (groups, _, S, n) in (("rest", r), ("stim", s_)):
            ci = groups.index("__ctrl__"); ctrl[(nm, st)] = np.stack([np.asarray(S[ci, h])[col[st]] / n[ci, h] for h in range(2)])
    ok = None
    for nm, (r, s_) in D.items():
        for groups, _, S, n in (r, s_):
            st_ = set(t for i, t in enumerate(groups) if t != "__ctrl__" and n[i].min() >= MIN_HALF)
            ok = st_ if ok is None else ok & st_
    cm = {st: ctrl[("cell", st)].mean(0) for st in ("rest", "stim")}
    targets = sorted(t for t in ok if t in gall and min(cm["rest"][gall[t]], cm["stim"][gall[t]]) >= 0.5)
    gk = np.where(np.minimum(cm["rest"], cm["stim"]) >= 0.1)[0]
    H, kd = {}, {}
    for nm, (r, s_) in D.items():
        arr = []
        for st, (groups, _, S, n) in (("rest", r), ("stim", s_)):
            gi = {g: i for i, g in enumerate(groups)}; idx = np.array([gi[t] for t in targets])
            X = np.stack([np.asarray(S[idx, h])[:, col[st]] / n[idx, h][:, None] - ctrl[(nm, st)][h] for h in range(2)])
            tcol = np.array([gall[t] for t in targets])
            kd[(nm, st)] = X[:, np.arange(len(targets)), tcol].T          # targets x 2 halves
            for k_, t in enumerate(targets):
                X[:, k_, gall[t]] = 0
            arr.append(X[..., gk])
        H[nm] = arr
    m = np.ones(len(targets), bool)
    for st in ("rest", "stim"):
        k2 = kd[("guide", st)]
        strong, weak = np.minimum(k2[:, 0], k2[:, 1]), np.maximum(k2[:, 0], k2[:, 1])     # more negative = stronger
        m &= (weak <= KD) & (weak / strong >= RATIO)
    res = {"rules": __doc__, "n_expressed": len(targets), "n_kd_matched": int(m.sum()), "KD": KD, "RATIO": RATIO,
           "kd_matched_median_kd": {st: [float(np.median(kd[("guide", st)][m, h])) for h in range(2)] for st in ("rest", "stim")}}
    for sname, mm in (("expressed", np.ones(len(targets), bool)), ("kd_matched", m)):
        for nm in H:
            R_, S_ = H[nm]
            c0 = lambda M: M - M.mean(0)
            XA, XB, YA, YB = (c0(M[mm]) for M in (R_[0], R_[1], S_[0], S_[1]))
            bA, gA, gB = (XA + YA) / 2, (XA - YA) / 2, (XB - YB) / 2
            den = float((gA * gB).sum())
            d = {"n": int(mm.sum()), "gamma_share": den / (den + float((bA * (XB + YB) / 2).sum())),
                 "gamma_cross_per_target": den / int(mm.sum())}
            for k in (10, 20):
                d["readout_k%d" % k] = expl(predict(bA, gB, k), gB) / den
            d["readout_k20_shuffled"] = expl(predict(bA, gB, 20, shuffle=True), gB) / den
            res["%s|%s" % (sname, nm)] = d
            print(sname, nm, json.dumps({k: round(v, 4) for k, v in d.items()}), flush=True)
    json.dump(res, open(out + "/metrics_cd4_kd_matched.json", "w"), indent=1)
    print("saved", flush=True)


if __name__ == "__main__":
    main()
