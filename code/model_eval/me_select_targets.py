"""Model evaluation step 1 (model_eval/plan.md section 7): test targets for one CD4 donor by the cell-eval selection rule
of code/a100/cd4_celleval_subset.py (the rule that chose the donor-1 pilot targets), applied to the donor's resting
aggregates (cd4_sc_extract_cells.py). Only the target names are written; cd4_pilot_extract.py then keeps the targets
whose two guides have >= 40 cells in both states.
  candidates  both guide halves (gp = 0 / 1) >= MINC cells, target gene measured
  inert       control mean log1p CP10K < 0.01: NI at random (seed 4)
  strong      control mean >= 0.5, both guides knock the target down (Delta log1p <= KD) and the weaker guide reaches
              >= RATIO of the stronger: the NS strongest by mean knockdown
Output npz: inert_names, strong_names, strong_kd (the keys cd4_pilot_extract.py reads).
Usage: me_select_targets.py <cells_rest_dir> <out_npz>"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json
import numpy as np

NI, NS, MINC, KD, RATIO = 150, 150, 40, -0.25, 0.7   # as cd4_celleval_subset.py


def main():
    cd, outp = sys.argv[1:3]
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
    np.savez(outp, inert_names=np.array([groups[i] for i in pick_i]), strong_names=np.array([groups[i] for i in pick_s]),
             strong_kd=np.array([kd[i] for i in pick_s]))
    print("candidates", len(cand), "inert", len(inert), "strong matched", len(kd), "picked", len(pick_i), len(pick_s), flush=True)


if __name__ == "__main__":
    main()
