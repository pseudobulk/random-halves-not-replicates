"""Model evaluation step 3: a pilot-format data folder with the TEST guides' cells capped at K_source / K_target per
state (random, seed 9), to compare clone gains at the cells per guide of another data set (Frangieh test sgRNAs:
median 48 Control, 75 IFN-gamma cells). Training guides and controls are unchanged; the design is copied.
A - B depends on cells per guide (fewer cells, more sampling noise against the same clone deviation), so data sets are
compared at matched cells. Usage (env PILOT_DONOR or PILOT_STATES): me_downsample.py <data_dir> <out_dir> <K_src> <K_tgt>"""
import os, sys, json, shutil
import numpy as np
from scipy.sparse import csr_matrix
from pilot_common import STATES, design, test_targets


def main():
    d, out, ks, kt = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    os.makedirs(out, exist_ok=True); shutil.copy(os.path.join(d, "pilot_design.json"), out)
    test = test_targets(design(d)); rng = np.random.default_rng(9)
    for st, k in zip(STATES, (ks, kt)):
        z = dict(np.load(os.path.join(d, "pilot_%s.npz" % st), allow_pickle=False))
        X = csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
        keep = np.ones(X.shape[0], bool)
        lab = np.array(["%s|%d" % (t, g) for t, g in zip(z["target"], z["gp"])])
        for u in np.unique(lab[np.isin(z["target"], sorted(test)) & (z["kind"] != "control")]):
            r = np.where(lab == u)[0]
            if len(r) > k:
                keep[np.setdiff1d(r, rng.choice(r, k, replace=False))] = False
        X = X[keep]
        z.update(data=X.data, indices=X.indices, indptr=X.indptr, shape=np.array(X.shape))
        for f in z:
            if f not in ("data", "indices", "indptr", "shape", "genes") and len(z[f]) == len(keep):
                z[f] = z[f][keep]
        np.savez(os.path.join(out, "pilot_%s.npz" % st), **z)
        print(st, "cells", X.shape[0], "of", len(keep), "cap", k)


if __name__ == "__main__":
    main()
