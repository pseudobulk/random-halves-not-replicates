"""Whole-brain: the independent replicate of this design is another set of mice. Copy wb_meta.npz with the half bit
replaced by a mouse group (mice split at random into two groups of equal size, seed 12), so that wb_aggregate.py gives
mouse halves: the nuclei of one guide in group-0 mice against its nuclei in group-1 mice.
Usage: wb_mouse_halves.py <wb_meta.npz> <out.npz>"""
import sys
import numpy as np

z = np.load(sys.argv[1], allow_pickle=True)
z = {k: (z[k].astype(str) if z[k].dtype == object else z[k]) for k in z.files}
mice = np.unique(z["mouse"]); rng = np.random.default_rng(12)
grp = dict(zip(mice, (rng.permutation(len(mice)) < len(mice) // 2).astype(np.int8)))
z["half"] = np.array([grp[m] for m in z["mouse"]], np.int8)
print("mice", len(mice), "group sizes", np.bincount(list(grp.values())).tolist(), "nuclei per half", np.bincount(z["half"]).tolist())
np.savez(sys.argv[2], **z)
