"""Sign flips and discrete-field equality between two npz files (per-target arrays). Usage: npz_signs.py <a.npz> <b.npz>"""
import sys, numpy as np
a, b = np.load(sys.argv[1]), np.load(sys.argv[2])
for k in a.files:
    x, y = a[k], b[k]
    if x.dtype.kind in "fc":
        flips = int(((x > 0) != (y > 0)).sum())
        print("%-28s sign flips %d  maxabs %.3g" % (k, flips, np.nanmax(np.abs(x - y)) if x.size else 0))
    else:
        print("%-28s equal %s" % (k, np.array_equal(x, y)))
