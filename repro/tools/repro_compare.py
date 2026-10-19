"""Compare reproduced outputs with the originals.
  repro_compare.py pairs <pairs.tsv>      each line: label <TAB> original path <TAB> reproduced path (file or directory)
Files: .json (every numeric leaf; strings except 'rules'; key sets), .npz / .npy (every array). Directories: all such
files present in both, plus files present in only one. Prints one line per file: status, n values, max abs diff,
max rel diff, first mismatching key."""
import sys, os, json, math
import numpy as np

TOL = 1e-9


def leaves(x, p=""):
    if isinstance(x, dict):
        for k, v in x.items():
            yield from leaves(v, p + "/" + str(k))
    elif isinstance(x, list):
        for i, v in enumerate(x):
            yield from leaves(v, p + "[%d]" % i)
    else:
        yield p, x


ROUNDFLIPS = []


def roundflip(a, b):
    """lowest number of significant digits (1-3, the printed precision of the paper) at which a and b round differently"""
    for sig in (1, 2, 3):
        if float("%.*g" % (sig, a)) != float("%.*g" % (sig, b)):
            return "%d sig: %.*g vs %.*g" % (sig, sig, a, sig, b)
    return None


def cmp_vals(pairs):
    n, mabs, mrel, first, bad = 0, 0.0, 0.0, None, 0
    for k, a, b in pairs:
        if isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
            n += 1
            if not (isinstance(a, float) and math.isnan(a)) and a != b:
                rf = roundflip(a, b)
                if rf:
                    ROUNDFLIPS.append("%s: %s (%r vs %r)" % (k, rf, a, b))
            if (isinstance(a, float) and math.isnan(a)) and (isinstance(b, float) and math.isnan(b)):
                continue
            d = abs(a - b); r = d / max(abs(a), abs(b), 1e-300)
            mabs, mrel = max(mabs, d), max(mrel, r if d > 0 else 0)
            if d > TOL * max(1, abs(a)):
                bad += 1; first = first or "%s: %r vs %r" % (k, a, b)
        elif a != b:
            bad += 1; first = first or "%s: %r vs %r" % (k, str(a)[:60], str(b)[:60])
    return n, mabs, mrel, bad, first


def cmp_json(fa, fb):
    del ROUNDFLIPS[:]
    A, B = json.load(open(fa)), json.load(open(fb))
    la = {k: v for k, v in leaves(A) if not k.startswith("/rules")}
    lb = {k: v for k, v in leaves(B) if not k.startswith("/rules")}
    only = sorted(set(la) ^ set(lb))
    n, mabs, mrel, bad, first = cmp_vals([(k, la[k], lb[k]) for k in sorted(set(la) & set(lb))])
    if only:
        bad += len(only); first = first or "key only in one: %s" % only[0]
    return n, mabs, mrel, bad, first


def arrays(f):
    if f.endswith(".npz"):
        z = np.load(f, allow_pickle=True); return {k: z[k] for k in z.files}
    return {"": np.load(f, mmap_mode="r", allow_pickle=True)}


def cmp_arr(fa, fb):
    A, B = arrays(fa), arrays(fb)
    n, mabs, mrel, bad, first = 0, 0.0, 0.0, 0, None
    for k in sorted(set(A) | set(B)):
        if k not in A or k not in B:
            bad += 1; first = first or "array only in one: %s" % k; continue
        a, b = np.asarray(A[k]), np.asarray(B[k])
        if a.shape != b.shape:
            bad += 1; first = first or "%s shape %s vs %s" % (k, a.shape, b.shape); continue
        if a.dtype.kind in "fc" or b.dtype.kind in "fc":
            af, bf = a.astype(np.float64), b.astype(np.float64)
            both = np.isnan(af) & np.isnan(bf); d = np.where(both, 0, np.abs(af - bf)); d = np.nan_to_num(d, nan=np.inf)
            n += d.size
            if d.size:
                m = float(d.max()); mabs = max(mabs, m)
                sc = np.maximum(np.maximum(np.abs(af), np.abs(bf)), 1e-300)
                mrel = max(mrel, float(np.nan_to_num(np.where(d > 0, d / sc, 0), nan=np.inf).max()))
                nb = int((d > TOL * np.maximum(1, np.abs(np.nan_to_num(af)))).sum())
                if nb:
                    bad += nb; first = first or "%s: %d values differ, max %.3g" % (k, nb, m)
        else:
            n += a.size
            if not np.array_equal(a, b):
                nb = int((a != b).sum()) if a.shape else 1
                bad += nb; first = first or "%s: %d entries differ" % (k, nb)
    return n, mabs, mrel, bad, first


def cmp_file(fa, fb):
    if fa.endswith(".json"):
        return cmp_json(fa, fb)
    if fa.endswith((".npz", ".npy")):
        return cmp_arr(fa, fb)
    return None


def report(label, fa, fb):
    r = cmp_file(fa, fb)
    if r is None:
        return
    n, mabs, mrel, bad, first = r
    rf = len(ROUNDFLIPS) if fa.endswith(".json") else "-"
    print("%-6s %-50s n=%-9d maxabs=%-10.3g maxrel=%-10.3g roundflips=%-4s %s" % ("OK" if bad == 0 else "DIFF", label, n, mabs, mrel, rf, first or ""), flush=True)
    if fa.endswith(".json"):
        for r in ROUNDFLIPS:
            print("       ROUNDFLIP %s" % r, flush=True)


def main():
    for line in open(sys.argv[2]):
        if not line.strip() or line.startswith("#"):
            continue
        label, a, b = line.rstrip("\n").split("\t")
        if not os.path.exists(a) or not os.path.exists(b):
            print("MISS   %-50s %s" % (label, a if not os.path.exists(a) else b)); continue
        if os.path.isdir(a):
            fa = {f for f in os.listdir(a) if f.endswith((".json", ".npz", ".npy"))}
            fb = {f for f in os.listdir(b) if f.endswith((".json", ".npz", ".npy"))}
            for f in sorted(fa ^ fb):
                print("ONLY   %-50s %s" % (label + "/" + f, "original" if f in fa else "reproduced"))
            for f in sorted(fa & fb):
                report(label + "/" + f, os.path.join(a, f), os.path.join(b, f))
        else:
            report(label, a, b)


if __name__ == "__main__":
    main()
