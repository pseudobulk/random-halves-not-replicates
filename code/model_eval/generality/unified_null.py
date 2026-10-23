"""Unified clone ruler for h5ad screens (null_ruler.py, both references). Reads only what the ruler needs: all
non-targeting cells and the cells of up to NTG randomly chosen targeting units, straight from the h5ad with h5py (dense
or CSR X; old-style obs with __categories or new-style categorical obs). Counts are log1p CP10K; if X is not integer
counts (max fractional part > 0) it is used as given.
Columns are given per dataset: unit (guide or guide pair), target (gene; non-targeting cells are those whose target
contains 'non-targeting' or equals 'control' / 'non_target' / 'NTC'), sample (10x lane or batch; reference group and
replicate label), optional quality column (cells with a true value are dropped).
Options: qc=<col> drops cells whose value is true; keep=<col>:<value> keeps only cells with that value (repeatable);
drop=<col>:<value> drops cells with that value (repeatable); ntg=<n> targeting units sampled (default 400).
Usage: unified_null.py <h5ad> <unit col> <target col | @guide_type> <sample col or -> <out.json> [options]"""
import sys, json
import numpy as np, h5py
from scipy.sparse import csr_matrix, vstack
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from null_ruler import ruler

NULL = ("non-targeting", "non_targeting", "nontargeting", "control", "non_target", "ntc", "safe-targeting")


def col(f, k):
    o = f["obs"]
    if "__categories" in o and k in o["__categories"]:
        cats = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in o["__categories"][k][:]])
        codes = o[k][:]; return np.where(codes >= 0, cats[np.maximum(codes, 0)], "")
    g = o[k]
    if isinstance(g, h5py.Group):
        cats = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]])
        codes = g["codes"][:]; return np.where(codes >= 0, cats[np.maximum(codes, 0)], "")
    v = g[:]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v]) if v.dtype.kind in "OS" else v.astype(str)


def rows(f, idx):
    X = f["X"]; idx = np.sort(idx); parts = []
    if isinstance(X, h5py.Group) and str(X.attrs.get("encoding-type", "")).startswith("csc"):
        from scipy.sparse import csc_matrix                       # column-compressed: read whole, then select rows
        M = csc_matrix((X["data"][:].astype(np.float32), X["indices"][:], X["indptr"][:]), shape=tuple(X.attrs["shape"]))
        return M.tocsr()[idx], idx
    if isinstance(X, h5py.Group):
        n = X.attrs["shape"][1]; N = X.attrs["shape"][0]; B = 20000
        for r0 in range(0, N, B):
            r1 = min(r0 + B, N); w = idx[(idx >= r0) & (idx < r1)]
            if not len(w):
                continue
            ip = X["indptr"][r0:r1 + 1]; a, b = int(ip[0]), int(ip[-1])
            M = csr_matrix((X["data"][a:b].astype(np.float32), X["indices"][a:b], ip - a), shape=(r1 - r0, n))
            parts.append(M[w - r0])
        return vstack(parts).tocsr(), idx
    for s in range(0, len(idx), 2000):
        parts.append(csr_matrix(X[idx[s:s + 2000]].astype(np.float32)))
    return vstack(parts).tocsr(), idx


def main():
    fn, ucol, tcol, scol, outp = sys.argv[1:6]
    opts = sys.argv[6:]; qcol = None; NTG = 400; keeps, drops = [], []
    for o in opts:
        k, v = o.split("=", 1)
        if k == "qc":
            qcol = v
        elif k == "ntg":
            NTG = int(v)
        elif k == "keep":
            keeps.append(v.split(":", 1))
        elif k == "drop":
            drops.append(v.split(":", 1))
    rng = np.random.default_rng(0)
    with h5py.File(fn, "r") as f:
        unit = col(f, ucol)
        if tcol == "@guide_type":            # CD4 (Zhu): null from guide_type, target from the guide name "GENE-1"
            gt = col(f, "guide_type")
            target = np.where(gt == "non-targeting", "non-targeting", np.array([x.rsplit("-", 1)[0] for x in unit]))
        else:
            target = col(f, tcol)
        sample = col(f, scol) if scol != "-" else np.array(["all"] * len(unit))
        keep = np.array([u != "" and t != "" for u, t in zip(unit, target)])
        if qcol:
            q = col(f, qcol); keep &= ~np.isin(np.char.lower(q.astype(str)), ["true", "1"])
        for c, v in keeps:
            keep &= col(f, c) == v
        for c, v in drops:
            keep &= col(f, c) != v
        isn = np.array([any(t.lower() == x or t.lower().startswith(x) for x in NULL) for t in target]) & keep
        tg_units = np.unique(unit[keep & ~isn])
        pick = set(rng.choice(tg_units, min(NTG, len(tg_units)), replace=False))
        sel = np.where(isn | (keep & np.isin(unit, list(pick))))[0]
        X, sel = rows(f, sel)
    frac = X.data[:10000] - np.floor(X.data[:10000])
    if frac.max() == 0:
        tot = np.asarray(X.sum(1)).ravel()
        X = X.multiply(1e4 / np.maximum(tot, 1)[:, None]).tocsr().astype(np.float32); X.data = np.log1p(X.data)
    u, t, s, n = unit[sel], target[sel], sample[sel], isn[sel]
    info = dict(file=fn, options=opts, cells_read=int(len(sel)), null_cells=int(n.sum()), null_units=int(len(np.unique(u[n]))),
                target_units_sampled=len(pick), samples=int(len(np.unique(s))))
    print(info, flush=True)
    out = dict(rules=__doc__, info=info)
    for mode in ("all", "null"):
        out["ruler_" + mode] = ruler(X, s, u.astype(object), t, n, s, ref_mode=mode)
        print(mode, json.dumps(out["ruler_" + mode]), flush=True)
    json.dump(out, open(outp, "w"), indent=1)


if __name__ == "__main__":
    main()
