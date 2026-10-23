"""Model evaluation step 3: clonality in expression, on one scale for CD4 and Frangieh (the expression analogue of the
TRBC vc of me_clonality.py; TRBC exists only in T cells). Null units: guides with no true effect (CD4: guides of
unexpressed targets, kind 'inert'; Frangieh: non-targeting sgRNAs, kind 'null').
For a unit, its cells are split into random halves A and B, the control cells into halves C1 and C2 (seed 11);
d_A = mean_A - mean_C1, d_B = mean_B - mean_C2 (log1p CP10K, the run's evaluation genes, target gene excluded).
  vc_expr   mean over units of sum_g d_A,g d_B,g / sum_g var_g, with var_g the per-cell variance of control cells.
            Sampling noise of different cells averages out of the cross product, so this estimates the squared
            deviation that all of the unit's cells share (its clone mix) relative to the per-cell variance, at any
            cell number; it is what the shared-clone truth A shares with the input.
  vc_guide  the same cross product between the two guides of a null target (independent transductions; ~0 expected)
  *_ds      the same with every unit first subsampled to NDS cells (mean of NREP draws), as a check on the
            cell-number independence
SEs over targets (a target's two guides averaged). Usage (env PILOT_DONOR or PILOT_STATES): me_vc_expr.py <data_dir>
<kind> <out_json>"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "4")
import sys, json
import numpy as np
from pilot_common import STATES, DONOR, load, lognorm, design

NDS, NREP, MINU = 40, 5, 20


def cross(L, ra, rb, c1, c2, keep):
    da = np.asarray(L[ra].mean(0)).ravel() - c1; db = np.asarray(L[rb].mean(0)).ravel() - c2
    return float((da * db)[keep].sum())


def main():
    d, kind, outp = sys.argv[1:4]
    des = design(d); res = {"donor": DONOR, "kind": kind, "rules": __doc__}
    for st in STATES:
        D = load(d, st); L = lognorm(D["X"], D["umi"]).tocsr(); genes = D["genes"]; gpos = {g: j for j, g in enumerate(genes)}
        rng = np.random.default_rng(11)
        ci = np.where(D["kind"] == "control")[0]; rng.shuffle(ci); C1, C2 = ci[: len(ci) // 2], ci[len(ci) // 2:]
        c1, c2 = np.asarray(L[C1].mean(0)).ravel(), np.asarray(L[C2].mean(0)).ravel()
        Lc = L[ci]; var = np.asarray(Lc.multiply(Lc).mean(0)).ravel() - np.asarray(Lc.mean(0)).ravel() ** 2
        units = {}
        for i in np.where(D["kind"] == kind)[0]:
            units.setdefault(D["lab"][i], []).append(i)
        per_t = {}
        for u, r in sorted(units.items()):
            r = np.array(r); t = u.split("|")[0]
            if len(r) < MINU:
                continue
            keep = np.ones(len(genes), bool)
            if t in gpos:
                keep[gpos[t]] = False
            den = var[keep].sum()
            p = rng.permutation(r); full = cross(L, p[: len(p) // 2], p[len(p) // 2:], c1, c2, keep) / den
            ds = []
            for _ in range(NREP):
                q = rng.choice(r, min(NDS, len(r)), replace=False)
                ds.append(cross(L, q[: len(q) // 2], q[len(q) // 2:], c1, c2, keep) / den)
            per_t.setdefault(t, {})[u] = (full, float(np.mean(ds)), r)
        rows = []
        for t, us in per_t.items():
            full = np.mean([v[0] for v in us.values()]); ds = np.mean([v[1] for v in us.values()])
            g = None
            if len(us) == 2:
                (ra, rb) = [v[2] for v in us.values()]
                keep = np.ones(len(genes), bool)
                if t in gpos:
                    keep[gpos[t]] = False
                g = cross(L, ra, rb, c1, c2, keep) / var[keep].sum()
            rows.append((full, ds, g))
        A = np.array([r[0] for r in rows]); B = np.array([r[1] for r in rows]); G = np.array([r[2] for r in rows if r[2] is not None])
        se = lambda x: float(x.std(ddof=1) / np.sqrt(len(x)))
        res[st] = dict(n_units=int(sum(len(v) for v in per_t.values())), n_targets=len(rows), vc_expr=float(A.mean()), vc_expr_se=se(A),
                       vc_expr_ds=float(B.mean()), vc_expr_ds_se=se(B), vc_guide=float(G.mean()) if len(G) else None,
                       vc_guide_se=se(G) if len(G) > 1 else None, n_genes=len(genes),
                       median_cells_per_unit=float(np.median([len(v[2]) for us in per_t.values() for v in us.values()])))
        print(DONOR, st, json.dumps(res[st]), flush=True)
    json.dump(res, open(outp, "w"), indent=1)


if __name__ == "__main__":
    main()
