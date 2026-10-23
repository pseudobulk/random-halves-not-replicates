"""Model evaluation, no-shared-clone control inside the CD4 data (A100). Each donor's predictions (me_export_preds.py;
models trained on that donor) are scored against the stimulated cells of ANOTHER donor, which was transduced separately,
so no clone is shared with the input:
  A'  the same guide in donor e (same guide sequence, independent transduction)
  B'  the other guide of the same target in donor e
A' - B' is what a prediction carries that is specific to the guide but not to its clones (e.g. guide strength); the
within-donor A - B minus this is the part due to shared clones. Same score as pilot_common.score (Pearson over the
source donor's pilot genes, target gene excluded); the truth effect in donor e = guide mean minus control mean of the
48 h stimulated aggregates (cd4_sc_extract_cells.py; log1p CP10K, all of the guide's cells, so the truth is less noisy
than the capped within-donor truth). Test guides of the source donor whose two guides both have >= MINC cells in donor
e. SEs clustered by target. Usage: me_cross_donor.py <out_json> <donor>:<preds.npz>:<pilot_design.json>:<cells_stim_dir> ..."""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "4")
import sys, json
import numpy as np

MINC = 40
MODELS = ("mean", "copy", "additive", "linear", "state_transition")


def truth(cd, genes, targets):
    meta = json.load(open(cd + "/meta.json")); gi = {g: i for i, g in enumerate(meta["groups"])}; ci = len(meta["groups"]) - 1
    gpos = {g: j for j, g in enumerate(meta["genes"])}; cols = np.array([gpos[g] for g in genes])
    S = np.load(cd + "/sum8.npy", mmap_mode="r"); n8 = np.load(cd + "/n8.npy")
    ctrl = np.asarray(S[ci].sum(0), dtype=np.float64)[cols] / n8[ci].sum()
    eff = {}
    for t in targets:
        i = gi.get(t)
        if i is None or min(n8[i, :4].sum(), n8[i, 4:].sum()) < MINC:
            continue
        Si = np.asarray(S[i], dtype=np.float64)[:, cols]
        for p in (0, 1):
            eff["%s|g%d" % (t, p)] = Si[p * 4:(p + 1) * 4].sum(0) / n8[i, p * 4:(p + 1) * 4].sum() - ctrl
    return eff


def summarise(rows):
    out = {}
    for grp in ("inert", "strong", "all"):
        rr = [r for r in rows if grp == "all" or r[1] == grp]
        if not rr:
            continue
        tg = sorted({r[0] for r in rr}); o = dict(n_guides=len(rr), n_targets=len(tg))
        for k, f in (("A", lambda r: r[2]), ("B", lambda r: r[3]), ("AB", lambda r: r[2] - r[3])):
            per_t = np.array([np.mean([f(r) for r in rr if r[0] == t]) for t in tg])
            o[k] = float(np.mean([f(r) for r in rr])); o[k + "_se"] = float(per_t.std(ddof=1) / np.sqrt(len(tg)))
        out[grp] = o
    return out


def main():
    outp = sys.argv[1]; specs = [a.split(":") for a in sys.argv[2:]]
    res = {"rules": __doc__, "pairs": {}}
    for d, pf, dj, _ in specs:
        z = np.load(pf); genes = [str(g) for g in z["genes"]]; gpos = {g: j for j, g in enumerate(genes)}
        des = json.load(open(dj)); kinds = {t: k for k in ("inert", "strong") for t in des[k]}
        for e, _, _, cd in specs:
            if e == d:
                continue
            T = truth(cd, genes, sorted(kinds)); r = {}
            for m in MODELS:
                rows = []
                for t, k in kinds.items():
                    keep = np.ones(len(genes), bool)
                    if t in gpos:
                        keep[gpos[t]] = False
                    for p in (0, 1):
                        g, o = "%s|g%d" % (t, p), "%s|g%d" % (t, 1 - p)
                        key = m + "__" + g
                        if key not in z.files or g not in T or o not in T:
                            continue
                        pr = z[key][keep]
                        if pr.std() == 0:
                            continue
                        rows.append((t, k, np.corrcoef(pr, T[g][keep])[0, 1], np.corrcoef(pr, T[o][keep])[0, 1]))
                r[m] = summarise(rows)
            res["pairs"]["%s->%s" % (d, e)] = r
            v = r["state_transition"]["all"]; c = r["copy"]["all"]
            print(d, "->", e, "n", v["n_guides"], "STATE A' %.3f B' %.3f diff %.3f (%.3f) | copy diff %.3f (%.3f)" % (
                v["A"], v["B"], v["AB"], v["AB_se"], c["AB"], c["AB_se"]), flush=True)
    json.dump(res, open(outp, "w"), indent=1)


if __name__ == "__main__":
    main()
