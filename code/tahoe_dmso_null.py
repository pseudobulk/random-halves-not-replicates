"""Tahoe null perturbation: the difference between a plate's two DMSO wells (preprint review 2.4; A100).
Both wells received vehicle only, so in expectation their difference is zero; any part of it that is reproducible
between random cell halves of the same wells is a well effect (spheroid, density, handling), the Tahoe analogue of
the inert-target control. For each plate with >= 2 DMSO wells (first two by sample order) and each line with >= MIN
cells in both halves of both wells:
  x_h = mean(well 1, half h) - mean(well 2, half h)        (log1p CP10K pseudobulk, genes with DMSO mean >= 0.1)
  within-well cross energy x_0 . x_1 (= about twice the per-well effect variance), Pearson(x_0, x_1), frac_pos
  split across lines into the well-common part (mean over lines, entering beta) and the line-specific part (gamma)
  cross-plate check: correlation of (well 1 - well 2) between plates for the same line (a fixed well position would
  give a consistent sign; independent well states give ~0)
  scale: within-well cross energy of drug effects on the same plate (effect vs the plate's DMSO wells, per half, 200
  random drug samples), so the well terms can be compared with what a drug condition shows.
Usage: tahoe_dmso_null.py <tahoe_agg dir> <out_dir>. Exploratory."""
import os, sys, json, ast
import numpy as np
import pyarrow.parquet as pq
import config

MIN = 20


def main():
    agg, out = sys.argv[1:3]; os.makedirs(out, exist_ok=True)
    idx = json.load(open(agg + "/index.json")); lines, samples, genes = idx["lines"], idx["samples"], idx["genes"]
    C, Ps = len(lines), len(samples)
    S = np.load(agg + "/sums64.npy", mmap_mode="r"); n = np.load(agg + "/counts.npy")[:, 0]
    sm = pq.read_table(config.tahoe_meta("sample_metadata.parquet")).to_pydict()
    info = {s: (d.strip(), p) for s, d, p in zip(sm["sample"], sm["drug"], sm["plate"])}
    row = lambda c, s, h: (c * Ps + s) * 2 + h
    plates = sorted({info[s][1] for s in samples})
    dmso = {p: [i for i, s in enumerate(samples) if info[s][1] == p and info[s][0].startswith("DMSO")] for p in plates}
    drugs = {p: [i for i, s in enumerate(samples) if info[s][1] == p and not info[s][0].startswith("DMSO")] for p in plates}
    rr = [row(c, s, h) for p in plates for s in dmso[p][:2] for c in range(C) for h in range(2)]
    gm = np.asarray(S[rr].sum(0)) / max(n[rr].sum(), 1); gk = np.where(gm >= 0.1)[0]
    mean = lambda r: np.asarray(S[r][gk]) / n[r]
    res = {"rules": __doc__, "n_genes": int(len(gk)), "plates": {}}
    rng = np.random.default_rng(5)
    diffs = {}
    for p in plates:
        if len(dmso[p]) < 2:
            continue
        w1, w2 = dmso[p][:2]
        okc = [c for c in range(C) if min(n[row(c, w, h)] for w in (w1, w2) for h in range(2)) >= MIN]
        if len(okc) < 5:
            continue
        X = np.array([[mean(row(c, w1, h)) - mean(row(c, w2, h)) for h in range(2)] for c in okc])   # lines x 2 x genes
        diffs[p] = (okc, X.mean(1))
        e = (X[:, 0] * X[:, 1]).sum(1)
        a, b = X[:, 0] - X[:, 0].mean(1, keepdims=True), X[:, 1] - X[:, 1].mean(1, keepdims=True)
        r = (a * b).sum(1) / np.sqrt((a ** 2).sum(1) * (b ** 2).sum(1))
        com = X.mean(0); spec = X - com
        d = dict(n_lines=len(okc), cross_median=float(np.median(e)), frac_pos=float((e > 0).mean()), pearson_median=float(np.median(r)),
                 common_cross=float((com[0] * com[1]).sum() * len(okc)), specific_cross=float((spec[:, 0] * spec[:, 1]).sum()))
        # drug scale: effect vs the plate's two DMSO wells (same half), within-well cross energy
        ref = {c: [np.mean([mean(row(c, w, h)) for w in (w1, w2)], 0) for h in range(2)] for c in okc}
        de = []
        for s in rng.choice(drugs[p], min(200, len(drugs[p])), replace=False):
            for c in okc:
                if min(n[row(c, s, h)] for h in range(2)) >= MIN:
                    x0, x1 = mean(row(c, s, 0)) - ref[c][0], mean(row(c, s, 1)) - ref[c][1]
                    de.append(float((x0 * x1).sum()))
        d["drug_cross_median"] = float(np.median(de)); d["drug_cross_q25"] = float(np.percentile(de, 25)); d["n_drug_pairs"] = len(de)
        d["well_per_well_over_drug_median"] = float(np.median(e) / 2 / np.median(de))
        res["plates"][str(p)] = d
        print(p, json.dumps(d), flush=True)
    # cross-plate: same line, (well 1 - well 2) correlation between plates
    pl = sorted(diffs); xs = []
    for i in range(len(pl)):
        for j in range(i + 1, len(pl)):
            (ci, Xi), (cj, Xj) = diffs[pl[i]], diffs[pl[j]]
            for c in set(ci) & set(cj):
                a, b = Xi[ci.index(c)], Xj[cj.index(c)]
                a, b = a - a.mean(), b - b.mean()
                xs.append(float((a * b).sum() / np.sqrt((a ** 2).sum() * (b ** 2).sum())))
    res["cross_plate_corr"] = dict(n=len(xs), mean=float(np.mean(xs)), median=float(np.median(xs)), sd=float(np.std(xs)))
    print("cross_plate", json.dumps(res["cross_plate_corr"]), flush=True)
    v = list(res["plates"].values())
    res["summary"] = {k: float(np.median([x[k] for x in v])) for k in ("cross_median", "frac_pos", "pearson_median", "drug_cross_median", "well_per_well_over_drug_median")}
    res["summary"]["common_share_of_cross"] = float(np.median([x["common_cross"] / (x["common_cross"] + x["specific_cross"]) for x in v]))
    print("summary", json.dumps(res["summary"]), flush=True)
    json.dump(res, open(out + "/metrics_tahoe_dmso_null.json", "w"), indent=1)
    print("saved", flush=True)


if __name__ == "__main__":
    main()
