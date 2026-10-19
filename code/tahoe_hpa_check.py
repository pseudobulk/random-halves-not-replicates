"""Spheroid confound check (gap_analysis §12; A100). In Tahoe, a few programs' gains (EMT, inflammatory, hypoxia,
glycolysis) rise with the line's basal expression of that program, beyond coherence-matched nulls (Z1 ab). Tahoe cells
grow as pooled 3D spheroids, so the basal hypoxia / glycolysis level may reflect the spheroid micro-environment rather
than the line. Test: use each line's basal expression in ordinary monolayer culture (HPA rna_celline nTPM) instead.
If the HPA basal level predicts the Tahoe gain about as well, the property is cell-intrinsic; if only the Tahoe basal
level does, the micro-environment is a likely cause.
Per program (Hallmark + curated as tahoe_arms): Spearman between Tahoe-basal and HPA-basal program scores across
matched lines; LOO R2 of Tahoe gain from Tahoe-basal and from HPA-basal; HPA null = expression-decile-matched random
gene sets (their HPA basal mean -> their own Tahoe gain).
Usage: tahoe_hpa_check.py <tahoe_program dir> <tahoe_arms dir> <hallmark txt> <hpa rna_celline.tsv> <out_dir>"""
import os, sys, json, re, csv
import numpy as np
import pyarrow.parquet as pq
import config
from tahoe_arms import kr_loo, HALL, ISR, HSF1, NRF2, P53
from tahoe_tone_check import gains

NRAND = 30
norm = lambda s: re.sub(r"[^A-Z0-9]", "", str(s).upper())


def main():
    src, arms, hpath, hpa, out = sys.argv[1:6]; os.makedirs(out, exist_ok=True)
    z = np.load(src + "/tahoe_arrays.npz", allow_pickle=True)
    lines = [str(x) for x in z["lines"]]; genes = [str(g) for g in z["genes"]]; gpos = {g: j for j, g in enumerate(genes)}
    resp = z["resp"]; XA = z["XA"][:, resp].astype(np.float64); XB = z["XB"][:, resp].astype(np.float64); base = z["base"].astype(np.float64)
    L = XA.shape[0]; sA, sB = XA.sum(0), XB.sum(0); refA = (sA[None] - XA) / (L - 1); refB = (sB[None] - XB) / (L - 1)
    tA, tB = XA - refA, XB - refB
    t = pq.read_table(config.tahoe_meta("cell_line_metadata.parquet")); cols = t.column_names; d = t.to_pydict()
    print("metadata columns", cols, flush=True)
    namecols = [c for c in cols if re.search(r"name|cell_line|Cell_Name|cell", c, re.I) and c != "Cell_ID_Cellosaur"]
    names = {}
    for i, cid in enumerate(d["Cell_ID_Cellosaur"]):
        for c in namecols:
            v = d[c][i]
            if isinstance(v, str) and v:
                names.setdefault(cid, set()).add(norm(v))
    want = {nm: cid for cid, ns in names.items() for nm in ns}
    H2 = {}                                                  # HPA: line -> gene -> log1p nTPM
    with open(hpa) as f:
        r = csv.reader(f, delimiter="\t"); hdr = next(r)
        gi, li, vi = hdr.index("Gene name"), hdr.index("Cell line"), hdr.index("nTPM")
        for row in r:
            nm = norm(row[li])
            if nm in want and row[gi] in gpos:
                H2.setdefault(want[nm], {})[row[gi]] = np.log1p(float(row[vi]))
    matched = [i for i, c in enumerate(lines) if c in H2]
    print("matched lines", len(matched), "of", L, [lines[i] for i in matched][:10], flush=True)
    hb = np.full((len(matched), len(genes)), np.nan)
    for k, i in enumerate(matched):
        for g, v in H2[lines[i]].items():
            hb[k, gpos[g]] = v
    okg = ~np.isnan(hb).any(0)
    H = {}
    for line in open(hpath):
        parts = [x for x in line.rstrip("\n").split("\t") if x]
        if len(parts) > 1:
            H[parts[0]] = set(parts[1:])
    S = {h: H[h] for h in HALL if h in H}; S.update(ISR_ATF4=set(ISR), HSF1_heat_shock=set(HSF1), NRF2_oxidative=set(NRF2), p53_core=set(P53))
    M = np.array(matched); tb = base[M] - base[M].mean(0); hbr = hb - np.nanmean(hb, 0)
    sd = lambda F: F / np.maximum(F.std(0), 1e-12)
    mu = base.mean(0); dec = np.digitize(mu, np.quantile(mu, np.linspace(0.1, 0.9, 9))); rng = np.random.default_rng(31)
    res = {"rules": __doc__, "n_matched": len(matched), "per_program": {}}
    for h, s in S.items():
        ix = np.array(sorted(gpos[g] for g in s if g in gpos and okg[gpos[g]]))
        if len(ix) < 10:
            continue
        g_ = gains(tA, tB, refA, refB, ix)[M]
        ts, hs = tb[:, ix].mean(1), hbr[:, ix].mean(1)
        from scipy.stats import spearmanr
        rr = []
        for _ in range(NRAND):
            rix = np.array([rng.choice(np.where((dec == dec[j]) & okg)[0]) for j in ix])
            rr.append(kr_loo(sd(hbr[:, rix].mean(1, keepdims=True)), gains(tA, tB, refA, refB, rix)[M][:, None]))
        e = dict(n_genes=int(len(ix)), rho_tahoe_vs_hpa_basal=float(spearmanr(ts, hs)[0]),
                 r2_tahoe_basal=kr_loo(sd(ts[:, None]), g_[:, None]), r2_hpa_basal=kr_loo(sd(hs[:, None]), g_[:, None]),
                 hpa_null_mean=float(np.mean(rr)), hpa_null_q95=float(np.quantile(rr, 0.95)))
        res["per_program"][h] = e
        print(h, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in e.items()}, flush=True)
    json.dump(res, open(out + "/metrics_tahoe_hpa_check.json", "w"), indent=1)
    print("saved", flush=True)


if __name__ == "__main__":
    main()
