"""Tone or scale? (Tahoe, A100.) tahoe_relay_check found that a program's own basal mean predicts the line's gain on that
program (LOO R2 up to 0.95), while upstream relay nodes add nothing. That could be biology (a program that runs higher
at rest is read more strongly when triggered) or scale (log1p effects grow with expression, so any gene set expressed
higher in a line shows larger effects there). Test: for each program, NRAND random gene sets of the same size,
matched gene by gene on mean basal expression decile, get their own gain vectors (same formula as tahoe_arms.py) and
the LOO R2 of their own basal mean -> own gain. Scale predicts random sets do as well as programs; tone predicts
programs exceed their matched random sets.
Usage: tahoe_tone_check.py <tahoe_program dir> <tahoe_arms dir> <hallmark txt> <out_dir>. Exploratory."""
import os, sys, json
import numpy as np
from tahoe_arms import kr_loo, HALL, ISR, HSF1, NRF2, P53

NRAND = 20


def gains(tA, tB, refA, refB, ix):
    num = ((tA[:, :, ix] * refB[:, :, ix]).sum((1, 2)) + (tB[:, :, ix] * refA[:, :, ix]).sum((1, 2))) / 2
    return num / np.maximum((refA[:, :, ix] * refB[:, :, ix]).sum((1, 2)), 1e-12)


def main():
    src, arms, hpath, out = sys.argv[1:5]; os.makedirs(out, exist_ok=True)
    z = np.load(src + "/tahoe_arrays.npz", allow_pickle=True)
    resp = z["resp"]; XA = z["XA"][:, resp].astype(np.float64); XB = z["XB"][:, resp].astype(np.float64)
    base = z["base"].astype(np.float64); genes = [str(g) for g in z["genes"]]; gpos = {g: j for j, g in enumerate(genes)}
    L = XA.shape[0]
    sA, sB = XA.sum(0), XB.sum(0); refA = (sA[None] - XA) / (L - 1); refB = (sB[None] - XB) / (L - 1)
    tA, tB = XA - refA, XB - refB
    programs = [str(p) for p in np.load(arms + "/program_gains.npz", allow_pickle=True)["programs"]]
    H = {}
    for line in open(hpath):
        parts = [x for x in line.rstrip("\n").split("\t") if x]
        if len(parts) > 1:
            H[parts[0]] = set(parts[1:])
    S = {h: H[h] for h in HALL if h in H}; S.update(ISR_ATF4=set(ISR), HSF1_heat_shock=set(HSF1), NRF2_oxidative=set(NRF2), p53_core=set(P53))
    relb = base - base.mean(0); mu = base.mean(0); dec = np.digitize(mu, np.quantile(mu, np.linspace(0.1, 0.9, 9)))
    rng = np.random.default_rng(13); res = {"rules": __doc__, "per_program": {}}
    sd = lambda F: F / np.maximum(F.std(0), 1e-12)
    allr = []
    for h in programs:
        ix = np.array(sorted(gpos[g] for g in S[h] if g in gpos))
        r_prog = kr_loo(sd(relb[:, ix].mean(1, keepdims=True)), gains(tA, tB, refA, refB, ix)[:, None])
        rr = []
        for _ in range(NRAND):
            rix = np.array([rng.choice(np.where(dec == dec[j])[0]) for j in ix])
            rr.append(kr_loo(sd(relb[:, rix].mean(1, keepdims=True)), gains(tA, tB, refA, refB, rix)[:, None]))
        allr += rr
        res["per_program"][h] = dict(n_genes=int(len(ix)), program=r_prog, random_mean=float(np.mean(rr)), random_q95=float(np.quantile(rr, 0.95)),
                                     frac_random_ge=float((np.array(rr) >= r_prog).mean()))
        print(h, len(ix), "program", round(r_prog, 3), "| random mean", round(float(np.mean(rr)), 3), "q95", round(float(np.quantile(rr, 0.95)), 3), flush=True)
    res["random_overall_mean"] = float(np.mean(allr))
    json.dump(res, open(out + "/metrics_tahoe_tone_check.json", "w"), indent=1)
    print("random overall mean", round(float(np.mean(allr)), 3), "saved", flush=True)


if __name__ == "__main__":
    main()
