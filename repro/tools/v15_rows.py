"""Rows for the v15 wording (imported by printed_numbers.py). Each extractor takes a folder holding the figure data
plus methods_numbers.json (intermediate results), celleval_header.txt (cell-eval run header) and plate_overlap.txt."""
import os, json
import numpy as np
from scipy.stats import spearmanr


def _J(d, f):
    return json.load(open(os.path.join(d, f)))


def _npz(d, donor):
    return np.load(os.path.join(d, "cd4_%s.npz" % donor))


def share(d, design, donor, st):
    z = _npz(d, donor); m = z["inert"]
    return float((z["cross_%s_%s" % (design, st)][m] > 0).mean())


def rng(vals, fmt="%.0f–%.0f%%", mul=100):
    return fmt % (mul * min(vals), mul * max(vals))


def skew_link(d, donor):
    z = _npz(d, donor); m = z["inert"] & (z["trbc_n_rest"] >= 30)
    return spearmanr(np.abs(z["trbc_z_rest"][m]), z["cross_random_rest"][m])[0], int(m.sum())


def n_trbc(d, donor, k):
    z = _npz(d, donor)
    if k == "hist":
        return int((z["trbc_n_rest"] >= 30).sum())
    a, b = z["trbc_nA_%s_rest" % k], z["trbc_nB_%s_rest" % k]
    return int(((a >= 15) & (b >= 15)).sum())


def clone(d, donor, path):
    x = _J(d, "clone_%s.json" % donor)
    for p in path.split("/"):
        x = x[p]
    return x


def ci_halfwidth(d, design):
    w = []
    for f in ("ci_cd4.json", "ci_cd4_D2.json"):
        c = _J(d, f)
        for k in ("inert_frac_pos_rest", "inert_frac_pos_stim"):
            v = c["random_r1" if design == "random" else "guide"][k]; w.append(max(v[0] - v[1], v[2] - v[0]))
    return 100 * max(w)


def M(d, k):
    return _J(d, "methods_numbers.json")[k]


def celleval_median(d):
    return float(open(os.path.join(d, "celleval_header.txt")).read().split("median cells per guide side")[1].split()[0])


def _r(d, donor):
    z = _npz(d, donor); a, b = z["trbc_nA_guide_rest"], z["trbc_nB_guide_rest"]; ok = (a >= 15) & (b >= 15)
    x, y = z["trbc_dA_guide_rest"][ok], z["trbc_dB_guide_rest"][ok]; k = np.isfinite(x) & np.isfinite(y)
    return float(np.corrcoef(x[k], y[k])[0, 1])


def _comp(d, donor, st, rk, part):
    c = _J(d, "clone_%s.json" % donor)[st]["depth"]; gd = _J(d, "guide_donor.json")
    rep = c["guide_half"]["expressed"]["cross_per_target"] / c["random_half"]["expressed"]["cross_per_target"]
    g = gd["guide_specific_over_diff_guide_" + rk] * rep
    return {"guide": g, "transduction": 1 - rep - g}[part]


def _tahoe_readout(d, which):
    t = _J(d, "ci_tahoe.json"); v = [t["readout_target_plate%s" % p][which][0] for p in ("6", "14")]
    return "%.2f–%.2f" % (min(v), max(v))


def _robust_tahoe(d):
    t = _J(d, "robust_tahoe.json")
    same = [x for k, v in t.items() if k != "rules" for kk, x in v.items() if kk.startswith("same_well") and "shuffled" not in kk]
    other = [x for k, v in t.items() if k != "rules" for kk, x in v.items() if kk.startswith("other_plate") and "shuffled" not in kk]
    return "%.2f–%.2f / %.2f–%.2f" % (min(same), max(same), min(other), max(other))


def _gamma_tahoe(d):
    t = _J(d, "ci_tahoe.json")
    v = [t[k][0] for k in ("gamma_share_across_plates", "gamma_share_within_well_plate6", "gamma_share_within_well_plate14")]
    return "%.2f–%.2f" % (min(v), max(v))


def _tahoe_ceil(d):
    t = _J(d, "ceiling_tahoe.json"); w = [t["within_well_plate%s" % p]["pearson_delta"]["ceiling"] for p in ("6", "14")]
    return "%.2f–%.2f against %.2f" % (min(w), max(w), t["across_plates"]["pearson_delta"]["ceiling"])


DS = (("D1", "rest", "rest"), ("D1", "stim48", "stim"), ("D2", "rest", "rest"))

ROWS15 = [
    ("Null positive, random halves, all donors and states", "96–99%", lambda d: rng([share(d, "random", D, s) for D in ("D1", "D2") for s in ("rest", "stim")])),
    ("Null positive, guide halves, all donors and states", "49–54%", lambda d: rng([share(d, "guide", D, s) for D in ("D1", "D2") for s in ("rest", "stim")])),
    ("Null positive, donor halves, rest and stim", "53–55%", lambda d: rng([_J(d, "donor_split.json")["C_null_baseline_donor_split"]["inert"][k] for k in ("rest_frac_pos", "stim_frac_pos")])),
    ("NT guides positive, random halves, CD4", "97–99%", lambda d: rng([v["random_frac_pos"] for k, v in _J(d, "ntc_null.json").items() if k != "rules"])),
    ("NT guides positive, different guides, CD4", "49–50%", lambda d: rng([v["pairs_frac_pos_mean"] for k, v in _J(d, "ntc_null.json").items() if k != "rules"])),
    ("TRBC var(z), D1 / D2", "34 / 81", lambda d: "%.0f / %.0f" % (clone(d, "D1", "rest/trbc/overdispersion/var_z"), clone(d, "D2", "rest/trbc/overdispersion/var_z"))),
    ("TRBC deviation r, guide halves, D1 / D2 (printed as r ≈ 0)", "≈ 0", lambda d: "%.2f / %.2f" % (_r(d, "D1"), _r(d, "D2"))),
    ("Clone-skew Spearman, D1 / D2 (Fig. 3d)", "0.21 and 0.24", lambda d: "%.2f and %.2f" % (skew_link(d, "D1")[0], skew_link(d, "D2")[0])),
    ("Fig. 3d targets, D1 / D2", "1,234 / 864", lambda d: "%s / %d" % (format(skew_link(d, "D1")[1], ","), skew_link(d, "D2")[1])),
    ("K562 unexpressed recurrence across screens", "85%", lambda d: "%.0f%%" % (100 * _J(d, "k562_inert_screens.json")["inert"]["guide_fraction"])),
    ("Expressed energy lost under guide halves, D1 / D2", "about 80% / 90%", lambda d: "%.0f%% / %.0f%%" % tuple(100 * (1 - clone(d, D, "rest/depth/guide_half/expressed/cross_per_target") / clone(d, D, "rest/depth/random_half/expressed/cross_per_target")) for D in ("D1", "D2"))),
    ("Guide-specific share, range", "8–15%", lambda d: rng([_comp(d, D, s, r, "guide") for D, s, r in DS])),
    ("Transduction-level share, range", "66–81%", lambda d: rng([_comp(d, D, s, r, "transduction") for D, s, r in DS])),
    ("Tahoe Pearson ceilings within wells vs across plates", "0.77–0.87 against 0.64", _tahoe_ceil),
    ("Context-specific share, CD4 random / guide (printed as about 0.3)", "about 0.3", lambda d: "%.2f / %.2f" % (_J(d, "ci_cd4.json")["random_r1"]["gamma_share"][0], _J(d, "ci_cd4.json")["guide"]["gamma_share"][0])),
    ("Context-specific share, Tahoe range", "0.36–0.39", _gamma_tahoe),
    ("Energy inflation, context-specific / shared (printed as 8–9 times)", "8–9", lambda d: "%.1f / %.1f" % tuple(_J(d, "ci_cd4.json")["random_r1"][k] / _J(d, "ci_cd4.json")["guide"][k] for k in ("gamma_cross_per_target", "beta_cross_per_target"))),
    ("Tahoe predictor, same wells, range", "0.30–0.33", lambda d: _tahoe_readout(d, "same_well")),
    ("Tahoe predictor, replicate plate, range", "0.04–0.09", lambda d: _tahoe_readout(d, "other_plate")),
    ("CD4 predictor, random halves, D1 / D2 (printed as about 0.5)", "about 0.5", lambda d: "%.2f / %.2f" % (_J(d, "ci_cd4.json")["random_r1"]["readout_k20"][0], _J(d, "ci_cd4_D2.json")["random_r1"]["readout_k20"][0])),
    ("CD4 predictor near 0: guide D1, guide D2, matched D1, matched D2, same guide across donors", "about 0", lambda d: "%.3f / %.3f / %.3f / %.3f / %.3f" % (_J(d, "ci_cd4.json")["guide"]["readout_k20"][0], _J(d, "ci_cd4_D2.json")["guide"]["readout_k20"][0], _J(d, "kd_D1.json")["kd_matched|guide"]["readout_k20"], _J(d, "kd_D2.json")["kd_matched|guide"]["readout_k20"], _J(d, "guide_donor.json")["same_guide"]["readout_k20"])),
    ("Null shares over extractions and halvings, D1 rest", "98.4–99.0%", lambda d: rng(M(d, "null_share_rest_all")[:-1], "%.1f–%.1f%%")),
    ("Null shares over extractions and halvings, D1 stim", "96.2–97.7%", lambda d: rng(M(d, "null_share_stim_all")[:-1], "%.1f–%.1f%%")),
    ("Second / main implementation, CD4 D1 random", "0.446 and 0.469", lambda d: "%.3f and %.3f" % (M(d, "inert2_cell_readout")["/B_template_removed/readout_k20"], M(d, "ci_D1_random_r1_readout"))),
    ("Both implementations, CD4 D1 guide (at most)", "0.004", lambda d: "%.3f" % max(M(d, "inert2_guide_readout")["/B_template_removed/readout_k20"], _J(d, "ci_cd4.json")["guide"]["readout_k20"][0])),
    ("Tahoe second implementation, same wells / replicate plate", "0.23–0.32 / 0.02–0.07", _robust_tahoe),
    ("Effective-clone bound: perturbations with z² > 1, D1 / D2 rest", "85–91%", lambda d: rng([_J(d, "trbc_state_keff.json")[D]["keff_rest"]["frac_z2_gt1"] for D in ("D1", "D2")])),
    ("cell-eval median cells per side", "62", lambda d: "%.0f" % celleval_median(d)),
    ("Fig. 2b bootstrap half-widths, guide / random (percentage points)", "about ±3 / about ±1", lambda d: "%.1f / %.1f" % (ci_halfwidth(d, "guide"), ci_halfwidth(d, "random"))),
    ("Fig. 3b histogram perturbations (D1 rest, >= 30 TRBC cells)", "6,610", lambda d: format(n_trbc(d, "D1", "hist"), ",")),
    ("Fig. 3b var(z) n, D1 / D2", "10,165 / 9,334", lambda d: "%s / %s" % (format(clone(d, "D1", "rest/trbc/overdispersion/n_targets"), ","), format(clone(d, "D2", "rest/trbc/overdispersion/n_targets"), ","))),
    ("Fig. 3c n, random / guide (>= 15 TRBC cells per half)", "6,568 / 6,072", lambda d: "%s / %s" % (format(n_trbc(d, "D1", "random"), ","), format(n_trbc(d, "D1", "guide"), ","))),
    ("Fig. 3c r, guide halves, D1", "−0.01", lambda d: "%.2f" % _r(d, "D1")),
    ("Fig. 3e targets: D1 rest / stim, D2 rest / stim", "909 / 856 / 520 / 523", lambda d: "%d / %d / %d / %d" % (clone(d, "D1", "rest/n_inert"), clone(d, "D1", "stim48/n_inert"), clone(d, "D2", "rest/n_inert"), clone(d, "D2", "stim48/n_inert"))),
    ("Fig. 4b / 5b n: all / unexpressed / expressed", "6,615 / 1,235 / 1,097", lambda d: "%s / %s / %s" % tuple(format(_J(d, "ci_cd4.json")[k], ",") for k in ("n_targets", "n_inert", "n_expressed"))),
    ("Fig. 5d n, CD4 donor 2", "5,003", lambda d: format(_J(d, "ci_cd4_D2.json")["n_targets"], ",")),
    ("Fig. 5d n, K562–RPE1", "1,820", lambda d: format(_J(d, "replogle_indep.json")["n_perturbations"], ",")),
    ("Fig. 2b n, donor halves", "1,479", lambda d: format(_J(d, "donor_split.json")["n_inert"], ",")),
    ("DMSO lines per plate", "47–48", lambda d: "%d–%d" % (min(p["n_lines"] for p in _J(d, "dmso.json")["plates"].values()), max(p["n_lines"] for p in _J(d, "dmso.json")["plates"].values()))),
    ("Tahoe lines × conditions", "48 lines × 91 conditions", lambda d: "%d lines × %d conditions" % (_J(d, "ci_tahoe.json")["n_lines"], _J(d, "ci_tahoe.json")["n_conditions"])),
    ("Tahoe ceilings, reverse direction (Pearson, discrimination)", "0.64 and 0.92", lambda d: "%.2f and %.2f" % (_J(d, "ceiling_tahoe.json")["across_plates_rev"]["pearson_delta"]["ceiling"], _J(d, "ceiling_tahoe.json")["across_plates_rev"]["discrimination_l1_within_line"]["ceiling"])),
    ("Model perturbations", "800", lambda d: "%d" % _J(d, "sim_v6.json")["params"]["P"]),
    ("Tahoe plates 6 / 14 shared conditions; largest other pair", "94; 13", lambda d: open(os.path.join(d, "plate_overlap.txt")).read().strip()),
]
