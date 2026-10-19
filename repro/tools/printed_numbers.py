"""Numbers printed in the manuscript (paper/manuscript_v2.md), each with the value it should equal in the result files.
Usage: printed_numbers.py <figdata_dir> [<figdata_dir_2>]  -> table: item | text | value(dir 1) | value(dir 2) | rounded match"""
import sys, os, json
import numpy as np


def J(d, f):
    return json.load(open(os.path.join(d, f)))


def npz_inert_share(d, donor, design, state):
    z = np.load(os.path.join(d, "cd4_%s.npz" % donor)); m = z["inert"]
    return float((z["cross_%s_%s" % (design, state)][m] > 0).mean())


def n_inert(d, donor):
    return int(np.load(os.path.join(d, "cd4_%s.npz" % donor))["inert"].sum())


def trbc_r(d, donor, design):
    z = np.load(os.path.join(d, "cd4_%s.npz" % donor)); ok = (z["trbc_nA_%s_rest" % design] >= 15) & (z["trbc_nB_%s_rest" % design] >= 15)
    a, b = z["trbc_dA_%s_rest" % design][ok], z["trbc_dB_%s_rest" % design][ok]; ok2 = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[ok2], b[ok2])[0, 1])


def dmso(d, key, agg):
    v = [p[key] for p in J(d, "dmso.json")["plates"].values() if key in p]
    return {"median": float(np.median(v)), "min": float(min(v)), "max": float(max(v))}[agg]


def comp(d, donor, st, rk, part):
    c = J(d, "clone_%s.json" % donor)[st]["depth"]; gd = J(d, "guide_donor.json")
    R, G = c["random_half"]["expressed"]["cross_per_target"], c["guide_half"]["expressed"]["cross_per_target"]
    rep = G / R; g = gd["guide_specific_over_diff_guide_" + rk] * rep
    return {"guide": g, "transduction": 1 - rep - g}[part]


def sim(d, k, m="readout_share"):
    return J(d, "sim_v6.json")[k][m]


def ratio_sim(d, sa, kind):
    s = J(d, "sim_v6.json"); a, b = s["atten_sa%s_%s|random" % (sa, kind)], s["atten_sa%s_%s|independent" % (sa, kind)]
    return a["inert_cross_per_target"] / b["expressed_cross_per_target"]


def unit_ratio(d, donor):
    c = J(d, "clone_%s.json" % donor)["rest"]["depth"]
    return c["random_half"]["inert"]["cross_per_target"] / c["guide_half"]["expressed"]["cross_per_target"]


def tcorr(d, p):
    t = J(d, "tahoe_recompute.json")
    return t["cross_plate_halves_twoway"][p + "_cross"] / (t["old_within_well_twoway"][p + "_cross"] * t["within_well_plate14_twoway"][p + "_cross"]) ** 0.5


def ce(d, sub, m, design, agg="mean", kind="ceiling"):
    c = J(d, "celleval.json")["shared"]
    ks = ["random_%d" % s for s in (0, 1, 2)] if design == "random" else ["guide_0real_1pred_0", "guide_1real_0pred_0"]
    x = [c[k][sub][kind].get(m) for k in ks]; x = [np.nan if v is None else v for v in x]
    return {"mean": float(np.nanmean(x)), "min": float(np.nanmin(x)), "max": float(np.nanmax(x))}[agg]


def ce_pair(d, sub, m):
    return "%.2f vs %.2f" % (ce(d, sub, m, "random"), ce(d, sub, m, "guide"))


def clnt(d, key):
    v = [x[key] for f in ("cl_ntc_replogle.json", "cl_ntc_nadig.json", "cl_ntc_h1.json") for k, x in J(d, f).items() if k != "rules"]
    return min(v), max(v)


def tr(d, donor, st, path):
    x = J(d, "clone_%s.json" % donor)[st]
    for p in path.split("/"):
        x = x[p]
    return x


def g(d, f, path):
    x = J(d, f)
    for p in path.split("/"):
        x = x[int(p)] if isinstance(x, list) else x[p]
    return x


ROWS = [
    # (item, text, extractor, format: 'pct1', 'pct0', 'f2', 'f3', 'f1', 'int')
    ("Unexpressed targets, donor 1 / donor 2 (n)", "1,239 / 872", lambda d: "%s / %s" % (format(n_inert(d, "D1"), ","), n_inert(d, "D2")), "str"),
    ("Null positive, random halves, D1 rest", "98.9%", lambda d: npz_inert_share(d, "D1", "random", "rest"), "pct1"),
    ("Null positive, random halves, D1 stim", "96.2%", lambda d: npz_inert_share(d, "D1", "random", "stim"), "pct1"),
    ("Null positive, random halves, D2 rest", "99.2%", lambda d: npz_inert_share(d, "D2", "random", "rest"), "pct1"),
    ("Null positive, random halves, D2 stim", "98.5%", lambda d: npz_inert_share(d, "D2", "random", "stim"), "pct1"),
    ("Null positive, guide halves, D1 rest", "48.6%", lambda d: npz_inert_share(d, "D1", "guide", "rest"), "pct1"),
    ("Null positive, guide halves, D1 stim", "53.7%", lambda d: npz_inert_share(d, "D1", "guide", "stim"), "pct1"),
    ("Null positive, guide halves, D2 rest", "52.2%", lambda d: npz_inert_share(d, "D2", "guide", "rest"), "pct1"),
    ("Null positive, guide halves, D2 stim", "49.9%", lambda d: npz_inert_share(d, "D2", "guide", "stim"), "pct1"),
    ("Donor halves, null positive, rest", "55%", lambda d: g(d, "donor_split.json", "C_null_baseline_donor_split/inert/rest_frac_pos"), "pct0"),
    ("Donor halves, null positive, stim", "53%", lambda d: g(d, "donor_split.json", "C_null_baseline_donor_split/inert/stim_frac_pos"), "pct0"),
    ("NTC random halves positive, min", "96.6%", lambda d: min(v["random_frac_pos"] for k, v in J(d, "ntc_null.json").items() if k != "rules"), "pct1"),
    ("NTC random halves positive, max", "98.8%", lambda d: max(v["random_frac_pos"] for k, v in J(d, "ntc_null.json").items() if k != "rules"), "pct1"),
    ("NTC guides per donor-state, min-max", "578–785", lambda d: "%d–%d" % (min(v["n_guides"] for k, v in J(d, "ntc_null.json").items() if k != "rules"), max(v["n_guides"] for k, v in J(d, "ntc_null.json").items() if k != "rules")), "str"),
    ("NTC different guides positive, min", "48.8%", lambda d: min(v["pairs_frac_pos_mean"] for k, v in J(d, "ntc_null.json").items() if k != "rules"), "pct1"),
    ("NTC different guides positive, max", "49.8%", lambda d: max(v["pairs_frac_pos_mean"] for k, v in J(d, "ntc_null.json").items() if k != "rules"), "pct1"),
    ("DMSO between-half Pearson, median over plates", "0.53", lambda d: dmso(d, "pearson_median", "median"), "f2"),
    ("DMSO Pearson range across plates, min", "0.32", lambda d: dmso(d, "pearson_median", "min"), "f2"),
    ("DMSO Pearson range across plates, max", "0.79", lambda d: dmso(d, "pearson_median", "max"), "f2"),
    ("DMSO well effect per well / drug, median", "0.19", lambda d: dmso(d, "well_per_well_over_drug_median", "median"), "f2"),
    ("DMSO well effect range, min", "0.08", lambda d: dmso(d, "well_per_well_over_drug_median", "min"), "f2"),
    ("DMSO well effect range, max", "0.76", lambda d: dmso(d, "well_per_well_over_drug_median", "max"), "f2"),
    ("DMSO share common to all lines", "76%", lambda d: g(d, "dmso.json", "summary/common_share_of_cross"), "pct0"),
    ("DMSO cross-plate correlation, median", "0.35", lambda d: g(d, "dmso.json", "cross_plate_corr/median"), "f2"),
    ("TRBC var(z), D1 rest", "34", lambda d: tr(d, "D1", "rest", "trbc/overdispersion/var_z"), "f0"),
    ("TRBC var(z), D2 rest", "81", lambda d: tr(d, "D2", "rest", "trbc/overdispersion/var_z"), "f0"),
    ("TRBC var(z) permuted, D1 rest", "1.0", lambda d: tr(d, "D1", "rest", "trbc/overdispersion/var_z_permuted"), "f1"),
    ("TRBC deviation r, random halves, D1 (npz, Fig. 3c)", "0.93", lambda d: trbc_r(d, "D1", "random"), "f2"),
    ("TRBC deviation r, random halves, D2", "0.96", lambda d: trbc_r(d, "D2", "random"), "f2"),
    ("TRBC deviation r, guide halves, D1", "−0.01", lambda d: trbc_r(d, "D1", "guide"), "f2"),
    ("TRBC deviation r, guide halves, D2", "0.02", lambda d: trbc_r(d, "D2", "guide"), "f2"),
    ("Clone skew Spearman (per-state def.), D1 rest", "0.229", lambda d: tr(d, "D1", "rest", "trbc/link_inert/spearman_skew_vs_expr_cross"), "f3"),
    ("Clone skew Spearman (per-state def.), D2 rest", "0.227", lambda d: tr(d, "D2", "rest", "trbc/link_inert/spearman_skew_vs_expr_cross"), "f3"),
    ("Clone skew Spearman (per-state def.), D1 stim", "0.159", lambda d: tr(d, "D1", "stim48", "trbc/link_inert/spearman_skew_vs_expr_cross"), "f3"),
    ("Clone skew Spearman (per-state def.), D2 stim", "0.160", lambda d: tr(d, "D2", "stim48", "trbc/link_inert/spearman_skew_vs_expr_cross"), "f3"),
    ("TRBC valid cells, D1 rest (millions)", "1.75", lambda d: tr(d, "D1", "rest", "trbc/cells_valid") / 1e6, "f2"),
    ("TRBC valid cells, D1 stim (millions)", "1.58", lambda d: tr(d, "D1", "stim48", "trbc/cells_valid") / 1e6, "f2"),
    ("Same guide across donors, null positive rest", "57%", lambda d: g(d, "guide_donor.json", "same_guide/inert_frac_pos_rest"), "pct0"),
    ("Same guide across donors, null positive stim", "54%", lambda d: g(d, "guide_donor.json", "same_guide/inert_frac_pos_stim"), "pct0"),
    ("Different guides across donors, null positive rest", "51%", lambda d: g(d, "guide_donor.json", "diff_guide/inert_frac_pos_rest"), "pct0"),
    ("Different guides across donors, null positive stim", "51%", lambda d: g(d, "guide_donor.json", "diff_guide/inert_frac_pos_stim"), "pct0"),
    ("K562 unexpressed recurrence across screens", "0.85", lambda d: g(d, "k562_inert_screens.json", "inert/guide_fraction"), "f2"),
    ("K562 unexpressed recurrence CI low", "0.75", lambda d: g(d, "k562_inert_screens.json", "inert/ci95/0"), "f2"),
    ("K562 unexpressed recurrence CI high", "0.92", lambda d: g(d, "k562_inert_screens.json", "inert/ci95/1"), "f2"),
    ("K562 unexpressed targets (n)", "178", lambda d: g(d, "k562_inert_screens.json", "inert/n"), "int"),
    ("K562 expressed targets (n)", "651", lambda d: g(d, "k562_inert_screens.json", "expressed/n"), "int"),
    ("Depth: random quarter cross, null, D1 rest", "6.7", lambda d: tr(d, "D1", "rest", "depth/random_quarter/inert/cross_per_target"), "f1"),
    ("Depth: random half cross, null, D1 rest", "6.8", lambda d: tr(d, "D1", "rest", "depth/random_half/inert/cross_per_target"), "f1"),
    ("Depth: random quarter cross, null, D2 rest", "12.6", lambda d: tr(d, "D2", "rest", "depth/random_quarter/inert/cross_per_target"), "f1"),
    ("Depth: random half cross, null, D2 rest", "12.6", lambda d: tr(d, "D2", "rest", "depth/random_half/inert/cross_per_target"), "f1"),
    ("Depth: Pearson quarter, null, D1 rest", "0.14", lambda d: tr(d, "D1", "rest", "depth/random_quarter/inert/pearson_mean"), "f2"),
    ("Depth: Pearson half, null, D1 rest", "0.24", lambda d: tr(d, "D1", "rest", "depth/random_half/inert/pearson_mean"), "f2"),
    ("Depth: Pearson quarter, null, D2 rest", "0.27", lambda d: tr(d, "D2", "rest", "depth/random_quarter/inert/pearson_mean"), "f2"),
    ("Depth: Pearson half, null, D2 rest", "0.40", lambda d: tr(d, "D2", "rest", "depth/random_half/inert/pearson_mean"), "f2"),
    ("Depth: SB from quarter, D1 rest", "0.25", lambda d: tr(d, "D1", "rest", "depth/random_quarter/inert/pearson_sb"), "f2"),
    ("Expressed energy, random halves, D1 rest", "8.9", lambda d: tr(d, "D1", "rest", "depth/random_half/expressed/cross_per_target"), "f1"),
    ("Expressed energy, guide halves, D1 rest", "1.7", lambda d: tr(d, "D1", "rest", "depth/guide_half/expressed/cross_per_target"), "f1"),
    ("Expressed energy, random halves, D2 rest", "16.3", lambda d: tr(d, "D2", "rest", "depth/random_half/expressed/cross_per_target"), "f1"),
    ("Expressed energy, guide halves, D2 rest", "1.7", lambda d: tr(d, "D2", "rest", "depth/guide_half/expressed/cross_per_target"), "f1"),
    ("Same guide, expressed signal, rest", "3.09", lambda d: g(d, "guide_donor.json", "same_guide/expressed_mean_rest"), "f2"),
    ("Different guides, expressed signal, rest", "1.71", lambda d: g(d, "guide_donor.json", "diff_guide/expressed_mean_rest"), "f2"),
    ("Same guide, expressed signal, stim", "2.85", lambda d: g(d, "guide_donor.json", "same_guide/expressed_mean_stim"), "f2"),
    ("Different guides, expressed signal, stim", "1.92", lambda d: g(d, "guide_donor.json", "diff_guide/expressed_mean_stim"), "f2"),
    ("Guide-donor targets (n)", "6,450", lambda d: format(g(d, "guide_donor.json", "n_targets"), ","), "str"),
    ("Guide-specific ratio, rest", "0.80", lambda d: g(d, "guide_donor.json", "guide_specific_over_diff_guide_rest"), "f2"),
    ("Guide-specific ratio, stim", "0.48", lambda d: g(d, "guide_donor.json", "guide_specific_over_diff_guide_stim"), "f2"),
    ("Guide-specific share of random-half signal, D1 rest", "15%", lambda d: comp(d, "D1", "rest", "rest", "guide"), "pct0"),
    ("Guide-specific share, D1 stim", "8%", lambda d: comp(d, "D1", "stim48", "stim", "guide"), "pct0"),
    ("Guide-specific share, D2 rest", "8%", lambda d: comp(d, "D2", "rest", "rest", "guide"), "pct0"),
    ("Transduction-level share, D1 rest", "66%", lambda d: comp(d, "D1", "rest", "rest", "transduction"), "pct0"),
    ("Transduction-level share, D1 stim", "77%", lambda d: comp(d, "D1", "stim48", "stim", "transduction"), "pct0"),
    ("Transduction-level share, D2 rest", "81%", lambda d: comp(d, "D2", "rest", "rest", "transduction"), "pct0"),
    ("K562 expressed recurrence", "0.84 (84%)", lambda d: g(d, "k562_inert_screens.json", "expressed/guide_fraction"), "f2"),
    ("K562 expressed CI", "83–86%", lambda d: "%.0f–%.0f%%" % (100 * g(d, "k562_inert_screens.json", "expressed/ci95/0"), 100 * g(d, "k562_inert_screens.json", "expressed/ci95/1")), "str"),
    ("K562 bound, one screen", "≤16%", lambda d: 1 - g(d, "k562_inert_screens.json", "expressed/guide_fraction"), "pct0"),
    ("K562 additive clonal part", "≈4%", lambda d: (1 - g(d, "k562_inert_screens.json", "inert/guide_fraction")) * (g(d, "k562_inert_screens.json", "inert/within_K6") * g(d, "k562_inert_screens.json", "inert/within_K8")) ** 0.5 / g(d, "k562_inert_screens.json", "inert/n") / ((g(d, "k562_inert_screens.json", "expressed/within_K6") * g(d, "k562_inert_screens.json", "expressed/within_K8")) ** 0.5 / g(d, "k562_inert_screens.json", "expressed/n")), "pct0"),
    ("Tahoe within-well / cross-plate energy (gamma / beta)", "1.4–1.5", lambda d: "%.2f / %.2f" % tuple(1 / tcorr(d, p) for p in ("gamma", "beta")), "str"),
    ("Tahoe cross-plate signal correlation gamma / beta / alpha (ledger)", "0.65 / 0.73 / 0.78", lambda d: "%.2f / %.2f / %.2f" % tuple(tcorr(d, p) for p in ("gamma", "beta", "alpha")), "str"),
    ("Wang reliable, random, unexpressed, D1", "18%", lambda d: g(d, "ci_cd4.json", "wang_reliable_share_rest/random/inert"), "pct0"),
    ("Wang reliable, random, expressed, D1", "24%", lambda d: g(d, "ci_cd4.json", "wang_reliable_share_rest/random/expressed"), "pct0"),
    ("Wang reliable, guide, unexpressed, D1", "0.6%", lambda d: g(d, "ci_cd4.json", "wang_reliable_share_rest/guide/inert"), "pct1"),
    ("Wang reliable, guide, expressed, D1", "2.2%", lambda d: g(d, "ci_cd4.json", "wang_reliable_share_rest/guide/expressed"), "pct1"),
    ("Wang reliable, random, unexpressed, D2", "42%", lambda d: g(d, "ci_cd4_D2.json", "wang_reliable_share_rest/random/inert"), "pct0"),
    ("Wang reliable, random, expressed, D2", "52%", lambda d: g(d, "ci_cd4_D2.json", "wang_reliable_share_rest/random/expressed"), "pct0"),
    ("Wang reliable, guide, unexpressed, D2", "0.9%", lambda d: g(d, "ci_cd4_D2.json", "wang_reliable_share_rest/guide/inert"), "pct1"),
    ("Wang reliable, guide, expressed, D2", "2.3%", lambda d: g(d, "ci_cd4_D2.json", "wang_reliable_share_rest/guide/expressed"), "pct1"),
    ("Tahoe Pearson ceiling, plate 6", "0.87", lambda d: g(d, "ceiling_tahoe.json", "within_well_plate6/pearson_delta/ceiling"), "f2"),
    ("Tahoe Pearson ceiling, plate 14", "0.77", lambda d: g(d, "ceiling_tahoe.json", "within_well_plate14/pearson_delta/ceiling"), "f2"),
    ("Tahoe Pearson ceiling, across plates", "0.64", lambda d: g(d, "ceiling_tahoe.json", "across_plates/pearson_delta/ceiling"), "f2"),
    ("Tahoe discrimination ceiling, within (6 / 14)", "0.99", lambda d: "%.2f / %.2f" % (g(d, "ceiling_tahoe.json", "within_well_plate6/discrimination_l1_within_line/ceiling"), g(d, "ceiling_tahoe.json", "within_well_plate14/discrimination_l1_within_line/ceiling")), "str"),
    ("Tahoe discrimination ceiling, across (both directions)", "0.89–0.92", lambda d: "%.2f–%.2f" % (g(d, "ceiling_tahoe.json", "across_plates/discrimination_l1_within_line/ceiling"), g(d, "ceiling_tahoe.json", "across_plates_rev/discrimination_l1_within_line/ceiling")), "str"),
    ("CD4 context-specific share, random, D1 (CI)", "0.30 (0.30–0.31)", lambda d: "%.2f (%.2f–%.2f)" % tuple(g(d, "ci_cd4.json", "random_r1/gamma_share")), "str"),
    ("CD4 context-specific share, guide, D1 (CI)", "0.34 (0.31–0.38)", lambda d: "%.2f (%.2f–%.2f)" % tuple(g(d, "ci_cd4.json", "guide/gamma_share")), "str"),
    ("Context-specific energy ratio random / guide", "7.6", lambda d: g(d, "ci_cd4.json", "random_r1/gamma_cross_per_target") / g(d, "ci_cd4.json", "guide/gamma_cross_per_target"), "f1"),
    ("Context-specific energy, random / guide", "2.58 / 0.34", lambda d: "%.2f / %.2f" % (g(d, "ci_cd4.json", "random_r1/gamma_cross_per_target"), g(d, "ci_cd4.json", "guide/gamma_cross_per_target")), "str"),
    ("Shared energy, random / guide", "5.88 / 0.65", lambda d: "%.2f / %.2f" % (g(d, "ci_cd4.json", "random_r1/beta_cross_per_target"), g(d, "ci_cd4.json", "guide/beta_cross_per_target")), "str"),
    ("Shared energy ratio", "9", lambda d: g(d, "ci_cd4.json", "random_r1/beta_cross_per_target") / g(d, "ci_cd4.json", "guide/beta_cross_per_target"), "f0"),
    ("Tahoe context-specific share, within wells (6 / 14)", "0.39", lambda d: "%.3f / %.3f" % (g(d, "ci_tahoe.json", "gamma_share_within_well_plate6/0"), g(d, "ci_tahoe.json", "gamma_share_within_well_plate14/0")), "str"),
    ("Tahoe context-specific share, across plates (CI)", "0.36 (0.25–0.40)", lambda d: "%.2f (%.2f–%.2f)" % tuple(g(d, "ci_tahoe.json", "gamma_share_across_plates")), "str"),
    ("Tahoe predictor, same well, target plate 14 (CI)", "0.33 (0.28–0.41)", lambda d: "%.2f (%.2f–%.2f)" % tuple(g(d, "ci_tahoe.json", "readout_target_plate14/same_well")), "str"),
    ("Tahoe predictor, same well, target plate 6 (CI)", "0.30 (0.24–0.36)", lambda d: "%.2f (%.2f–%.2f)" % tuple(g(d, "ci_tahoe.json", "readout_target_plate6/same_well")), "str"),
    ("Tahoe predictor, other plate, target 14 (CI)", "0.04 (−0.05 to 0.07)", lambda d: "%.2f (%.2f to %.2f)" % tuple(g(d, "ci_tahoe.json", "readout_target_plate14/other_plate")), "str"),
    ("Tahoe predictor, other plate, target 6 (CI)", "0.09 (0.03–0.14)", lambda d: "%.2f (%.2f–%.2f)" % tuple(g(d, "ci_tahoe.json", "readout_target_plate6/other_plate")), "str"),
    ("CD4 predictor, random, D1 (CI)", "0.47 (0.46–0.48)", lambda d: "%.2f (%.2f–%.2f)" % tuple(g(d, "ci_cd4.json", "random_r1/readout_k20")), "str"),
    ("CD4 predictor, random, D2 (CI)", "0.52 (0.51–0.52)", lambda d: "%.2f (%.2f–%.2f)" % tuple(g(d, "ci_cd4_D2.json", "random_r1/readout_k20")), "str"),
    ("CD4 predictor, guide, D1 (CI)", "0.004 (0.000–0.007)", lambda d: "%.3f (%.3f–%.3f)" % tuple(g(d, "ci_cd4.json", "guide/readout_k20")), "str"),
    ("CD4 predictor, guide, D2 (CI)", "−0.002 (−0.007 to 0.002)", lambda d: "%.3f (%.3f to %.3f)" % tuple(g(d, "ci_cd4_D2.json", "guide/readout_k20")), "str"),
    ("CD4 predictor, three random halvings, D1", "0.468–0.471", lambda d: "%.3f–%.3f" % (min(g(d, "ci_cd4.json", k + "/readout_k20/0") for k in ("random_r1", "random_r2", "random_xor")), max(g(d, "ci_cd4.json", k + "/readout_k20/0") for k in ("random_r1", "random_r2", "random_xor"))), "str"),
    ("Strength-matched predictor, random, D1 (k=20)", "0.24", lambda d: g(d, "kd_D1.json", "kd_matched|cell/readout_k20"), "f2"),
    ("Strength-matched predictor, random, D2 (k=20)", "0.34", lambda d: g(d, "kd_D2.json", "kd_matched|cell/readout_k20"), "f2"),
    ("Strength-matched targets D1 / D2", "615 / 470", lambda d: "%d / %d" % (g(d, "kd_D1.json", "n_kd_matched"), g(d, "kd_D2.json", "n_kd_matched")), "str"),
    ("Strength-matched predictor k=10, D1 / D2", "0.20 / 0.29", lambda d: "%.2f / %.2f" % (g(d, "kd_D1.json", "kd_matched|cell/readout_k10"), g(d, "kd_D2.json", "kd_matched|cell/readout_k10")), "str"),
    ("Strength-matched predictor, guide, D1 / D2", "≈0", lambda d: "%.3f / %.3f" % (g(d, "kd_D1.json", "kd_matched|guide/readout_k20"), g(d, "kd_D2.json", "kd_matched|guide/readout_k20")), "str"),
    ("Same-guide predictor across donors", "0.012", lambda d: g(d, "guide_donor.json", "same_guide/readout_k20"), "f3"),
    ("Different-guide predictor across donors", "0.0002", lambda d: g(d, "guide_donor.json", "diff_guide/readout_k20"), "f4"),
    ("K562-RPE1 predictor, same pair (targets K8, K6; k 10-20)", "0.57–0.63", lambda d: "%.2f–%.2f" % (min(g(d, "replogle_indep.json", "target_%s/pair_same_k%d" % (t, k)) for t in ("K8", "K6") for k in (10, 20)), max(g(d, "replogle_indep.json", "target_%s/pair_same_k%d" % (t, k)) for t in ("K8", "K6") for k in (10, 20))), "str"),
    ("K562-RPE1 predictor, independent K562", "0.57–0.62", lambda d: "%.2f–%.2f" % (min(g(d, "replogle_indep.json", "target_%s/pair_indepK_k%d" % (t, k)) for t in ("K8", "K6") for k in (10, 20)), max(g(d, "replogle_indep.json", "target_%s/pair_indepK_k%d" % (t, k)) for t in ("K8", "K6") for k in (10, 20))), "str"),
    ("Model: fabricated predictor, shared units", "0.10", lambda d: sim(d, "fab_low_shared|random"), "f2"),
    ("Model: fabricated, independent / separate", "0.00", lambda d: "%.2f / %.2f" % (sim(d, "fab_low_shared|independent"), sim(d, "fab_low_separate|random")), "str"),
    ("Model: true explained share (attenuation)", "0.70", lambda d: J(d, "sim_v6.json")["truth"]["readout_share"], "f2"),
    ("Unit-to-signal ratio, D1", "4.0", lambda d: unit_ratio(d, "D1"), "f1"),
    ("Unit-to-signal ratio, D2", "7.5", lambda d: unit_ratio(d, "D2"), "f1"),
    ("Model ratios at sigma_a 1.5 (low / iso)", "4.1–4.2", lambda d: "%.1f / %.1f" % (ratio_sim(d, "1.5", "low"), ratio_sim(d, "1.5", "iso")), "str"),
    ("Model ratios at sigma_a 2.0 (low / iso)", "7.2–7.6", lambda d: "%.1f / %.1f" % (ratio_sim(d, "2.0", "low"), ratio_sim(d, "2.0", "iso")), "str"),
    ("Model explained share at sigma_a 1.5, along response", "0.10", lambda d: sim(d, "atten_sa1.5_low|independent"), "f2"),
    ("Model explained share at sigma_a 2.0, along response", "0.03", lambda d: sim(d, "atten_sa2.0_low|independent"), "f2"),
    ("Tahoe signal correlation across plates (beta)", "0.73", lambda d: tcorr(d, "beta"), "f2"),
    ("K562 template artefact, inert-referenced -> template removed", "83% → 89%", lambda d: "%.0f%% → %.0f%%" % (100 * g(d, "k562_template.json", "random/inert_referenced/inert/frac_pos"), 100 * g(d, "k562_template.json", "random/inert/frac_pos")), "str"),
    ("cell-eval, all 300, Pearson ceiling random vs guide", "0.32 vs 0.06", lambda d: ce_pair(d, "all", "pearson_delta"), "str"),
    ("cell-eval, all 300, DE overlap", "0.31 vs 0.11", lambda d: ce_pair(d, "all", "overlap_at_N"), "str"),
    ("cell-eval, all 300, DE direction match", "0.84 vs 0.63", lambda d: ce_pair(d, "all", "de_direction_match"), "str"),
    ("cell-eval, all 300, discrimination (cosine)", "0.98 vs 0.69", lambda d: ce_pair(d, "all", "discrimination_score_cosine"), "str"),
    ("cell-eval, strong knockdowns, Pearson", "0.38 vs 0.13", lambda d: ce_pair(d, "strong", "pearson_delta"), "str"),
    ("cell-eval, unexpressed, Pearson random", "0.26", lambda d: ce(d, "inert", "pearson_delta", "random"), "f2"),
    ("cell-eval, unexpressed, Pearson guide split-half r (ceiling n.d.)", "−0.01", lambda d: ce(d, "inert", "pearson_delta", "guide", kind="raw_mean"), "f2"),
    ("cell-eval, unexpressed, direction match random", "0.81", lambda d: ce(d, "inert", "de_direction_match", "random"), "f2"),
    ("cell-eval, unexpressed, direction match guide (range)", "0.57–0.61", lambda d: "%.2f–%.2f" % (ce(d, "inert", "de_direction_match", "guide", "min"), ce(d, "inert", "de_direction_match", "guide", "max")), "str"),
    ("cell-eval, unexpressed, discrimination (cosine) random", "0.97", lambda d: ce(d, "inert", "discrimination_score_cosine", "random"), "f2"),
    ("cell-eval, unexpressed, discrimination (cosine) guide (range)", "0.63–0.65", lambda d: "%.2f–%.2f" % (ce(d, "inert", "discrimination_score_cosine", "guide", "min"), ce(d, "inert", "discrimination_score_cosine", "guide", "max")), "str"),
    ("cell-eval, perturbations in all runs", "300", lambda d: J(d, "celleval.json")["n_shared_perturbations"], "int"),
    ("Cell-line and H1 NT constructs, random halves positive, min-max", "60–76%", lambda d: "%.0f–%.0f%%" % tuple(100 * x for x in clnt(d, "random_frac_pos")), "str"),
    ("Cell-line and H1 NT constructs, different constructs positive, min-max", "49–51%", lambda d: "%.0f–%.0f%%" % tuple(100 * x for x in clnt(d, "pairs_frac_pos_mean")), "str"),
    ("Cell-line and H1 NT / targeting random-half signal, min-max", "0.6–7%", lambda d: "%.1f–%.0f%%" % tuple(100 * x for x in clnt(d, "ratio_nt_to_target_mean")), "str"),
    ("H1 NT / targeting random-half signal", "0.6%", lambda d: 100 * g(d, "cl_ntc_h1.json", "h1/ratio_nt_to_target_mean"), "f1p"),
    ("H1 cells / targets / NT constructs", "221,273 / 150 / 31", lambda d: "%s / %d / %d" % (format(g(d, "cl_ntc_h1.json", "h1/n_cells"), ","), g(d, "cl_ntc_h1.json", "h1/n_targets"), g(d, "cl_ntc_h1.json", "h1/n_constructs")), "str"),
    ("CD4 NT / targeting random-half signal, D1 / D2 rest", "1.2–1.3", lambda d: "%.2f / %.2f" % (g(d, "cd4_ntc_ratio.json", "D2_rest/ratio_nt_to_target_mean"), g(d, "cd4_ntc_ratio.json", "D1_rest/ratio_nt_to_target_mean")), "str"),
    ("TRBC1 deviation, rest vs stim, D1", "0.89", lambda d: g(d, "trbc_state_keff.json", "D1/pearson_rest_vs_stim"), "f2"),
    ("TRBC1 deviation, rest vs stim, D2", "0.94", lambda d: g(d, "trbc_state_keff.json", "D2/pearson_rest_vs_stim"), "f2"),
    ("TRBC rest-vs-stim perturbations, D1 / D2", "9,830 / 9,037", lambda d: "%s / %s" % (format(g(d, "trbc_state_keff.json", "D1/n_common"), ","), format(g(d, "trbc_state_keff.json", "D2/n_common"), ",")), "str"),
    ("Effective clones, median upper bound, D1 rest", "about 7", lambda d: g(d, "trbc_state_keff.json", "D1/keff_rest/keff_upper_median"), "f0"),
    ("Effective clones, median upper bound, D2 rest", "4", lambda d: g(d, "trbc_state_keff.json", "D2/keff_rest/keff_upper_median"), "f0"),
]


sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v15_rows import ROWS15
ROWS += [(a, b, c, "str") for a, b, c in ROWS15]


def fmt(v, f):
    if isinstance(v, str):
        return v
    if f == "pct1": return "%.1f%%" % (100 * v)
    if f == "pct0": return "%.0f%%" % (100 * v)
    if f == "f0": return "%.0f" % v
    if f == "f1": return "%.1f" % v
    if f == "f1p": return "%.1f%%" % v
    if f == "f2": return "%.2f" % v
    if f == "f3": return "%.3f" % v
    if f == "f4": return "%.4f" % v
    if f == "int": return "%d" % v
    return json.dumps(v)[:80]


def main():
    dirs = sys.argv[1:]
    print("| Item | text | " + " | ".join(os.path.basename(os.path.normpath(d)) for d in dirs) + " |")
    print("|---|---|" + "---|" * len(dirs))
    for item, txt, fn, f in ROWS:
        vals = []
        for d in dirs:
            try:
                vals.append(fmt(fn(d), f))
            except Exception as e:
                vals.append("ERR %s" % str(e)[:40])
        print("| %s | %s | %s |" % (item, txt, " | ".join(vals)))


if __name__ == "__main__":
    main()
