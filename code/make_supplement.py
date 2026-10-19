"""Supplementary tables and figures of the clone-effect preprint. Reads the result files in <data_dir> (data/ in this
repository; data/README.md lists the script behind each file), fills the placeholders of the text template and draws the
Supplementary Figures in the style of the main figures. Numbering follows first citation in the manuscript (v10);
placeholder names are internal and keep their original labels.
Table S1  non-targeting guides as null perturbations (cd4_ntc_null.py; cell lines cl_ntc_null.py, cd4_ntc_ratio.py)  {TABLE_NTC}, {TABLE_NTC_CL}
Table S2  same vs different guide across donors (cd4_guide_donor.py)              {TABLE_S5}, {S5_META}, {S5_RATIO}
Table S3  cell-eval 0.8.2 ceilings, corrected and uncorrected (run_celleval_ceiling.py)  {TABLE_S3_CORR}, {TABLE_S3_RAW}
Table S4  template-direction pass test (readout_inert_line.py)                    {TABLE_S1}, {S1_SHUF}
Table S5  robustness (robust_cd4_tahoe.py), CD4 and Tahoe                         {TABLE_S2_CD4}, {TABLE_S2_TAHOE}
Table S6  donor-2 expectations (static text in the template)
Fig S1    three-term model scenarios (clone_model_sim.py v5)                      figS1_model
Fig S2    Tahoe DMSO wells per plate (tahoe_dmso_null.py)                         figS2_dmso
Fig S3    K562: unexpressed-target structure across two independent screens; template artifact   figS3_k562
Usage: make_supplement.py <data_dir> <template.md> <out.md> <fig_dir>"""
import sys, os, json
import numpy as np
import plot_main_figs as P
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

RED, BLUE, GREY = P.RED, P.BLUE, P.GREY


def f2(x, nd=2):
    if x is not None and not (isinstance(x, float) and np.isnan(x)) and abs(x) < 0.5 * 10 ** -nd:
        x = 0.0
    return "–" if x is None or (isinstance(x, float) and np.isnan(x)) else ("%.*f" % (nd, x)).replace("-", "−")


def span(lo, hi):
    """a range as an en dash, or with "to" when an end is negative, where a dash would read as a minus sign"""
    return "%s to %s" % (lo, hi) if "−" in lo + hi else "%s–%s" % (lo, hi)


def ci(v, nd=2):
    return "%s (%s)" % (f2(v[0], nd), span(f2(v[1], nd), f2(v[2], nd)))


def rng_(vals, nd=2):
    vals = [v for v in vals if v is not None and not np.isnan(v)]
    if not vals:
        return "–"
    lo, hi = f2(min(vals), nd), f2(max(vals), nd)
    return lo if lo == hi else span(lo, hi)


# ------------------------------------------------------------------------------------------------------------ tables
def table_s1(d):
    rows = []
    def add(name, halves, js, lab):
        for k in ("k20",):
            ii, ee = js[k].get(lab + "_inert", {}), js[k].get(lab + "_expressed", {})
            if "pass_share" not in ii:
                continue
            k10i, k10e = js["k10"][lab + "_inert"]["pass_share"], js["k10"][lab + "_expressed"]["pass_share"]
            rows.append("| %s | %s | %d | %s | %s | %d | %s | %s |" % (
                name, halves, ii["n"], ci([ii["pass_share"]] + ii["ci95"]), f2(k10i), ee["n"],
                ci([ee["pass_share"]] + ee["ci95"]), f2(k10e)))
    for donor in ("D1", "D2"):
        for halves, f in (("random", "pass_cd4_%s_cell" % donor), ("guide", "pass_cd4_%s_guide" % donor)):
            js = P.J(d, f + ".json")
            for lab, nm in (("rest", "resting"), ("stim", "stimulated")):
                add("CD4 donor %s, %s" % (donor[1], nm), halves, js, lab)
    k = P.J(d, "pass_k562_rpe1.json"); add("K562 (K562–RPE1)", "random", k, "k562"); add("RPE1 (K562–RPE1)", "random", k, "rpe1")
    j = P.J(d, "pass_jh.json"); add("Jurkat (Jurkat–HepG2)", "random", j, "jurkat"); add("HepG2 (Jurkat–HepG2)", "random", j, "hepg2")
    shuf = []
    for f in ("pass_cd4_D1_cell", "pass_cd4_D1_guide", "pass_cd4_D2_cell", "pass_cd4_D2_guide", "pass_k562_rpe1", "pass_jh"):
        js = P.J(d, f + ".json")
        for kk in ("k10", "k20"):
            for key, v in js[kk].items():
                if isinstance(v, dict) and "pass_share_shuffled" in v:
                    shuf.append(abs(v["pass_share_shuffled"]))
    head = ("| Context | Halves | Unexpressed, n | Pass, k = 20 (95% CI) | k = 10 | Expressed, n | Pass, k = 20 (95% CI) | k = 10 |\n"
            "|------|---|---|-----|--|---|-----|--|")
    return head + "\n" + "\n".join(rows), max(shuf)


def table_s2(d):
    out = []
    head = ("| Donor | Targets | Genes | Halves | Unexpressed positive, resting / stimulated | Context-specific share | "
            "Explained share, k = 5 / 10 / 20 / 50 | Shuffled, k = 20 |\n|--|----|---|---|----|---|-----|--|")
    out.append(head)
    for donor in ("D1", "D2"):
        r = P.J(d, "robust_%s.json" % donor)
        for ts, tsn in (("all_common", "common (%s)" % format(r["n_targets"], ",")), ("guide_balanced", "guide-balanced (%s)" % format(r["n_balanced"], ","))):
            for gs, gsn in (("all_expressed", "all kept"), ("top2000", "top 2,000")):
                for h, hn in (("cell", "random"), ("guide", "guide")):
                    v = r["%s|%s|%s" % (ts, gs, h)]
                    out.append("| %s | %s | %s | %s | %.1f%% / %.1f%% | %s | %s / %s / %s / %s | %s |" % (
                        donor[1], tsn, gsn, hn, 100 * v["inert_frac_pos_rest"], 100 * v["inert_frac_pos_stim"], f2(v["gamma_share"]),
                        f2(v["readout_k5"], 3), f2(v["readout_k10"], 3), f2(v["readout_k20"], 3), f2(v["readout_k50"], 3), f2(v["readout_k20_shuffled"], 3)))
    cd4 = "\n".join(out)
    t = P.J(d, "robust_tahoe.json")
    rows = ["| Target | Genes | Same wells, k = 5 / 10 / 20 | Replicate plate, k = 5 / 10 / 20 | Largest shuffled value |\n|--|----|----|----|--|"]
    for key, v in t.items():
        if key == "rules":
            continue
        tp, gs = key.split("|")
        sh = max(abs(x) for kk, x in v.items() if kk.endswith("shuffled"))
        rows.append("| plate %s | %s | %s / %s / %s | %s / %s / %s | %s |" % (
            tp.replace("target_plate", ""), "all kept" if gs == "all" else "top 2,000 on the other plate",
            f2(v["same_well_k5"]), f2(v["same_well_k10"]), f2(v["same_well_k20"]),
            f2(v["other_plate_k5"]), f2(v["other_plate_k10"]), f2(v["other_plate_k20"]), "%.4f" % sh))
    return cd4, "\n".join(rows)


METRICS = [("overlap_at_N", "DE overlap (top N)"), ("overlap_at_50", "DE overlap (top 50)"), ("overlap_at_100", "DE overlap (top 100)"),
           ("overlap_at_200", "DE overlap (top 200)"), ("overlap_at_500", "DE overlap (top 500)"), ("precision_at_N", "DE precision (top N)"),
           ("precision_at_50", "DE precision (top 50)"), ("precision_at_100", "DE precision (top 100)"), ("precision_at_200", "DE precision (top 200)"),
           ("precision_at_500", "DE precision (top 500)"), ("de_sig_genes_recall", "DE recall"), ("de_direction_match", "DE direction match*"),
           ("de_spearman_sig", "DE Spearman (significant genes)†"), ("de_spearman_lfc_sig", "DE Spearman of log fold changes†"),
           ("pr_auc", "PR-AUC"), ("roc_auc", "ROC-AUC*"), ("pearson_delta", "Pearson Δ"),
           ("discrimination_score_l1", "Discrimination (L1)*"), ("discrimination_score_l2", "Discrimination (L2)*"),
           ("discrimination_score_cosine", "Discrimination (cosine)*")]


def table_s3(d):
    c = P.J(d, "celleval.json")["shared"]
    rnd = ["random_0", "random_1", "random_2"]; gd = ["guide_0real_1pred_0", "guide_1real_0pred_0"]
    def block(kind):
        rows = ["| Metric | All: random | All: guide | Unexpressed: random | Unexpressed: guide | Strong: random | Strong: guide |",
                "|------|---|---|---|---|---|---|"]
        for m, name in METRICS:
            cells = []
            for sub in ("all", "inert", "strong"):
                for keys in (rnd, gd):
                    cell = rng_([c[k][sub][kind].get(m, float("nan")) for k in keys])
                    raw = [c[k][sub]["raw_mean"].get(m) for k in keys]
                    raw = [v for v in raw if v is not None and not np.isnan(v)]
                    if kind == "ceiling" and cell == "–" and raw and max(raw) <= 0:
                        cell = "n.d. (r = %s)" % rng_(raw)          # split-half correlation <= 0: no ceiling
                    cells.append(cell)
            if all(x == "–" for x in cells):
                continue
            rows.append("| %s | %s |" % (name, " | ".join(cells)))
        return "\n".join(rows)
    n = {sub: c["random_0"][sub]["n"] for sub in ("all", "inert", "strong")}
    return block("ceiling"), block("raw_mean"), n


def table_s5(d):
    g = P.J(d, "guide_donor.json")
    rows = ["| Design | Unexpressed positive, resting / stimulated | Unexpressed mean signal, resting / stimulated | "
            "Expressed mean signal, resting / stimulated (pairings) | Context-specific share | Explained share, k = 20 (shuffled) |",
            "|----|---|---|-----|--|---|"]
    for key, name in (("same_guide", "same guide, donors 1–2 against 3–4"), ("diff_guide", "different guides, donors 1–2 against 3–4")):
        v = g[key]; pr = g[key + "_pairings"]
        rows.append("| %s | %.1f%% / %.1f%% | %s / %s | %s / %s (%s; %s) | %s | %s (%s) |" % (
            name, 100 * v["inert_frac_pos_rest"], 100 * v["inert_frac_pos_stim"], f2(v["inert_mean_rest"]), f2(v["inert_mean_stim"]),
            f2(v["expressed_mean_rest"]), f2(v["expressed_mean_stim"]),
            " / ".join(f2(p["expressed_mean_rest"]) for p in pr), " / ".join(f2(p["expressed_mean_stim"]) for p in pr),
            f2(v["gamma_share"]), f2(v["readout_k20"], 4), f2(v["readout_k20_shuffled"], 4)))
    meta = "%s targets (%s unexpressed, %s expressed), %s genes" % (format(g["n_targets"], ","), format(g["n_inert"], ","), format(g["n_expressed"], ","), format(g["n_genes"], ","))
    bt = g.get("bootstrap_over_targets", {})
    ci_ = lambda k: " (95%% interval over expressed targets %.2f–%.2f)" % tuple(bt[k]) if k in bt else ""
    ratio = "%.2f%s in resting and %.2f%s in stimulated cells" % (g["guide_specific_over_diff_guide_rest"], ci_("ratio_rest"),
                                                                 g["guide_specific_over_diff_guide_stim"], ci_("ratio_stim"))
    return "\n".join(rows), meta, ratio


# ----------------------------------------------------------------------------------------------------------- figures
def fig_s1(d, out):
    k = P.J(d, "k562_inert_screens.json"); t = P.J(d, "k562_template.json")
    W, H = 183, 56
    fig = plt.figure(figsize=(W * P.MM, H * P.MM))
    ax = P.axmm(fig, W, H, 12, 10, 34, 36); P.letter(fig, W, H, 2, 3, "a")
    for i, (key, nm, col) in enumerate((("inert", "unexpressed\n(n = %d)" % k["inert"]["n"], GREY), ("expressed", "expressed\n(n = %d)" % k["expressed"]["n"], BLUE))):
        v = k[key]["guide_fraction"]; lo, hi = k[key]["ci95"]
        ax.bar(i, v, 0.6, color=col); ax.errorbar(i, v, yerr=[[v - lo], [hi - v]], color="k", lw=0.6, capsize=1.5)
        ax.text(i, hi + 0.03, "%.2f" % v, ha="center", fontsize=5.8)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["unexpressed\n(n = %d)" % k["inert"]["n"], "expressed\n(n = %d)" % k["expressed"]["n"]])
    ax.set_ylim(0, 1.1); ax.set_ylabel("share recurring across\nindependent K562 screens")
    P.title(ax, "K562 essential library")
    ax = P.axmm(fig, W, H, 62, 10, 50, 36); P.letter(fig, W, H, 52, 3, "b")
    x = np.arange(3); w = 0.36
    for j, (key, nm, col) in enumerate((("inert", "unexpressed", GREY), ("expressed", "expressed", BLUE))):
        e = k[key]; n = e["n"]
        vals = [e["within_K6"] / n, e["within_K8"] / n, e["cross"] / n]
        ax.bar(x + (j - 0.5) * w, vals, w, color=col, label=nm)
    ax.set_xticks(x); ax.set_xticklabels(["within day 6\n(essential)", "within day 8\n(genome-wide)", "across the\ntwo screens"])
    ax.set_ylabel("cross-half signal per target"); ax.legend(fontsize=5.6, loc="upper right"); P.title(ax, "Energy per target")
    ax = P.axmm(fig, W, H, 128, 10, 50, 36); P.letter(fig, W, H, 118, 3, "c")
    for j, (key, nm, col) in enumerate((("random", "random halves", RED), ("lane", "lane halves", "#e6b0aa"))):
        v = t[key]
        vals = [100 * v["inert"]["frac_pos"], 100 * v["inert_referenced"]["inert"]["frac_pos"]]
        ax.bar(np.arange(2) + (j - 0.5) * w, vals, w, color=col, label=nm)
        for xi, vv in zip(np.arange(2) + (j - 0.5) * w, vals):
            ax.text(xi, vv + 1.5, "%.0f" % vv, ha="center", fontsize=5.6)
    ax.axhline(50, color="k", lw=0.6, ls="--"); ax.set_ylim(0, 122)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["template\nremoved", "null\nbaseline"]); ax.set_ylabel("% unexpressed targets\nwith positive signal")
    ax.legend(fontsize=5.6, loc="upper center", ncol=2, handlelength=1); P.title(ax, "K562 genome-wide screen: reference")
    P.save(fig, out, "figS3_k562")


def fig_s2(d, out):
    dm = P.J(d, "dmso.json")["plates"]
    pl = sorted(dm, key=lambda s: int(s[5:])); x = np.array([int(p[5:]) for p in pl])
    W, H = 183, 52
    fig = plt.figure(figsize=(W * P.MM, H * P.MM))
    ax = P.axmm(fig, W, H, 12, 10, 50, 34); P.letter(fig, W, H, 2, 3, "a")
    hi = np.array([p in ("plate6", "plate14") for p in pl])
    ax.scatter(x, [dm[p]["drug_cross_median"] for p in pl], s=10, color="#bdc3c7", label="drug conditions (median)", zorder=2)
    ax.scatter(x, [dm[p]["cross_median"] for p in pl], s=10, color=np.where(hi, RED, BLUE), label="DMSO well difference (median)", zorder=3)
    ax.set_ylim(0, 48); ax.set_xticks(x); ax.set_xticklabels(x, fontsize=5); ax.set_xlabel("Tahoe plate"); ax.set_ylabel("within-well cross-half signal")
    ax.legend(fontsize=5.4, loc="upper left", handletextpad=0.2); P.title(ax, "Well difference against drug effects")
    ax = P.axmm(fig, W, H, 74, 10, 48, 34); P.letter(fig, W, H, 66, 3, "b")
    v = [dm[p]["well_per_well_over_drug_median"] for p in pl]
    ax.bar(x, v, 0.7, color=np.where(hi, RED, BLUE)); ax.axhline(np.median(v), color="k", lw=0.6, ls="--")
    ax.text(6.6, np.median(v) + 0.02, "median %.2f" % np.median(v), fontsize=5.4, ha="left")
    ax.set_xticks(x); ax.set_xticklabels(x, fontsize=5); ax.set_xlabel("Tahoe plate"); ax.set_ylabel("well effect per well /\nmedian drug signal")
    P.title(ax, "Size of the well effect")
    ax = P.axmm(fig, W, H, 134, 10, 46, 34); P.letter(fig, W, H, 126, 3, "c")
    v = [dm[p]["common_cross"] / (dm[p]["common_cross"] + dm[p]["specific_cross"]) for p in pl]
    ax.bar(x, v, 0.7, color=np.where(hi, RED, BLUE)); ax.axhline(np.median(v), color="k", lw=0.6, ls="--")
    ax.text(0.4, 0.93, "dashed: median %.2f" % np.median(v), fontsize=5.4, ha="left")
    ax.set_ylim(0, 1); ax.set_xticks(x); ax.set_xticklabels(x, fontsize=5); ax.set_xlabel("Tahoe plate"); ax.set_ylabel("share common to all lines")
    P.title(ax, "Line-shared part")
    P.save(fig, out, "figS2_dmso")


SCEN = [("none", "no unit\neffect"), ("additive_separate", "additive,\nseparate"), ("additive_shared", "additive,\nshared"),
        ("gain_separate", "gain,\nseparate"), ("gain_shared", "gain,\nshared"), ("both_shared", "both,\nshared")]


def fig_s3(d, out):
    s = P.J(d, "sim_v5.json")
    W, H = 183, 100
    fig = plt.figure(figsize=(W * P.MM, H * P.MM))
    panels = [("inert_frac_pos", "% null perturbations positive", 100, None, "Null perturbations"),
              ("ceiling_pearson_all", "Pearson ceiling, all", 1, None, "Ceiling, all perturbations"),
              ("ceiling_pearson_inert", "Pearson ceiling, null", 1, None, "Ceiling, null perturbations"),
              ("gamma_share", "context-specific share", 1, s["truth"]["gamma_share"], "Context-specific share"),
              ("readout_share", "explained share", 1, s["truth"]["readout_share"], "Shared → context-specific predictor")]
    pos = [(12, 8), (72, 8), (132, 8), (12, 58), (72, 58)]
    x = np.arange(len(SCEN)); w = 0.38
    for i, ((key, ylab, mul, truth, ttl), (px, py)) in enumerate(zip(panels, pos)):
        ax = P.axmm(fig, W, H, px, py, 48, 30); P.letter(fig, W, H, px - 10, py - 5, "abcde"[i])
        for j, (dsg, col) in enumerate((("random", RED), ("independent", BLUE))):
            vals = [mul * s["%s|%s" % (k, dsg)][key] for k, _ in SCEN]
            vals = [0 if np.isnan(v) else v for v in vals]
            ax.bar(x + (j - 0.5) * w, vals, w, color=col, label="random halves" if dsg == "random" else "independent units")
        if truth is not None:
            ax.axhline(truth, color=GREY, lw=0.6, ls="--"); ax.text(len(SCEN) - 0.5, truth, " truth", fontsize=5.2, color=GREY, va="center")
        if key == "inert_frac_pos":
            ax.axhline(50, color="k", lw=0.5, ls=":")
        ax.set_xticks(x); ax.set_xticklabels([l for _, l in SCEN], fontsize=5); ax.set_ylabel(ylab); P.title(ax, ttl)
    ax = P.axmm(fig, W, H, 132, 58, 48, 30); ax.axis("off"); P.letter(fig, W, H, 122, 53, "f")
    ax.legend(handles=[Patch(color=RED, label="random halves (one transduction)"), Patch(color=BLUE, label="independent units (two transductions)")],
              fontsize=5.6, loc="lower left", bbox_to_anchor=(0, -0.05))
    p = s["params"]
    txt = ("cell = μ(1 + g) + a + e\n"
           "a ~ N(0, σa²I),  g ~ N(0, σg²),  e ~ N(0, σe²I)\n\n"
           "G = %d genes, %d perturbations (20%% null)\n"
           "rank %d responses, %d cells per perturbation\n"
           "K = %d units, Dirichlet(%.1f) weights\n"
           "σa = %.1f, σg = %.1f, σe = %.0f, ρ = %.1f\n"
           "predictor k = %d; seeds %s" % (p["G"], p["P"], p["R"], p["NCELL"], p["K"], p["ALPHA"], p["SIG_A"], p["SIG_G"], p["SIG_E"], p["RHO"], p["KREAD"], ", ".join(map(str, p["seeds"]))))
    ax.text(0, 1, txt, fontsize=5.8, va="top", family="DejaVu Sans")
    P.save(fig, out, "figS1_model")


def table_ntc(d):
    """Table S1: non-targeting guides as null perturbations (ntc_null.json from cd4_ntc_null.py)"""
    v = P.J(d, "ntc_null.json")
    rows = ["| Donor, state | Guides | Cells per half (median) | Positive, random halves | Positive, different guides | Mean signal, random / different | TRBC var(z) (permuted) | TRBC deviation, between random halves |",
            "|---|---|---|---|---|---|---|---|"]
    for lab, k in (("Donor 1, resting", "D1_rest"), ("Donor 1, stimulated", "D1_stim"), ("Donor 2, resting", "D2_rest"), ("Donor 2, stimulated", "D2_stim")):
        r = v[k]
        rows.append("| %s | %d | %d | %.1f%% (%.1f–%.1f) | %.1f%% (%.1f–%.1f) | %s / %s | %.1f (%.2f) | %.2f |" % (
            lab, r["n_guides"], r["cells_per_half_median"], 100 * r["random_frac_pos"], 100 * r["random_frac_pos_ci95"][0],
            100 * r["random_frac_pos_ci95"][1], 100 * r["pairs_frac_pos_mean"], 100 * r["pairs_frac_pos_range"][0],
            100 * r["pairs_frac_pos_range"][1], f2(r["random_mean_cross"], 1), f2(r["pairs_mean_cross"], 2), r["trbc_var_z"],
            r["trbc_var_z_permuted"], r["trbc_random_half_pearson"]))
    return "\n".join(rows)


def table_ntc_cl(d):
    """Table S1, second block: non-targeting constructs in cell-line screens (cl_ntc_null.py) and the CD4 ratio
    (cd4_ntc_ratio.py)"""
    cl = P.J(d, "cl_ntc_replogle.json"); cl.update(P.J(d, "cl_ntc_nadig.json")); cl.update(P.J(d, "cl_ntc_h1.json")); rat = P.J(d, "cd4_ntc_ratio.json")
    rows = ["| Screen | Constructs | Cells per half (median) | Positive, random halves | Positive, different constructs | Mean signal, non-targeting / targeting | Ratio |",
            "|---|---|---|---|---|---|---|"]
    for lab, k in (("CD4 donor 1, resting", "D1_rest"), ("CD4 donor 2, resting", "D2_rest")):
        r = rat[k]
        rows.append("| %s | %d | – | see above | see above | %s / %s | %.2f |" % (lab, r["n_constructs"], f2(r["random_mean_cross"], 1),
                    f2(r["target_random_mean_cross"], 1), r["ratio_nt_to_target_mean"]))
    for lab, k in (("K562 essential (day 6)", "k562_ess"), ("K562 genome-wide (day 8, part 1)", "k562_gw1"), ("RPE1 essential", "rpe1"),
                   ("Jurkat essential", "jurkat"), ("HepG2 essential", "hepg2"), ("H1 hESC, VCC 2025 training data", "h1")):
        r = cl[k]
        rows.append("| %s | %d | %d | %.1f%% (%.1f–%.1f) | %.1f%% (%.1f–%.1f) | %s / %s | %.3f |" % (
            lab, r["n_constructs"], r["cells_per_half_median"], 100 * r["random_frac_pos"], 100 * r["random_frac_pos_ci95"][0],
            100 * r["random_frac_pos_ci95"][1], 100 * r["pairs_frac_pos_mean"], 100 * r["pairs_frac_pos_range"][0],
            100 * r["pairs_frac_pos_range"][1], f2(r["random_mean_cross"], 2), f2(r["target_random_mean_cross"], 1), r["ratio_nt_to_target_mean"]))
    return "\n".join(rows)


def main():
    d, tpl, out_md, fig_dir = sys.argv[1:5]
    os.makedirs(fig_dir, exist_ok=True)
    ntc = table_ntc(d); ntc_cl = table_ntc_cl(d); s1, shuf1 = table_s1(d); s2a, s2b = table_s2(d); s3a, s3b, n3 = table_s3(d); s5, s5meta, s5ratio = table_s5(d)
    text = open(tpl, encoding="utf-8").read()
    for k, v in (("{TABLE_NTC}", ntc), ("{TABLE_NTC_CL}", ntc_cl), ("{TABLE_S1}", s1), ("{S1_SHUF}", "%.4f" % shuf1), ("{TABLE_S2_CD4}", s2a), ("{TABLE_S2_TAHOE}", s2b),
                 ("{TABLE_S3_CORR}", s3a), ("{TABLE_S3_RAW}", s3b), ("{TABLE_S5}", s5), ("{S5_META}", s5meta), ("{S5_RATIO}", s5ratio)):
        assert k in text, k
        text = text.replace(k, v)
    open(out_md, "w", encoding="utf-8").write(text)
    fig_s1(d, fig_dir); fig_s2(d, fig_dir); fig_s3(d, fig_dir)
    print("ok", n3)


if __name__ == "__main__":
    main()
