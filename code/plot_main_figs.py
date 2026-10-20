"""Main data figures of the clone-effect preprint. All numbers come from the small result files in <data_dir>
(data/ in this repository; data/README.md lists the script behind each file); per-perturbation arrays from
export_fig_data.py.
fig2_null.{pdf,png}       null perturbations look reproducible under random halves (CD4 unexpressed targets, donor
                          halves, non-targeting guides across screens)
fig3_clones.{pdf,png}     direct clone evidence (TRBC marker), depth, composition of the random-half signal
                          (guide_donor.json from cd4_guide_donor.py)
fig4_ceilings.{pdf,png}   consequences for ceilings and reliability labels (cell-eval itself, Wang-style labels,
                          Spearman-Brown at chance)
fig5_contexts.{pdf,png}   cross-context structure (shared vs separate units, absolute energy, predictors, model)
Panel titles are short labels; the message of each figure is in its legend title. Axis labels in sentence case with
units in parentheses; error bars wherever an interval exists (defined in the legends).
Usage: plot_main_figs.py <data_dir> <out_dir>"""
import sys, os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import plot_schematics as PS

plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"], "font.size": 6.5, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.spines.top": False,
                     "axes.spines.right": False, "pdf.fonttype": 42, "legend.frameon": False, "axes.titlesize": 7,
                     "axes.labelsize": 6.5, "xtick.labelsize": 6, "ytick.labelsize": 6})
RED, BLUE, GREY, ORANGE = PS.RED, PS.BLUE, PS.GREY, PS.ORANGE
MM = 1 / 25.4


def J(d, f):
    return json.load(open(os.path.join(d, f)))


def axmm(fig, W, H, x, y, w, h):
    """axes at (x, y) mm from top-left, size w x h mm"""
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])


def letter(fig, W, H, x, y, s):
    fig.text(x / W, 1 - y / H, s, fontsize=8, fontweight="bold", va="top", ha="left")


def title(ax, s, sub=None):
    """short black panel label; the message of a figure is in its legend title (journal style)"""
    ax.set_title(s, loc="left", fontsize=6.5, pad=3)


MIN_PT = 5.6   # smallest type at native size (183 mm); printed at 170 mm this stays above the 5-pt minimum


def floor_type(fig, minimum=MIN_PT):
    import matplotlib.text as mtext
    for t in fig.findobj(mtext.Text):
        if t.get_text().strip() and t.get_fontsize() < minimum:
            t.set_fontsize(minimum)


def sentence_case(fig):
    """axis labels start with a capital letter (journal style); tick and key labels are left as written"""
    for ax in fig.axes:
        for get, put in ((ax.get_xlabel, ax.set_xlabel), (ax.get_ylabel, ax.set_ylabel)):
            s = get()
            if s and s[0].islower():
                put(s[0].upper() + s[1:])


def share_ci(p, n, n_boot=1000, seed=0):
    """bootstrap 95% interval of a share over n items; resampling a 0/1 vector is a binomial draw, so this is exact"""
    b = np.random.default_rng(seed).binomial(int(n), p, n_boot) / n
    return np.percentile(b, [2.5, 97.5])


def err(ax, x, v, lo, hi, scale=1.0):
    ax.errorbar(x, scale * v, yerr=[[scale * (v - lo)], [scale * (hi - v)]], color="k", lw=0.6, capsize=1.5, ls="none")


def save(fig, out, name):
    sentence_case(fig)
    floor_type(fig)
    fig.savefig(os.path.join(out, name + ".pdf")); fig.savefig(os.path.join(out, name + ".png"), dpi=300); plt.close(fig)


def edge_ticks(ax, lo, hi, ticks):
    """clipped histogram: label the outer bins as <= lo and >= hi"""
    lab = []
    for t in ticks:
        s = ("%g" % t).replace("-", "−")
        lab.append("≤" + s if t == lo else ("≥" + s if t == hi else s))
    ax.set_xticks(ticks); ax.set_xticklabels(lab)


def patches(items):
    return [Patch(color=c, label=l) for c, l in items]


# ---------------------------------------------------------------------------------------------------- figure 2
def fig2(d, out):
    D1, D2 = np.load(os.path.join(d, "cd4_D1.npz")), np.load(os.path.join(d, "cd4_D2.npz"))
    W, H = 183, 112
    fig = plt.figure(figsize=(W * MM, H * MM))
    # a: distribution of cross-half signal, unexpressed targets, donor 1 rest
    ax = axmm(fig, W, H, 12, 14, 70, 33); letter(fig, W, H, 2, 3, "a")
    m = D1["inert"]
    for k, col, nm in (("random", RED, "random halves"), ("guide", BLUE, "guide halves")):
        v = D1["cross_%s_rest" % k][m]
        ax.hist(np.clip(v, -15, 30), bins=np.linspace(-15, 30, 61), color=col, alpha=0.75, histtype="stepfilled", lw=0,
                label="%s: %d%% above 0" % (nm, round(100 * (v > 0).mean())))
    ax.axvline(0, color="k", lw=0.6, ls="--")
    edge_ticks(ax, -15, 30, [-15, -10, 0, 10, 20, 30])
    ax.set_xlabel("cross-half signal (half A · half B)"); ax.set_ylabel("targets")
    title(ax, "CD4 T cells, unexpressed targets")
    ax.set_ylim(0, ax.get_ylim()[1] * 1.3); ax.legend(fontsize=5.6, loc="upper right", handlelength=1.0)
    # b: share positive across donors / states / designs
    ax = axmm(fig, W, H, 102, 14, 78, 33); letter(fig, W, H, 92, 3, "b")
    ds = J(d, "donor_split.json")["C_null_baseline_donor_split"]["inert"]
    groups = [("D1\nrest", D1, "rest"), ("D1\nstim", D1, "stim"), ("D2\nrest", D2, "rest"), ("D2\nstim", D2, "stim")]
    x = np.arange(4); w = 0.36
    for j, (k, col) in enumerate((("random", RED), ("guide", BLUE))):
        for xi, (_, dd, s) in zip(x, groups):
            ind = dd["cross_%s_%s" % (k, s)][dd["inert"]] > 0
            ax.bar(xi + (j - 0.5) * w, 100 * ind.mean(), w, color=col); err(ax, xi + (j - 0.5) * w, ind.mean(), *share_ci(ind.mean(), ind.size), scale=100)
    for xi, s in ((4.15, "rest"), (4.85, "stim")):   # donor halves: 1,479 unexpressed targets (donors 1-2 against 3-4)
        p = ds["%s_frac_pos" % s]; ax.bar(xi, 100 * p, w, color=BLUE); err(ax, xi, p, *share_ci(p, 1479), scale=100)
    ax.axhline(50, color="k", lw=0.6, ls="--")
    ax.set_ylim(0, 118); ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xticks(list(x) + [4.15, 4.85]); ax.set_xticklabels([g for g, _, _ in groups] + ["donors\nrest", "donors\nstim"], fontsize=5.6)
    ax.set_ylabel("Targets with positive\ncross-half signal (%)")
    ax.legend(handles=patches(((RED, "random halves"), (BLUE, "guide or donor halves"))), fontsize=5.6, loc="upper center",
              ncol=2, handlelength=1, columnspacing=1)
    title(ax, "Unexpressed targets")
    # c: non-targeting guides across screens (findings cc): positive share and size relative to targeting perturbations
    letter(fig, W, H, 2, 62, "c")
    ntc, rat = J(d, "ntc_null.json"), J(d, "cd4_ntc_ratio.json")
    cl = dict(J(d, "cl_ntc_replogle.json")); cl.update(J(d, "cl_ntc_nadig.json")); cl.update(J(d, "cl_ntc_h1.json"))
    rows = [("CD4\nD1", ntc["D1_rest"], rat["D1_rest"]["ratio_nt_to_target_mean"]), ("CD4\nD2", ntc["D2_rest"], rat["D2_rest"]["ratio_nt_to_target_mean"])]
    for k, lab in (("k562_ess", "K562\nessential"), ("k562_gw1", "K562\nGW"), ("rpe1", "RPE1"), ("jurkat", "Jurkat"), ("hepg2", "HepG2"), ("h1", "H1\nVCC 2025")):
        rows.append((lab, cl[k], cl[k]["ratio_nt_to_target_mean"]))
    x = np.arange(len(rows)); w = 0.36
    axl = axmm(fig, W, H, 13, 72, 78, 30)
    axl.bar(x - w / 2, [100 * r[1]["random_frac_pos"] for r in rows], w, color=RED)
    axl.bar(x + w / 2, [100 * r[1]["pairs_frac_pos_mean"] for r in rows], w, color=BLUE)
    for xi, r in zip(x, rows):   # random halves: bootstrap 95% interval over guides; different guides: range over the 50 pairings
        err(axl, xi - w / 2, r[1]["random_frac_pos"], *r[1]["random_frac_pos_ci95"], scale=100)
        err(axl, xi + w / 2, r[1]["pairs_frac_pos_mean"], *r[1]["pairs_frac_pos_range"], scale=100)
    axl.axhline(50, color="k", lw=0.6, ls="--"); axl.set_ylim(0, 118); axl.set_yticks([0, 25, 50, 75, 100])
    axl.set_xticks(x); axl.set_xticklabels([r[0] for r in rows], fontsize=5.4); axl.set_ylabel("Guides with positive\ncross-half signal (%)")
    axl.legend(handles=patches(((RED, "random halves"), (BLUE, "different guides"))), fontsize=5.6, loc="upper right", ncol=2, handlelength=1, columnspacing=1)
    title(axl, "Non-targeting guides")
    axr = axmm(fig, W, H, 113, 72, 65, 30)
    axr.bar(x, [r[2] for r in rows], 0.6, color=RED)
    axr.set_yscale("log"); axr.set_ylim(0.003, 3); axr.axhline(1, color="k", lw=0.6, ls="--")
    axr.set_yticks([0.01, 0.1, 1]); axr.set_yticklabels(["0.01", "0.1", "1"]); axr.minorticks_off()
    for xi, r in zip(x, rows):
        axr.text(xi, r[2] * 1.15, ("%.3f" % r[2]) if r[2] < 0.01 else (("%.2f" % r[2]) if r[2] < 1 else ("%.1f" % r[2])), ha="center", va="bottom", fontsize=5.2)
    axr.set_xticks(x); axr.set_xticklabels([r[0] for r in rows], fontsize=5.4)
    axr.set_ylabel("non-targeting / targeting\ncross-half signal")
    title(axr, "Relative to targeting perturbations")
    save(fig, out, "fig2_null")


# ---------------------------------------------------------------------------------------------------- figure 3
def fig3(d, out):
    D1, D2 = np.load(os.path.join(d, "cd4_D1.npz")), np.load(os.path.join(d, "cd4_D2.npz"))
    c1, c2 = J(d, "clone_D1.json"), J(d, "clone_D2.json")
    W, H = 183, 167
    fig = plt.figure(figsize=(W * MM, H * MM))
    ax = axmm(fig, W, H, 0, 0, 89, 52); ax.set_xlim(0, 89); ax.set_ylim(0, 52); ax.axis("off"); PS.draw_fig3a(ax)
    # b: overdispersion
    ax = axmm(fig, W, H, 104, 14, 74, 33); letter(fig, W, H, 95, 3, "b")
    z, zp = D1["trbc_z_rest"], D1["trbc_zperm_rest"]; ok = D1["trbc_n_rest"] >= 30
    bins = np.linspace(-25, 25, 81)
    ax.hist(np.clip(zp[ok], -25, 25), bins=bins, color=GREY, alpha=0.8, label="perturbation labels permuted", histtype="stepfilled", lw=0)
    ax.hist(np.clip(z[ok], -25, 25), bins=bins, color=ORANGE, alpha=0.75, label="observed", histtype="stepfilled", lw=0)
    ax.set_yscale("log"); edge_ticks(ax, -25, 25, [-25, -10, 0, 10, 25])
    ax.set_yticks([1, 10, 100, 1000]); ax.set_yticklabels(["1", "10", "100", "1,000"]); ax.minorticks_off()   # no superscripts (type floor)
    ax.set_xlabel("TRBC1 share of a perturbation (z-score)"); ax.set_ylabel("perturbations")
    ax.legend(fontsize=5.6, loc="upper left")
    title(ax, "TRBC1 share per perturbation")
    # c: scatter random vs guide halves
    letter(fig, W, H, 2, 58, "c")
    for j, (k, col, nm) in enumerate((("random", RED, "random halves"), ("guide", BLUE, "guide halves"))):
        ax = axmm(fig, W, H, 12 + j * 50, 70, 36, 32)
        a, b = D1["trbc_dA_%s_rest" % k], D1["trbc_dB_%s_rest" % k]
        okk = (D1["trbc_nA_%s_rest" % k] >= 15) & (D1["trbc_nB_%s_rest" % k] >= 15)
        ax.scatter(a[okk], b[okk], s=1.2, c=col, alpha=0.35, lw=0, rasterized=True)
        lim = 0.45; ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.axhline(0, color="k", lw=0.4); ax.axvline(0, color="k", lw=0.4)
        r = np.corrcoef(a[okk], b[okk])[0, 1]
        ax.set_xlabel("TRBC1 share deviation, half A"); ax.set_ylabel("TRBC1 share deviation, half B" if j == 0 else "")
        ax.set_title("%s (r = %s)" % (nm[0].upper() + nm[1:], ("%.2f" % r).replace("-", "−")), loc="left", fontsize=6.5, pad=3)
    # d: link to expression on unexpressed targets
    ax = axmm(fig, W, H, 122, 70, 56, 32); letter(fig, W, H, 112, 58, "d")
    for j, (dd, nm, mk, ls) in enumerate(((D1, "donor 1", "o", "-"), (D2, "donor 2", "s", "--"))):
        m = dd["inert"] & (dd["trbc_n_rest"] >= 30)
        sk = np.abs(dd["trbc_z_rest"][m]); ec = dd["cross_random_rest"][m]
        q = np.quantile(sk, [0.25, 0.5, 0.75]); bins_ = np.digitize(sk, q)
        mu = [ec[bins_ == b].mean() for b in range(4)]; se = [ec[bins_ == b].std() / np.sqrt((bins_ == b).sum()) for b in range(4)]
        ax.errorbar(np.arange(4) + (j - 0.5) * 0.12, mu, yerr=se, marker=mk, ms=3.2, color=RED, lw=0.8, ls=ls, capsize=1.5, label=nm)
    ax.set_xticks(range(4)); ax.set_xticklabels(["lowest", "2nd", "3rd", "highest"]); ax.set_xlabel("clonal skew (|z|), quartile")
    ax.set_ylabel("cross-half signal,\nrandom halves"); ax.legend(fontsize=5.6, loc="upper left")
    title(ax, "Unexpressed targets")
    # e: depth (resting solid, stimulated dashed; donor 1 circles, donor 2 squares)
    letter(fig, W, H, 2, 112, "e")
    axl = axmm(fig, W, H, 13, 124, 40, 28); axr = axmm(fig, W, H, 68, 124, 40, 28)
    for dd, mk in ((c1, "o"), (c2, "s")):
        for st, ls in (("rest", "-"), ("stim48", "--")):
            for des, col in (("random", RED), ("guide", BLUE)):
                q, h = dd[st]["depth"]["%s_quarter" % des], dd[st]["depth"]["%s_half" % des]
                xq = [q["median_cells_per_side"][0], h["median_cells_per_side"][0]]
                axl.plot(xq, [q["inert"]["cross_per_target"], h["inert"]["cross_per_target"]], ls=ls, marker=mk, color=col, ms=2.6, lw=0.7)
                axr.plot(xq, [q["inert"]["pearson_mean"], h["inert"]["pearson_mean"]], ls=ls, marker=mk, color=col, ms=2.6, lw=0.7)
    for a in (axl, axr):
        a.axhline(0, color="k", lw=0.5); a.set_xlabel("cells per half (median)"); a.set_xlim(38, 125)
    axl.set_ylabel("cross-half signal\nper target"); axr.set_ylabel("agreement between\nhalves (Pearson r)")
    title(axl, "Unexpressed targets")
    title(axr, "Unexpressed targets")
    from matplotlib.lines import Line2D
    key = [Line2D([], [], color=RED, lw=0.8, label="random halves"), Line2D([], [], color=BLUE, lw=0.8, label="guide halves"),
           Line2D([], [], color="k", marker="o", ls="none", ms=2.6, label="donor 1"), Line2D([], [], color="k", marker="s", ls="none", ms=2.6, label="donor 2"),
           Line2D([], [], color="k", lw=0.7, label="resting"), Line2D([], [], color="k", lw=0.7, ls="--", label="stimulated")]
    fig.legend(handles=key, loc="upper left", bbox_to_anchor=(11 / W, 1 - 160.5 / H), ncol=6, fontsize=5.4, handlelength=1.6,
               columnspacing=1.2, frameon=False)
    # f: how much of the random-half signal is shared only within one transduction (expressed targets)
    ax = axmm(fig, W, H, 128, 124, 50, 28); letter(fig, W, H, 118, 112, "f")
    gd = J(d, "guide_donor.json")
    bars = []
    for lab, dd, st, rk in (("CD4 D1\nrest", c1, "rest", "rest"), ("CD4 D1\nstim", c1, "stim48", "stim"), ("CD4 D2\nrest", c2, "rest", "rest")):
        R = dd[st]["depth"]["random_half"]["expressed"]["cross_per_target"]; G = dd[st]["depth"]["guide_half"]["expressed"]["cross_per_target"]
        rep = G / R; gsp = gd["guide_specific_over_diff_guide_" + rk] * rep
        bars.append((lab, 1 - rep - gsp, gsp, rep))
    x = np.arange(len(bars)); w = 0.6
    for xi, (lab, tr, gsp, rep) in enumerate(bars):
        ax.bar(xi, 100 * tr, w, color=RED); ax.bar(xi, 100 * gsp, w, bottom=100 * tr, color="#e6b0aa")
        ax.bar(xi, 100 * rep, w, bottom=100 * (tr + gsp), color=BLUE)
        ax.text(xi, 100 * tr / 2, "%d%%" % round(100 * tr), ha="center", va="center", fontsize=5.6, color="white")
    ax.set_xticks(x); ax.set_xticklabels([b[0] for b in bars], fontsize=5.4)
    ax.set_ylim(0, 155); ax.set_yticks([0, 25, 50, 75, 100]); ax.set_ylabel("Share of random-half signal (%)")
    ax.legend(handles=patches(((RED, "shared within one transduction"), ("#e6b0aa", "specific to the guide"), (BLUE, "reproduced by the other guide"))),
              fontsize=4.8, loc="upper left", ncol=1, handlelength=0.9, borderaxespad=0.1, labelspacing=0.2)
    title(ax, "Expressed targets")
    save(fig, out, "fig3_clones")


# ---------------------------------------------------------------------------------------------------- figure 4
def fig4(d, out):
    ce = J(d, "celleval.json")["shared"]
    W, H = 183, 127
    fig = plt.figure(figsize=(W * MM, H * MM))
    mets = [("overlap_at_N", "DE overlap"), ("precision_at_N", "DE precision"), ("de_sig_genes_recall", "DE recall"),
            ("de_direction_match", "DE direction match*"), ("pearson_delta", "Pearson Δ"), ("discrimination_score_l1", "discrimination (L1)*"),
            ("discrimination_score_cosine", "discrimination (cosine)*")]
    letter(fig, W, H, 2, 3, "a")
    for j, (sub, nm) in enumerate((("all", "all 300 perturbations"), ("inert", "150 unexpressed targets"), ("strong", "150 strength-matched knockdowns"))):
        ax = axmm(fig, W, H, 36 + j * 49, 11, 44, 44)
        for i, (m, lab) in enumerate(mets):
            r = np.nanmean([ce["random_%d" % s][sub]["ceiling"].get(m) or np.nan for s in (0, 1, 2)])
            gv = [x for x in (ce[k][sub]["ceiling"].get(m) for k in ("guide_0real_1pred_0", "guide_1real_0pred_0")) if x is not None and not np.isnan(x)]
            g = float(np.mean(gv)) if gv else np.nan
            y = len(mets) - 1 - i
            if np.isnan(g):
                # split-half correlation <= 0: cell-eval returns no ceiling; drawn as an open marker at 0 and labelled
                ax.plot([0, r], [y, y], color="#d5dbdb", lw=1.6, zorder=1)
                ax.scatter([0], [y], s=14, facecolor="white", edgecolor=BLUE, zorder=3, lw=0.8)
                ax.text(0.03, y + 0.28, "r ≤ 0", fontsize=5.4, va="bottom")
            else:
                ax.plot([g, r], [y, y], color="#d5dbdb", lw=1.6, zorder=1)
                ax.scatter([g], [y], s=14, color=BLUE, zorder=3, lw=0)
            ax.scatter([r], [y], s=14, color=RED, zorder=3, lw=0)
        ax.axvline(2 * 0.5 / 1.5, color=GREY, lw=0.6, ls=":")
        ax.set_xlim(-0.03, 1.02); ax.set_ylim(-0.6, len(mets) - 0.4); ax.set_yticks(range(len(mets)))
        ax.set_yticklabels([l for _, l in mets][::-1] if j == 0 else [])
        ax.set_xlabel("noise ceiling"); ax.set_title(nm, loc="left", fontsize=6.5, pad=3)
    fig.text(140 / W, 1 - 63.6 / H, "●", color=RED, fontsize=7, va="top"); fig.text(143 / W, 1 - 64 / H, "random halves", fontsize=6, va="top")
    fig.text(160 / W, 1 - 63.6 / H, "●", color=BLUE, fontsize=7, va="top"); fig.text(163 / W, 1 - 64 / H, "guide halves", fontsize=6, va="top")
    # b: Wang-style reliability labels, both donors (the Spearman-Brown curve is Supplementary Fig. S3)
    ax = axmm(fig, W, H, 12, 85, 120, 30); letter(fig, W, H, 2, 73, "b")
    w = 0.36; cats = ("all", "inert", "expressed")
    for g, (fname, nn, x0) in enumerate((("ci_cd4.json", {"all": 6615, "inert": 1235, "expressed": 1097}, 0.0),
                                          ("ci_cd4_D2.json", {"all": 5003, "inert": 869, "expressed": 859}, 3.6))):
        ci = J(d, fname)["wang_reliable_share_rest"]; x = x0 + np.arange(3)
        for j, (k, col, nm) in enumerate((("random", RED, "random halves"), ("guide", BLUE, "guide halves"))):
            ax.bar(x + (j - 0.5) * w, [100 * ci[k][s] for s in cats], w, color=col, label=nm if g == 0 else None)
            for xi, s in zip(x, cats):
                err(ax, xi + (j - 0.5) * w, ci[k][s], *share_ci(ci[k][s], nn[s]), scale=100)
        ax.text(x0 + 1, -15.5, "donor %d" % (g + 1), ha="center", va="top", fontsize=6.5)
    ax.set_xticks([0, 1, 2, 3.6, 4.6, 5.6]); ax.set_xticklabels(["all", "unexpressed", "expressed"] * 2); ax.set_ylabel("Labelled reliable (%)")
    ax.set_ylim(0, 62); ax.legend(fontsize=5.6, loc="upper left")
    title(ax, "Reliability labels")
    save(fig, out, "fig4_ceilings")


# ---------------------------------------------------------------------------------------------------- figure 5
def fig5(d, out):
    W, H = 183, 160
    fig = plt.figure(figsize=(W * MM, H * MM))
    ax = axmm(fig, W, H, 0, 0, 183, 50); ax.set_xlim(0, 183); ax.set_ylim(0, 50); ax.axis("off"); PS.draw_fig5a(ax)
    ci = J(d, "ci_cd4.json")
    # b: the same clones in both states: TRBC1 deviation (mean f - f0) / sd_f = z / sqrt(n), resting against stimulated
    ax = axmm(fig, W, H, 134, 8, 42, 34); letter(fig, W, H, 124, 1, "b")
    D1 = np.load(os.path.join(d, "cd4_D1.npz"))
    nr, ns = D1["trbc_n_rest"], D1["trbc_n_stim"]
    ok = (nr >= 30) & (ns >= 30) & np.isfinite(D1["trbc_z_rest"]) & np.isfinite(D1["trbc_z_stim"])
    xr, ys = D1["trbc_z_rest"][ok] / np.sqrt(nr[ok]), D1["trbc_z_stim"][ok] / np.sqrt(ns[ok])
    ax.scatter(xr, ys, s=1.2, c=RED, alpha=0.35, lw=0, rasterized=True)
    lim = 0.85; ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.axhline(0, color="k", lw=0.4); ax.axvline(0, color="k", lw=0.4)
    ax.set_xlabel("TRBC1 deviation, resting"); ax.set_ylabel("TRBC1 deviation, stimulated")
    ax.set_title("One clone, two states (r = %.2f)" % np.corrcoef(xr, ys)[0, 1], loc="left", fontsize=6.5, pad=3)
    # c: absolute energy vs share
    ax = axmm(fig, W, H, 12, 67, 44, 31); letter(fig, W, H, 2, 55, "c")
    x = np.arange(2); w = 0.36
    for j, (k, col, nm) in enumerate((("random_r1", RED, "random halves"), ("guide", BLUE, "guide halves"))):
        ax.bar(x + (j - 0.5) * w, [ci[k]["beta_cross_per_target"], ci[k]["gamma_cross_per_target"]], w, color=col, label=nm)
    ax.set_xticks(x); ax.set_xticklabels(["shared", "context-\nspecific"]); ax.set_ylabel("reproducible energy per target")
    for j, (k, col) in enumerate((("random_r1", RED), ("guide", BLUE))):
        g = ci[k]["gamma_share"]; ax.text(1 + (j - 0.5) * w, ci[k]["gamma_cross_per_target"] + 0.25, "%.2f" % g[0], ha="center", fontsize=5.5)
    ax.legend(fontsize=5.6)
    title(ax, "CD4 T cells, donor 1")
    # d: CD4 predictor and K562
    ax = axmm(fig, W, H, 76, 67, 104, 31); letter(fig, W, H, 66, 55, "d")
    kd1, kd2, ci2, rp = J(d, "kd_D1.json"), J(d, "kd_D2.json"), J(d, "ci_cd4_D2.json"), J(d, "replogle_indep.json")
    bars = [("CD4 D1", ci["random_r1"]["readout_k20"][0], ci["guide"]["readout_k20"][0]),
            ("CD4 D1\nmatched", kd1["kd_matched|cell"]["readout_k20"], kd1["kd_matched|guide"]["readout_k20"]),
            ("CD4 D2", ci2["random_r1"]["readout_k20"][0], ci2["guide"]["readout_k20"][0]),
            ("CD4 D2\nmatched", kd2["kd_matched|cell"]["readout_k20"], kd2["kd_matched|guide"]["readout_k20"]),
            ("K562–RPE1\n(separate)", rp["target_K8"]["pair_same_k20"], rp["target_K8"]["pair_indepK_k20"])]
    x = np.arange(len(bars)); w = 0.36
    ax.bar(x - w / 2, [b[1] for b in bars], w, color=RED, label="same cells or screen")
    ax.bar(x + w / 2, [b[2] for b in bars], w, color=BLUE, label="independent guides or screen")
    # bootstrap 95% intervals over targets, available for all CD4 targets (ci_cd4*.json); not computed for the subsets
    for xi, cc in ((0, ci), (2, ci2)):
        for dx, k in ((-w / 2, "random_r1"), (w / 2, "guide")):
            err(ax, xi + dx, *cc[k]["readout_k20"])
    for xi, b in zip(x, bars):
        if abs(b[2]) < 0.01:
            ax.text(xi + w / 2, 0.04, "≈0", ha="center", fontsize=5.5)
    ax.set_xticks(x); ax.set_xticklabels([b[0] for b in bars], fontsize=5.5); ax.set_ylabel("held-out explained share"); ax.set_ylim(0, 0.9)
    ax.legend(fontsize=5.4, loc="upper center", ncol=2, handlelength=1)
    title(ax, "Shared-to-specific predictor")
    # e: model
    sim = J(d, "sim_v6.json")
    letter(fig, W, H, 2, 110, "e")
    ax = axmm(fig, W, H, 12, 122, 62, 30)
    sc = [("fab_none", "no unit effect"), ("fab_low_separate", "separate units"), ("fab_low_shared", "shared units")]
    x = np.arange(3); w = 0.36
    ax.bar(x - w / 2, [sim["%s|random" % k]["readout_share"] for k, _ in sc], w, color=RED, label="random halves")
    ax.bar(x + w / 2, [sim["%s|independent" % k]["readout_share"] for k, _ in sc], w, color=BLUE, label="independent units")
    for xi, (k, _) in zip(x, sc):
        for dx, dsg, col in ((-w / 2, "random", RED), (w / 2, "independent", BLUE)):
            v = sim["%s|%s" % (k, dsg)]["readout_share"]
            ax.text(xi + dx, max(v, 0) + 0.003, "%.2f" % abs(v) if abs(v) >= 0.005 else "0", ha="center", fontsize=5.2)
    ax.axhline(0, color="k", lw=0.5); ax.set_xticks(x); ax.set_xticklabels([l for _, l in sc]); ax.set_ylabel("explained share")
    ax.set_ylim(0, 0.125); ax.legend(fontsize=5.6, loc="upper left")
    title(ax, "Model, no true structure")
    ax = axmm(fig, W, H, 104, 122, 74, 30)
    for kind, mk, lab in (("iso", "o", "unit effects in random directions"), ("low", "s", "unit effects along the response")):
        ratio, ri = [], []
        for sa in ("0.2", "0.6", "1.2", "1.5", "2.0"):
            a = sim["atten_sa%s_%s|random" % (sa, kind)]; b = sim["atten_sa%s_%s|independent" % (sa, kind)]
            ratio.append(a["inert_cross_per_target"] / b["expressed_cross_per_target"]); ri.append(b["readout_share"])
        ax.plot(ratio, ri, "-" + mk, color=BLUE, ms=3.2, lw=0.8, label=lab, alpha=1 if kind == "low" else 0.45)
    ax.axvspan(4.0, 7.5, color="#eaecee", zorder=0); ax.text(5.75, 0.93, "CD4 donors", fontsize=5.5, ha="center")
    ax.axhline(0.70, color=GREY, lw=0.6, ls="--"); ax.text(2.0, 0.72, "true value", fontsize=5.5, ha="left")
    ax.set_xlabel("unit-to-signal energy ratio"); ax.set_ylabel("explained share,\nindependent units"); ax.set_ylim(0, 1.0)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8]); ax.legend(fontsize=5.4, loc="upper right", bbox_to_anchor=(0.99, 0.5))
    title(ax, "Model, true explained share 0.70")
    save(fig, out, "fig5_contexts")


if __name__ == "__main__":
    d, out = sys.argv[1:3]
    fig2(d, out); fig3(d, out); fig4(d, out); fig5(d, out); print("ok")
