"""Figures of the next version of the clone preprint that differ from v1 (local). Style and helpers from
plot_main_figs.py; v1 per-perturbation arrays and result files from <v1_data> (the public repository's data/ folder);
new results from <results> (cross_context/model_eval/results).
fig2_next   null perturbations look reproducible under random halves (v1 Fig. 2a-c without the ratio panel; d = depth,
            v1 Fig. 3e)
fig3_next   the shared state is clonal (v1 Fig. 3a-c; d = TRBC halves against random halves, donors 1-2, both states)
fig4_next   clone halves and the effective number of clones (organoid clone barcodes; unified clone strength across
            datasets; CD4 composition of the random-half signal with the TRBC-halves retention)
fig7_next   transfer-model scores against same-clone and separate-transduction answers
Usage: plot_next_figs.py <v1_data> <results> <out_dir>"""
import sys, os, json, re
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_main_figs as P
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

RED, BLUE, GREY, ORANGE, MM = P.RED, P.BLUE, P.GREY, P.ORANGE, P.MM
COL = {"copy": "#8a8a8a", "linear": "#000000", "state_transition": "#4d4d4d"}
LS = {"copy": ":", "linear": "-", "state_transition": "--"}
MK = {"copy": "o", "linear": "^", "state_transition": "D"}
NAMES = {"copy": "copy resting effect", "linear": "linear (ridge)", "state_transition": "STATE"}


def J(*p):
    return json.load(open(os.path.join(*p), encoding="utf-8"))


def depth_panels(fig, W, H, c1, c2, x0, y0, w, h, gap):
    # c1, c2: depth_grid_D1/D2.json (1/8, 1/4, 3/8, 1/2 of a target's cells per side; code/a100/depth_grid_cd4.py)
    axl = P.axmm(fig, W, H, x0, y0, w, h); axr = P.axmm(fig, W, H, x0 + w + gap, y0, w, h)
    for dd, mk in ((c1, "o"), (c2, "s")):
        for st, ls in (("rest", "-"), ("stim48", "--")):
            for des, col in (("random", RED), ("guide", BLUE)):
                g = dd[st]["grid"][des]; xq = [r["median_cells_per_side"] for r in g]
                axl.plot(xq, [r["inert"]["cross_per_target"] for r in g], ls=ls, marker=mk, color=col, ms=2.4, lw=0.7)
                axr.plot(xq, [r["inert"]["pearson_mean"] for r in g], ls=ls, marker=mk, color=col, ms=2.4, lw=0.7)
    for a in (axl, axr):
        a.axhline(0, color="k", lw=0.5); a.set_xlabel("cells per half (median)"); a.set_xlim(0, 125)
    axl.set_ylabel("cross-half signal\nper target"); axr.set_ylabel("agreement between\nhalves (Pearson r)")
    P.title(axl, "Unexpressed targets"); P.title(axr, "Unexpressed targets")
    return axl, axr


def fig2(v1, res, out):
    D1, D2 = np.load(os.path.join(v1, "cd4_D1.npz")), np.load(os.path.join(v1, "cd4_D2.npz"))
    W, H = 183, 122
    fig = plt.figure(figsize=(W * MM, H * MM))
    ax = P.axmm(fig, W, H, 12, 14, 70, 33); P.letter(fig, W, H, 2, 3, "a")
    m = D1["inert"]
    for k, col, nm in (("random", RED, "random halves"), ("guide", BLUE, "guide halves")):
        v = D1["cross_%s_rest" % k][m]
        ax.hist(np.clip(v, -15, 30), bins=np.linspace(-15, 30, 61), color=col, alpha=0.75, histtype="stepfilled", lw=0,
                label="%s: %d%% above 0" % (nm, round(100 * (v > 0).mean())))
    ax.axvline(0, color="k", lw=0.6, ls="--"); P.edge_ticks(ax, -15, 30, [-15, -10, 0, 10, 20, 30])
    ax.set_xlabel("cross-half signal (half A · half B)"); ax.set_ylabel("targets")
    P.title(ax, "CD4 T cells, unexpressed targets"); ax.set_ylim(0, ax.get_ylim()[1] * 1.3); ax.legend(fontsize=5.6, loc="upper right", handlelength=1.0)
    ax = P.axmm(fig, W, H, 102, 14, 78, 33); P.letter(fig, W, H, 92, 3, "b")
    ds = J(v1, "donor_split.json")["C_null_baseline_donor_split"]["inert"]
    groups = [("D1\nrest", D1, "rest"), ("D1\nstim", D1, "stim"), ("D2\nrest", D2, "rest"), ("D2\nstim", D2, "stim")]
    x = np.arange(4); w = 0.36
    for j, (k, col) in enumerate((("random", RED), ("guide", BLUE))):
        for xi, (_, dd, s) in zip(x, groups):
            ind = dd["cross_%s_%s" % (k, s)][dd["inert"]] > 0
            ax.bar(xi + (j - 0.5) * w, 100 * ind.mean(), w, color=col); P.err(ax, xi + (j - 0.5) * w, ind.mean(), *P.share_ci(ind.mean(), ind.size), scale=100)
    for xi, s in ((4.15, "rest"), (4.85, "stim")):
        p = ds["%s_frac_pos" % s]; ax.bar(xi, 100 * p, w, color=BLUE); P.err(ax, xi, p, *P.share_ci(p, 1479), scale=100)
    ax.axhline(50, color="k", lw=0.6, ls="--"); ax.set_ylim(0, 118); ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xticks(list(x) + [4.15, 4.85]); ax.set_xticklabels([g for g, _, _ in groups] + ["donors\nrest", "donors\nstim"], fontsize=5.6)
    ax.set_ylabel("Targets with positive\ncross-half signal (%)")
    ax.legend(handles=P.patches(((RED, "random halves"), (BLUE, "guide or donor halves"))), fontsize=5.6, loc="upper center", ncol=2, handlelength=1, columnspacing=1)
    P.title(ax, "Unexpressed targets")
    P.letter(fig, W, H, 2, 62, "c")
    ntc = J(v1, "ntc_null.json")
    cl = dict(J(v1, "cl_ntc_replogle.json")); cl.update(J(v1, "cl_ntc_nadig.json")); cl.update(J(v1, "cl_ntc_h1.json"))
    rows = [("CD4\nD1", ntc["D1_rest"]), ("CD4\nD2", ntc["D2_rest"])]
    for k, lab in (("k562_ess", "K562\nessential"), ("k562_gw1", "K562\nGW"), ("rpe1", "RPE1"), ("jurkat", "Jurkat"), ("hepg2", "HepG2"), ("h1", "H1\nVCC 2025")):
        rows.append((lab, cl[k]))
    x = np.arange(len(rows)); w = 0.36
    axl = P.axmm(fig, W, H, 13, 72, 82, 30)
    axl.bar(x - w / 2, [100 * r[1]["random_frac_pos"] for r in rows], w, color=RED)
    axl.bar(x + w / 2, [100 * r[1]["pairs_frac_pos_mean"] for r in rows], w, color=BLUE)
    for xi, r in zip(x, rows):
        P.err(axl, xi - w / 2, r[1]["random_frac_pos"], *r[1]["random_frac_pos_ci95"], scale=100)
        P.err(axl, xi + w / 2, r[1]["pairs_frac_pos_mean"], *r[1]["pairs_frac_pos_range"], scale=100)
    axl.axhline(50, color="k", lw=0.6, ls="--"); axl.set_ylim(0, 118); axl.set_yticks([0, 25, 50, 75, 100])
    axl.set_xticks(x); axl.set_xticklabels([r[0] for r in rows], fontsize=5.4); axl.set_ylabel("Guides with positive\ncross-half signal (%)")
    axl.legend(handles=P.patches(((RED, "random halves"), (BLUE, "different guides"))), fontsize=5.6, loc="upper right", ncol=2, handlelength=1, columnspacing=1)
    P.title(axl, "Non-targeting guides")
    P.letter(fig, W, H, 104, 62, "d")
    c1, c2 = J(res, "depth_grid", "depth_grid_D1.json"), J(res, "depth_grid", "depth_grid_D2.json")
    depth_panels(fig, W, H, c1, c2, 114, 72, 27, 30, 13)
    key = [Line2D([], [], color=RED, lw=0.8, label="random"), Line2D([], [], color=BLUE, lw=0.8, label="guide"),
           Line2D([], [], color="k", marker="o", ls="none", ms=2.6, label="D1"), Line2D([], [], color="k", marker="s", ls="none", ms=2.6, label="D2"),
           Line2D([], [], color="k", lw=0.7, label="rest"), Line2D([], [], color="k", lw=0.7, ls="--", label="stim")]
    fig.legend(handles=key, loc="upper left", bbox_to_anchor=(112 / W, 1 - 113.5 / H), ncol=6, fontsize=5.2, handlelength=1.4, columnspacing=0.8, frameon=False)
    P.save(fig, out, "fig2_next")


def draw_halves(ax):
    """schematic: one perturbation's 20 cells from 5 clones (9, 5, 3, 2, 1 cells); random halves split every clone,
    clone halves keep each clone whole. Each clone is a dashed ring, as in Fig. 1; the splits are fixed, not drawn at random."""
    ax.set_xlim(0, 60); ax.set_ylim(0, 44); ax.axis("off")
    cols = ["#2e0f3f", "#b07cc6", "#5b2c6f", "#ead8f0", "#8e44ad"]   # clones of sizes 9, 5, 3, 2, 1
    d = 1.9; s3 = np.sqrt(3) / 2
    OFF = {1: [(0, 0)], 2: [(-d / 2, 0), (d / 2, 0)], 3: [(-d / 2, -d * s3 / 3), (d / 2, -d * s3 / 3), (0, 2 * d * s3 / 3)],
           4: [(-d / 2, -d / 2), (d / 2, -d / 2), (-d / 2, d / 2), (d / 2, d / 2)],
           5: [(1.62 * np.sin(2 * np.pi * i / 5), 1.62 * np.cos(2 * np.pi * i / 5)) for i in range(5)],
           9: [(x * d, y * d) for y in (1, 0, -1) for x in (-1, 0, 1)]}

    def clone(cx, cy, k, n):
        off = OFF[n]; rr = max(np.hypot(*o) for o in off) + 1.25
        ax.add_patch(plt.Circle((cx, cy), rr, fill=False, ec="#aab7b8", lw=0.5, ls=(0, (2, 1.5))))
        for ox, oy in off:
            ax.add_patch(plt.Circle((cx + ox, cy + oy), 0.9, fc=cols[k], ec="white", lw=0.4))
        return rr

    splits = (((0, 5), (1, 3), (2, 2), (3, 1), (4, 1)), ((0, 4), (1, 2), (2, 1), (3, 1))), \
             (((0, 9), (3, 2)), ((1, 5), (2, 3), (4, 1)))
    for row, (y, lab, rule, col, (A, B)) in enumerate(((25, "Random halves", "coin per cell", RED, splits[0]),
                                                        (2, "Clone halves", "coin per clone", BLUE, splits[1]))):
        ax.text(0, y + 17.2, lab, fontsize=6.3, fontweight="bold", va="bottom")
        ax.text(21 if row == 0 else 19.5, y + 17.3, "(" + rule + ")", fontsize=5.6, va="bottom")
        ax.plot([0, 59], [y + 16.6, y + 16.6], color=col, lw=1.0)
        for grp, x0, nm in ((A, 0.5, "half A"), (B, 31, "half B")):
            ax.add_patch(plt.Rectangle((x0, y), 28, 12.5, fill=False, lw=0.6, ec="#777777"))
            ax.text(x0, y + 13.0, "%s, %d cells" % (nm, sum(n for _, n in grp)), fontsize=5.6, va="bottom")
            widths = [2 * (max(np.hypot(*o) for o in OFF[n]) + 1.25) for _, n in grp]
            gap = (28 - sum(widths)) / (len(grp) + 1); x = x0 + gap
            for (k, n), w in zip(grp, widths):
                clone(x + w / 2, y + 6.25, k, n); x += w + gap
    ax.text(59, 24.3, "halves share clones", fontsize=5.6, ha="right", va="top")
    ax.text(59, 1.3, "halves share no clone", fontsize=5.6, ha="right", va="top")


def fig3(v1, res, out):
    D1 = np.load(os.path.join(v1, "cd4_D1.npz"))
    st = J(res, "choose_stats.json"); per = st["perturbations"]
    W, H = 183, 168
    fig = plt.figure(figsize=(W * MM, H * MM))
    ax = P.axmm(fig, W, H, 0, 0, 89, 52); ax.set_xlim(0, 89); ax.set_ylim(0, 52); ax.axis("off"); P.PS.draw_fig3a(ax)
    ax = P.axmm(fig, W, H, 104, 14, 74, 33); P.letter(fig, W, H, 95, 3, "b")
    z, zp = D1["trbc_z_rest"], D1["trbc_zperm_rest"]; ok = D1["trbc_n_rest"] >= 30
    bins = np.linspace(-25, 25, 81)
    ax.hist(np.clip(zp[ok], -25, 25), bins=bins, color=GREY, alpha=0.8, label="perturbation labels permuted", histtype="stepfilled", lw=0)
    ax.hist(np.clip(z[ok], -25, 25), bins=bins, color=ORANGE, alpha=0.75, label="observed", histtype="stepfilled", lw=0)
    ax.set_yscale("log"); P.edge_ticks(ax, -25, 25, [-25, -10, 0, 10, 25])
    ax.set_yticks([1, 10, 100, 1000]); ax.set_yticklabels(["1", "10", "100", "1,000"]); ax.minorticks_off()
    ax.set_xlabel("TRBC1 share of a perturbation (z-score)"); ax.set_ylabel("perturbations"); ax.legend(fontsize=5.6, loc="upper left")
    P.title(ax, "TRBC1 share per perturbation")
    P.letter(fig, W, H, 2, 58, "c")
    for j, (k, col, nm) in enumerate((("random", RED, "random halves"), ("guide", BLUE, "guide halves"))):
        ax = P.axmm(fig, W, H, 12 + j * 46, 70, 33, 30)
        a, b = D1["trbc_dA_%s_rest" % k], D1["trbc_dB_%s_rest" % k]
        okk = (D1["trbc_nA_%s_rest" % k] >= 15) & (D1["trbc_nB_%s_rest" % k] >= 15)
        ax.scatter(a[okk], b[okk], s=1.2, c=col, alpha=0.35, lw=0, rasterized=True)
        lo, hi = -0.35, 0.56; ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.axhline(0, color="k", lw=0.4); ax.axvline(0, color="k", lw=0.4)   # covers all points
        r = np.corrcoef(a[okk], b[okk])[0, 1]
        ax.set_xlabel("TRBC1 share deviation, half A"); ax.set_ylabel("TRBC1 share deviation, half B" if j == 0 else "")
        ax.set_title("%s (r = %s)" % (nm[0].upper() + nm[1:], ("%.2f" % r).replace("-", "−")), loc="left", fontsize=6.5, pad=3)
    # d: TRBC halves (CD4: the strong system)
    ax = P.axmm(fig, W, H, 116, 70, 62, 30); P.letter(fig, W, H, 106, 58, "d")
    labs = [("D1_Rest", "D1\nrest"), ("D1_Stim48", "D1\nstim"), ("D2_Rest", "D2\nrest"), ("D2_Stim48", "D2\nstim")]
    x = np.arange(4); w = 0.36
    for xi, (k, lab) in zip(x, labs):
        g = J(res, "trbc_split", k + ".json")["groups"]["inert"]; n = g["n_targets"]
        for j, (key, col) in enumerate((("share_pos_random_r1", RED), ("share_pos_trbc", BLUE))):
            p = g[key]; ax.bar(xi + (j - 0.5) * w, 100 * p, w, color=col); P.err(ax, xi + (j - 0.5) * w, p, *P.share_ci(p, n), scale=100)
        ax.text(xi, 106, "%d%%" % round(100 * g["mean_trbc"] / g["mean_random_r1"]), ha="center", fontsize=5.6)
    ax.axhline(50, color="k", lw=0.6, ls="--"); ax.set_ylim(0, 118); ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xticks(x); ax.set_xticklabels([l for _, l in labs], fontsize=5.6); ax.set_ylabel("Targets with positive\ncross-half signal (%)")
    P.title(ax, "CD4, unexpressed targets")
    fig.legend(handles=P.patches(((RED, "random halves"), (BLUE, "TRBC halves (no shared clone)"))), fontsize=5.6, loc="upper left",
               bbox_to_anchor=(116 / W, 1 - 110 / H), ncol=2, handlelength=1, columnspacing=1, frameon=False)
    # e: clone-barcode halves, schematic; f: organoid perturbations (the weak system)
    ax = P.axmm(fig, W, H, 4, 122, 58, 44); draw_halves(ax); P.letter(fig, W, H, 2, 119, "e")
    ax = P.axmm(fig, W, H, 110, 124, 60, 36); P.letter(fig, W, H, 96, 119, "f")
    r = np.array([v["random"] for v in per.values()]); c = np.array([v["clone"] for v in per.values()])
    sh = np.array([v["perm_mean"] for v in per.values()])
    # shuffled clone labels: the "clone" halves are random groupings that still share clones (red, open)
    ax.scatter(r, sh, s=9, facecolors="none", edgecolors=RED, lw=0.6, zorder=3, label="clone labels shuffled")
    ax.scatter(r, c, s=8, color=BLUE, lw=0, zorder=4, label="true clones")
    ax.plot([-1, 18], [-1, 18], color="k", lw=0.5, ls="--"); ax.axhline(0, color="k", lw=0.4)
    ax.set_xlim(0, 17); ax.set_ylim(-1, 17); ax.set_xlabel("random halves"); ax.set_ylabel("clone halves")
    ax.legend(fontsize=5.6, loc="upper left", handletextpad=0.3, borderaxespad=0.2)
    P.title(ax, "Organoids, cross-half signal")
    P.save(fig, out, "fig3_next")


def fig4(v1, res, out):
    st = J(res, "choose_stats.json"); per = st["perturbations"]; V = st["V"]
    W, H = 183, 100
    fig = plt.figure(figsize=(W * MM, H * MM))
    # a: organoids, strength against sum p^2, prediction from the measured ICC
    ax = P.axmm(fig, W, H, 14, 10, 40, 30); P.letter(fig, W, H, 2, 3, "a")
    r = np.array([v["random"] for v in per.values()]); c = np.array([v["clone"] for v in per.values()])
    p2 = np.array([v["sum_p2"] for v in per.values()]); s = (r - c) / V
    icc = st["slope_vs_icc"]["median_icc"]; lo, hi = st["slope_vs_icc"]["icc_ci"]
    xx = np.linspace(0, p2.max() * 1.05, 50)
    ax.fill_between(xx, lo * xx, hi * xx, color="#d5dbdb", lw=0, zorder=1, label="ICC × Σp² (95% interval)")
    ax.plot(xx, icc * xx, color="k", lw=0.7, zorder=2)
    ax.scatter(p2, s, s=8, color="k", lw=0, zorder=3, label="perturbations")
    ax.set_xlabel("Σ p_k² (= 1 / effective clones)"); ax.set_ylabel("clone strength")
    ax.legend(fontsize=5.6, loc="upper left"); ax.set_xlim(0, p2.max() * 1.05)
    P.title(ax, "Organoid perturbations")
    # b: CD4 donors, two clone readouts of the same K_eff
    cl = J(res, "clonality.json")
    ax = P.axmm(fig, W, H, 14, 62, 40, 30); P.letter(fig, W, H, 2, 55, "b")
    xs, ys = [], []
    for d in ("D1", "D2", "D3", "D4"):
        for st_, key, mk in (("rest", "Rest", "o"), ("stim48", "Stim48hr", "s")):
            # expression readout: the unified clone strength of Fig. 4c (geometric mean of the two references)
            u = J(res, "unified", "cd4_%s_%s.json" % (d, key))
            yv = float(np.sqrt(u["ruler_all"]["null_clone_scaled"] * u["ruler_null"]["null_clone_scaled"]))
            xv = cl[d][st_]["vc"]; xs.append(xv); ys.append(yv)
            ax.scatter([xv], [yv], s=12, marker=mk, color="k", lw=0, zorder=3)
            if st_ == "rest":
                ax.text(xv - 0.025, yv, d, ha="right", va="center", fontsize=5.6)
    k = float(np.sum(np.array(xs) * np.array(ys)) / np.sum(np.array(xs) ** 2))
    xx = np.linspace(0, 0.6, 20); ax.plot(xx, k * xx, color="k", lw=0.6, ls="--")
    ax.set_xlim(0, 0.6); ax.set_ylim(0, 0.015); ax.set_yticks([0, 0.005, 0.010, 0.015]); ax.set_yticklabels(["0", "0.005", "0.010", "0.015"])
    ax.set_xlabel("TRBC readout (clone variance of a guide mean)"); ax.set_ylabel("expression readout\n(clone strength, as in c)")
    ax.legend(handles=[Line2D([], [], marker="o", ls="none", color="k", ms=3, label="resting"), Line2D([], [], marker="s", ls="none", color="k", ms=3, label="stimulated")],
              fontsize=5.6, loc="upper left")
    P.title(ax, "CD4, four donors")
    # c: unified clone strength across datasets
    ax = P.axmm(fig, W, H, 112, 10, 66, 82); P.letter(fig, W, H, 66, 3, "c")
    U = os.path.join(res, "unified"); G = os.path.join(res, "generality")
    def rng_of(files, key="null_clone_scaled"):
        v = [J(f)["ruler_" + m][key] for f in files for m in ("all", "null")]
        return min(v), max(v)
    rows = [("Human CAR-T, in mice 21 d", rng_of([os.path.join(G, "knudsen.json"), os.path.join(G, "knudsen_lineage.json")]), "p"),
            ("Mouse P14 T, LCMV mice 8 d", rng_of([os.path.join(G, "fagerberg.json")]), "p"),
            ("CD4 T, genome-wide (4 donors)", rng_of([os.path.join(U, "cd4_%s_%s.json" % (d, s_)) for d in ("D1", "D2", "D3", "D4") for s_ in ("Rest", "Stim48hr")]), "p"),
            ("Primary T, CRISPRa, 146 guides", rng_of([os.path.join(G, "schmidt_resting.json"), os.path.join(G, "schmidt_restim.json")]), "p"),
            ("Mouse OT-I T, tumour mice 7 d", rng_of([os.path.join(G, "zhou.json")]), "p"),
            ("Primary AML, in vitro", rng_of([os.path.join(G, "cao_invitro.json")]), "p")]
    jin = J(res, "jin2020_clone.json")["classes"]; jv = [v["summary"]["clone_strength"] for v in jin.values()]
    rows += [("Embryonic mouse cortex*", (min(jv), max(jv)), "p"), ("Brain organoids*", (st["permutation"]["T"],) * 2, "p")]
    rows += [("HepG2", rng_of([os.path.join(U, "hepg2.json")]), "l"), ("Jurkat", rng_of([os.path.join(U, "jurkat.json")]), "l"),
             ("RPE1", rng_of([os.path.join(U, "rpe1.json")]), "l"),
             ("Melanoma (3 treatments)", rng_of([os.path.join(U, "frangieh_%s.json" % t) for t in ("Control", "IFNγ", "Co-culture")]), "l"),
             ("K562 essential", rng_of([os.path.join(U, "k562ess.json")]), "l"), ("K562 genome-wide", rng_of([os.path.join(U, "k562gw.json")]), "l"),
             ("H1 (VCC 2025)", rng_of([os.path.join(U, "vcc2025_h1.json")]), "l")]
    wb = J(res, "wholebrain_wb_results.json")
    wbv = [(wb[c_]["within"]["mean_random"] - wb[c_]["within"]["mean_diff_guides"]) * wb[c_]["within"]["mean_random_scaled"] / wb[c_]["within"]["mean_random"] for c_ in ("Glut", "GABA")]
    rows.append(("Mouse neurons, AAV at P16", (min(wbv), max(wbv)), "n"))
    for i, (lab, (a, b), g) in enumerate(rows):
        y = len(rows) - 1 - i; open_ = lab.endswith("*")
        ax.plot([max(a, 5e-5), b], [y, y], color="k", lw=1.6, solid_capstyle="round")
        ax.scatter([np.sqrt(max(a, 5e-5) * b)], [y], s=10, color="white" if open_ else "k", edgecolor="k", lw=0.8, zorder=3)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r_[0] for r_ in rows][::-1], fontsize=5.6)
    n_line = sum(1 for r_ in rows if r_[2] == "l")
    for yb in (0.5, 0.5 + n_line):
        ax.axhline(yb, color="#bbbbbb", lw=0.5)
    for y0, y1, lab in ((0, 0, "non-dividing\nneurons"), (1, n_line, "cell lines"), (n_line + 1, len(rows) - 1, "primary, in vivo\nor organoid")):
        ax.text(0.17, (y0 + y1) / 2, lab, ha="right", va="center", fontsize=5.6)
    ax.set_xscale("log"); ax.set_xlim(5e-5, 0.2); ax.minorticks_off()
    ax.set_xticks([1e-4, 1e-3, 1e-2, 1e-1]); ax.set_xticklabels(["0.0001", "0.001", "0.01", "0.1"])
    ax.set_xlabel("clone strength (shared deviation / per-cell variance)")
    P.title(ax, "16 screens")
    P.save(fig, out, "fig4_next")


def fig5(v1, res, out):
    """v1 Fig. 4 (ceilings and reliability labels) preceded by the CD4 composition of the random-half signal"""
    ce = J(v1, "celleval.json")["shared"]
    W, H = 183, 126
    fig = plt.figure(figsize=(W * MM, H * MM))
    # a: composition of the random-half signal (the TRBC-halves check is reported in the text, on other cells and targets)
    ax = P.axmm(fig, W, H, 14, 14, 34, 42); P.letter(fig, W, H, 2, 3, "a")
    c1, c2 = J(v1, "clone_D1.json"), J(v1, "clone_D2.json"); gd = J(v1, "guide_donor.json")
    bars = []
    for lab, dd, stt, rk, tk in (("D1\nrest", c1, "rest", "rest", "D1_Rest"), ("D1\nstim", c1, "stim48", "stim", "D1_Stim48"), ("D2\nrest", c2, "rest", "rest", "D2_Rest"), ("D2\nstim", c2, "stim48", "stim", "D2_Stim48")):
        R = dd[stt]["depth"]["random_half"]["expressed"]["cross_per_target"]; Gg = dd[stt]["depth"]["guide_half"]["expressed"]["cross_per_target"]
        rep = Gg / R; gsp = gd["guide_specific_over_diff_guide_" + rk] * rep
        bars.append((lab, 1 - rep - gsp, gsp, rep))
    x = np.arange(len(bars)); w = 0.6
    for xi, (lab, tr, gsp, rep) in enumerate(bars):
        ax.bar(xi, 100 * rep, w, color=BLUE); ax.bar(xi, 100 * gsp, w, bottom=100 * rep, color="#e6b0aa")
        ax.bar(xi, 100 * tr, w, bottom=100 * (rep + gsp), color=RED)
        ax.text(xi, 100 * (rep + gsp) + 100 * tr / 2, "%d%%" % round(100 * tr), ha="center", va="center", fontsize=5.6, color="white")
    ax.set_xticks(x); ax.set_xticklabels([b[0] for b in bars], fontsize=5.6)
    ax.set_ylim(0, 100); ax.set_yticks([0, 25, 50, 75, 100]); ax.set_ylabel("Share of random-half signal (%)")
    P.title(ax, "CD4, expressed targets")
    fig.legend(handles=P.patches(((RED, "within one transduction"), ("#e6b0aa", "specific to the guide"), (BLUE, "reproduced by the other guide")))
               , fontsize=5.6, loc="upper left", bbox_to_anchor=(4 / W, 1 - 66 / H),
               handlelength=0.9, labelspacing=0.25, frameon=False)
    # b: cell-eval ceilings (v1 Fig. 4a)
    mets = [("overlap_at_N", "DE overlap"), ("precision_at_N", "DE precision"), ("de_sig_genes_recall", "DE recall"),
            ("de_direction_match", "DE direction match*"), ("pearson_delta", "Pearson Δ"), ("discrimination_score_l1", "Discrimination (L1)*"),
            ("discrimination_score_cosine", "Discrimination (cosine)*")]
    P.letter(fig, W, H, 58, 3, "b")
    for j, (sub, nm) in enumerate((("all", "All 300"), ("inert", "150 unexpressed"), ("strong", "150 strength-matched"))):
        ax = P.axmm(fig, W, H, 92 + j * 30, 14, 27, 42)
        for i, (m, lab) in enumerate(mets):
            r_ = np.nanmean([ce["random_%d" % s_][sub]["ceiling"].get(m) or np.nan for s_ in (0, 1, 2)])
            gv = [x_ for x_ in (ce[k][sub]["ceiling"].get(m) for k in ("guide_0real_1pred_0", "guide_1real_0pred_0")) if x_ is not None and not np.isnan(x_)]
            g = float(np.mean(gv)) if gv else np.nan
            y = len(mets) - 1 - i
            if np.isnan(g):
                ax.plot([0, r_], [y, y], color="#d5dbdb", lw=1.6, zorder=1); ax.scatter([0], [y], s=12, facecolor="white", edgecolor=BLUE, zorder=3, lw=0.8)
            else:
                ax.plot([g, r_], [y, y], color="#d5dbdb", lw=1.6, zorder=1); ax.scatter([g], [y], s=12, color=BLUE, zorder=3, lw=0)
            ax.scatter([r_], [y], s=12, color=RED, zorder=3, lw=0)
        ax.axvline(2 * 0.5 / 1.5, color=GREY, lw=0.6, ls=":")
        ax.set_xlim(-0.03, 1.02); ax.set_ylim(-0.6, len(mets) - 0.4); ax.set_yticks(range(len(mets)))
        ax.set_yticklabels([l for _, l in mets][::-1] if j == 0 else []); ax.set_xticks([0, 0.5, 1]); ax.set_xticklabels(["0", "0.5", "1"])
        ax.set_xlabel("noise ceiling"); ax.set_title(nm, loc="left", fontsize=6.5, pad=3)
    fig.legend(handles=[Line2D([], [], marker="o", ls="none", color=RED, ms=3, label="random halves"),
                        Line2D([], [], marker="o", ls="none", color=BLUE, ms=3, label="guide halves")], fontsize=5.6, loc="upper left",
               bbox_to_anchor=(120 / W, 1 - 66 / H), ncol=2, frameon=False)
    # c: reliability labels (v1 Fig. 4b)
    ax = P.axmm(fig, W, H, 14, 92, 164, 22); P.letter(fig, W, H, 2, 84, "c")
    w = 0.36; cats = ("all", "inert", "expressed")
    for g_, (fname, nn, x0) in enumerate((("ci_cd4.json", {"all": 6615, "inert": 1235, "expressed": 1097}, 0.0),
                                           ("ci_cd4_D2.json", {"all": 5003, "inert": 869, "expressed": 859}, 3.6))):
        ci = J(v1, fname)["wang_reliable_share_rest"]; xg = x0 + np.arange(3)
        for j, (k, col, nm) in enumerate((("random", RED, "random halves"), ("guide", BLUE, "guide halves"))):
            ax.bar(xg + (j - 0.5) * w, [100 * ci[k][s_] for s_ in cats], w, color=col)
            for xi, s_ in zip(xg, cats):
                P.err(ax, xi + (j - 0.5) * w, ci[k][s_], *P.share_ci(ci[k][s_], nn[s_]), scale=100)
        ax.text(x0 + 1, -15, "donor %d" % (g_ + 1), ha="center", va="top", fontsize=6.5)
    ax.set_xticks([0, 1, 2, 3.6, 4.6, 5.6]); ax.set_xticklabels(["all", "unexpressed", "expressed"] * 2); ax.set_ylabel("Labelled reliable (%)")
    ax.set_ylim(0, 62); P.title(ax, "Reliability labels, resting CD4 T cells")
    P.save(fig, out, "fig5_next")


def regroup(pg, groups):
    rr = [r for r in pg if r[1] in groups]; tg = sorted({r[0] for r in rr})
    d = np.array([np.mean([r[2] - r[3] for r in rr if r[0] == t]) for t in tg])
    return float(np.mean([r[2] - r[3] for r in rr])), float(d.std(ddof=1) / np.sqrt(len(tg)))


def fig7(res, notes, out):
    W, H = 183, 74
    fig = plt.figure(figsize=(W * MM, H * MM))
    # a: donor-1 model scores, unexpressed targets, from the pilot summary table
    rows = []
    for line in open(notes, encoding="utf-8"):
        c = [x.strip() for x in line.strip().strip("|").split("|")]
        if len(c) > 3 and re.match(r"^-?\d", c[1]):
            rows.append((c[0], float(c[1]), float(c[2])))
    keep = [("Context mean", "baseline (context mean)"), ("Copy resting effect", "copy resting effect"), ("Linear (ridge)", "linear (ridge)"),
            ("Arc perturb_mean", "Arc perturb_mean"), ("Arc pseudobulk", "Arc pseudobulk"), ("STATE Transition (HVG)", "STATE"),
            ("State (SE-600M + ST)", "State (SE-600M + STATE)"), ("PIE (expression shift)", "PIE"), ("Stack (zero-shot)", "Stack")]
    rd = {r[0]: r for r in rows}
    ax = P.axmm(fig, W, H, 36, 12, 40, 52); P.letter(fig, W, H, 2, 3, "a")
    for i, (k, lab) in enumerate(keep):
        y = len(keep) - 1 - i; a, b = rd[k][1], rd[k][2]
        ax.plot([b, a], [y, y], color="#d5dbdb", lw=1.6, zorder=1)
        ax.scatter([a], [y], s=12, color=RED, zorder=3, lw=0); ax.scatter([b], [y], s=12, color=BLUE, zorder=3, lw=0)
    base = rd["Context mean"][1]; ax.axvline(base, color=GREY, lw=0.6, ls=":")
    ax.set_yticks(range(len(keep))); ax.set_yticklabels([l for _, l in keep][::-1], fontsize=5.4)
    ax.set_xlim(-0.03, 0.5); ax.set_xlabel("Pearson r, predicted vs answer")
    ax.legend(handles=[Line2D([], [], marker="o", ls="none", color=RED, ms=3, label="answer: same guide"),
                       Line2D([], [], marker="o", ls="none", color=BLUE, ms=3, label="answer: other guide")], fontsize=5.0, loc="lower right")
    P.title(ax, "Donor 1, unexpressed targets")
    # b: clone gain against donor clonality
    cl = J(res, "clonality.json"); det = {d: J(res, "detail_%s.json" % d) for d in ("D1", "D2", "D3", "D4")}
    vc = {d: (cl[d]["rest"]["vc"] + cl[d]["stim48"]["vc"]) / 2 for d in det}
    ax = P.axmm(fig, W, H, 92, 12, 34, 44); P.letter(fig, W, H, 84, 3, "b")
    order = sorted(det, key=lambda d: vc[d])
    for m in ("copy", "linear", "state_transition"):
        xx = [vc[d] for d in order]; yy = [det[d]["models"][m]["all"]["AB"] for d in order]; ee = [det[d]["models"][m]["all"]["AB_se"] for d in order]
        ax.plot(xx, yy, color=COL[m], lw=0.9, ls=LS[m]); ax.errorbar(xx, yy, yerr=ee, fmt=MK[m], ms=3.2, color=COL[m], mec="white", mew=0.5, elinewidth=0.7, capsize=0)
        ax.text(xx[-1] + 0.015, yy[-1] + {"copy": -0.012, "linear": 0, "state_transition": 0.012}[m], {"copy": "copy", "linear": "linear", "state_transition": "STATE"}[m], va="center", fontsize=5.6)
    for d in order:
        ax.text(vc[d], 0.56, d, ha="center", fontsize=5.6)
    ax.set_xlim(0.1, 0.66); ax.set_ylim(0, 0.6); ax.set_xlabel("TRBC readout (as in Fig. 4b)"); ax.set_ylabel("score difference,\nsame − other guide")
    P.title(ax, "Four donors")
    # c: answer from another donor (same guide sequence, separate transduction)
    xd = J(res, "cross_donor.json")["pairs"]
    ax = P.axmm(fig, W, H, 136, 12, 22, 44); P.letter(fig, W, H, 130, 3, "c")
    for j, m in enumerate(("copy", "linear", "state_transition")):
        v = [r[m]["all"]["AB"] for r in xd.values()]
        ax.errorbar([j], [np.mean(v)], yerr=[[np.mean(v) - min(v)], [max(v) - np.mean(v)]], fmt=MK[m], ms=3.2, color=COL[m], mec="white", mew=0.5, elinewidth=0.7, capsize=0)
    ax.set_xlim(-0.6, 2.6); ax.set_ylim(0, 0.6); ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["copy", "linear", "STATE"], fontsize=5.6); ax.set_yticklabels([])
    ax.axhline(0, color="k", lw=0.4); P.title(ax, "Answer: other donor")
    # d: matched cells per guide, CD4 vs melanoma
    ax = P.axmm(fig, W, H, 166, 12, 15, 44); P.letter(fig, W, H, 161, 3, "d")
    for j, m in enumerate(("copy", "linear")):
        cd = []
        for d in ("D1", "D2", "D3", "D4"):
            b = J(res, "baselines_ds48_75_%s.json" % d)["models"][m]; n = b["inert"]["n"] + b["strong"]["n"]
            cd.append((b["inert"]["n"] * b["inert"]["leak"] + b["strong"]["n"] * b["strong"]["leak"]) / n)
        fr = sum((J(res, "detail_fr_IFNg_f%d.json" % k)["models"][m]["per_guide"] for k in (0, 1)), [])
        fv, fe = regroup(fr, ("null", "gene"))
        ax.plot([j - 0.15] * 4, cd, MK[m], ms=3, color=COL[m], mec="white", mew=0.4)
        ax.errorbar([j + 0.15], [fv], yerr=[fe], fmt=MK[m], ms=3.2, mfc="white", color=COL[m], mew=0.8, elinewidth=0.7, capsize=0)
    ax.set_xlim(-0.6, 1.6); ax.set_ylim(0, 0.6); ax.set_xticks([0, 1]); ax.set_xticklabels(["copy", "linear"], fontsize=5.6); ax.set_yticklabels([])
    ax.legend(handles=[Line2D([], [], marker="o", ls="none", color="k", ms=3, label="CD4 donors"),
                       Line2D([], [], marker="o", ls="none", mfc="white", color="k", ms=3, label="melanoma")], fontsize=5.6, loc="upper center",
              bbox_to_anchor=(0.5, 1.0), handletextpad=0.2, borderaxespad=0)
    P.title(ax, "Equal cells")
    P.save(fig, out, "fig7_next")


def draw_fig6a(ax):
    """one clone in both states and both halves (CD4) against separately screened lines; neutral cards, clones purple"""
    from plot_schematics import cell, arrow, card, rule, label, cluster, G1, G2
    rng = np.random.default_rng(11)
    label(ax, 32, 36.5, "same clones in both states (CD4)", ha="center", fontsize=6.5, fontweight="bold"); rule(ax, 15, 49, 34.4, RED)
    label(ax, 19, 30.5, "resting", ha="center", fontsize=6); label(ax, 31, 30.5, "stimulated", ha="center", fontsize=6)
    for yy, nm in ((22.5, "half A"), (11.5, "half B")):
        label(ax, 6.5, yy, nm, fontsize=6)
        for xx in (19, 31):
            for j in range(3):
                cell(ax, xx - 2.6 + j * 2.6, yy, G1[0], r=1.05)
    ax.plot([25, 25], [8, 28], color="#d5d8dc", lw=0.5, zorder=0); ax.plot([13, 37], [17, 17], color="#d5d8dc", lw=0.5, zorder=0)
    arrow(ax, 37.5, 22.5, 41.5, 22.5, lw=0.8, ms=5); arrow(ax, 37.5, 11.5, 41.5, 11.5, lw=0.8, ms=5)
    card(ax, 42.2, 18.6, 17.5, 7.8, fc="white", ec="#4d4d4d", r=1, z=1)
    label(ax, 50.95, 22.5, "common part\n(mean of states)", ha="center", fontsize=5.6)
    card(ax, 42.2, 7.6, 17.5, 7.8, fc="white", ec="#4d4d4d", r=1, z=1)
    label(ax, 50.95, 11.5, "specific part\n(half-difference)", ha="center", fontsize=5.6)
    label(ax, 91.5, 36.5, "separate screens (K562, RPE1)", ha="center", fontsize=6.5, fontweight="bold"); rule(ax, 75, 108, 34.4, BLUE)
    cluster(ax, 77, 24, G1[1], 5, rng); cluster(ax, 77, 13, G2[1], 5, rng)
    label(ax, 83, 24, "K562 clones", fontsize=6); label(ax, 83, 13, "RPE1 clones", fontsize=6)


def fig6(v1, res, out):
    """v1 Fig. 5 redrawn: 'common' for the part shared across contexts ('shared' is kept for clones), no numbers in data
    areas, model panels split into e and f"""
    W, H = 183, 160
    fig = plt.figure(figsize=(W * MM, H * MM))
    ax = P.axmm(fig, W, H, 0, 0, 120, 50); ax.set_xlim(0, 120); ax.set_ylim(0, 50); ax.axis("off"); draw_fig6a(ax)
    P.letter(fig, W, H, 2, 1, "a")
    ci = J(v1, "ci_cd4.json")
    # b: TRBC1 deviation (mean f - f0) / sd_f = z / sqrt(n), resting against stimulated
    ax = P.axmm(fig, W, H, 136, 8, 40, 34); P.letter(fig, W, H, 124, 1, "b")
    D1 = np.load(os.path.join(v1, "cd4_D1.npz"))
    nr, ns = D1["trbc_n_rest"], D1["trbc_n_stim"]
    ok = (nr >= 30) & (ns >= 30) & np.isfinite(D1["trbc_z_rest"]) & np.isfinite(D1["trbc_z_stim"])
    # raw deviation (mean f - f0), as in Fig. 3c: z / sqrt(n) is in units of the per-cell s.d. of f, so multiply by it
    sd = J(res, "clonality.json")["D1"]
    xr = D1["trbc_z_rest"][ok] / np.sqrt(nr[ok]) * sd["rest"]["sd_f"]; ys = D1["trbc_z_stim"][ok] / np.sqrt(ns[ok]) * sd["stim48"]["sd_f"]
    ax.scatter(xr, ys, s=1.2, c=RED, alpha=0.35, lw=0, rasterized=True)
    lo, hi = -0.35, 0.56; ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.axhline(0, color="k", lw=0.4); ax.axvline(0, color="k", lw=0.4)
    ax.set_xlabel("TRBC1 share deviation, resting"); ax.set_ylabel("TRBC1 share deviation,\nstimulated")
    P.title(ax, "Perturbations, donor 1 (r = %.2f)" % np.corrcoef(xr, ys)[0, 1])
    # c: reproducible energy of the two parts; the specific share goes into the legend labels
    ax = P.axmm(fig, W, H, 14, 67, 44, 31); P.letter(fig, W, H, 2, 55, "c")
    x = np.arange(2); w = 0.36
    for j, (k, col, nm) in enumerate((("random_r1", RED, "random halves"), ("guide", BLUE, "guide halves"))):
        g = ci[k]["gamma_share"][0]
        ax.bar(x + (j - 0.5) * w, [ci[k]["beta_cross_per_target"], ci[k]["gamma_cross_per_target"]], w, color=col,
               label="%s (specific share %.2f)" % (nm, g))
    ax.set_xticks(x); ax.set_xticklabels(["common", "context-\nspecific"]); ax.set_ylabel("reproducible energy\nper target")
    ax.set_ylim(0, 9); ax.legend(fontsize=5.6, loc="upper right", handlelength=1)
    P.title(ax, "CD4 T cells, donor 1")
    # d: common-to-specific predictor
    ax = P.axmm(fig, W, H, 78, 67, 100, 31); P.letter(fig, W, H, 66, 55, "d")
    kd1, kd2, ci2, rp = J(v1, "kd_D1.json"), J(v1, "kd_D2.json"), J(v1, "ci_cd4_D2.json"), J(v1, "replogle_indep.json")
    bars = [("D1", ci["random_r1"]["readout_k20"][0], ci["guide"]["readout_k20"][0]),
            ("D1,\nmatched", kd1["kd_matched|cell"]["readout_k20"], kd1["kd_matched|guide"]["readout_k20"]),
            ("D2", ci2["random_r1"]["readout_k20"][0], ci2["guide"]["readout_k20"][0]),
            ("D2,\nmatched", kd2["kd_matched|cell"]["readout_k20"], kd2["kd_matched|guide"]["readout_k20"]),
            ("K562 to\nRPE1", rp["target_K8"]["pair_same_k20"], rp["target_K8"]["pair_indepK_k20"])]
    x = np.array([0, 1, 2, 3, 4.6]); w = 0.36
    ax.bar(x - w / 2, [b[1] for b in bars], w, color=RED, label="random halves")
    ax.bar(x + w / 2, [b[2] for b in bars], w, color=BLUE, label="independent replicates")
    ax.plot(x + w / 2, [b[2] for b in bars], "_", color=BLUE, ms=7, mew=1.4, zorder=4)   # keeps near-zero bars visible
    for xi, cc in ((0, ci), (2, ci2)):
        for dx, k in ((-w / 2, "random_r1"), (w / 2, "guide")):
            P.err(ax, xi + dx, *cc[k]["readout_k20"])
    ax.axvline(3.8, color="#bfbfbf", lw=0.5)
    ax.set_xticks(x); ax.set_xticklabels([b[0] for b in bars], fontsize=5.6); ax.set_ylabel("held-out explained share"); ax.set_ylim(0, 0.8)
    ax.axhline(0, color="k", lw=0.5)
    ax.legend(fontsize=5.6, loc="upper left", ncol=2, handlelength=1)
    P.title(ax, "Predicting the specific part from the common part")
    ax.text(1.5, -0.33, "CD4 donors (states share clones)", ha="center", fontsize=5.6, transform=ax.get_xaxis_transform())
    ax.text(4.6, -0.33, "separate screens", ha="center", fontsize=5.6, transform=ax.get_xaxis_transform())
    # e: model without true structure (points, so that zero values stay visible)
    sim = J(v1, "sim_v6.json")
    ax = P.axmm(fig, W, H, 14, 122, 58, 30); P.letter(fig, W, H, 2, 110, "e")
    sc = [("fab_none", "no clone\neffect"), ("fab_low_separate", "separate\nclones"), ("fab_low_shared", "shared\nclones")]
    x = np.arange(3)
    for dx, dsg, col, nm in ((-0.12, "random", RED, "random halves"), (0.12, "independent", BLUE, "independent replicates")):
        ax.plot(x + dx, [sim["%s|%s" % (k, dsg)]["readout_share"] for k, _ in sc], "o", color=col, ms=3.6, ls="none", label=nm)
    ax.axhline(0, color="k", lw=0.5); ax.set_xticks(x); ax.set_xticklabels([l for _, l in sc]); ax.set_xlim(-0.5, 2.5)
    ax.set_ylabel("explained share"); ax.set_ylim(-0.01, 0.13); ax.legend(fontsize=5.6, loc="upper left", handletextpad=0.2)
    P.title(ax, "Model, no true structure")
    # f: model with true structure 0.70, measured across independent replicates
    ax = P.axmm(fig, W, H, 100, 122, 64, 30); P.letter(fig, W, H, 88, 110, "f")
    for kind, mk, lab, al in (("iso", "o", "clone effects in\nrandom directions", 0.45), ("low", "s", "clone effects along\nthe response", 1.0)):
        ratio, ri = [], []
        for sa in ("0.2", "0.6", "1.2", "1.5", "2.0"):
            a = sim["atten_sa%s_%s|random" % (sa, kind)]; b = sim["atten_sa%s_%s|independent" % (sa, kind)]
            ratio.append(a["inert_cross_per_target"] / b["expressed_cross_per_target"]); ri.append(b["readout_share"])
        ax.plot(ratio, ri, "-" + mk, color=BLUE, ms=3.0, lw=0.8, alpha=al)
        ax.text(7.8, ri[-1], lab, fontsize=5.6, va="center", color="k")
    ax.axvspan(4.0, 7.5, color="#eaecee", zorder=0)
    ax.text(5.75, 1.02, "CD4 range", fontsize=5.6, ha="center", va="bottom", transform=ax.get_xaxis_transform())
    ax.axhline(0.70, color=GREY, lw=0.6, ls="--"); ax.text(7.8, 0.72, "true value", fontsize=5.6, va="center")
    ax.set_xlim(0, 7.6); ax.set_xlabel("clone-to-signal energy ratio"); ax.set_ylabel("explained share,\nindependent replicates")
    ax.set_ylim(0, 0.9); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
    P.title(ax, "Model, true explained share 0.70")
    P.save(fig, out, "fig6_next")


# ---- Fig. 1 (schematic) with panel a from plot_schematics.fig1 and new b, c ----
from plot_schematics import canvas, card, rule, cell, arrow, label, panel_letter, cluster, G1, G2, save as PS_save
from matplotlib.patches import Circle, Rectangle
def fig1_next(out):
    W, H = 183, 122
    fig, ax = canvas(W, H); rng = np.random.default_rng(4)
    # ---------- a ----------
    panel_letter(ax, 2, H - 1.5, "a")
    label(ax, 8, H - 4, "Pooled screen: random halves and independent replicates", fontsize=7)
    label(ax, 29, 108.5, "transduction and expansion", ha="center", fontsize=8, fontweight="bold")
    for gi, (gy, cols, gname) in enumerate(((97, G1, "guide 1"), (80, G2, "guide 2"))):
        label(ax, 12.5, gy, gname, ha="center", fontsize=6.5, color=cols[0], fontweight="bold")
        arrow(ax, 18.8, gy, 23.5, gy, c=cols[0], lw=0.7, ms=5)
        for k, c in enumerate(cols):
            # one dashed ring per clone (the descendants of one transduced cell); expansion is uneven: 9, 3 and 1 cells
            cx = (29.5, 39.5, 47.5)[k]; cy = gy
            n, spread = (9, 3.2) if k == 0 else ((3, 2.0) if k == 1 else (1, 0.0))
            cluster(ax, cx, cy, c, n, rng, spread=spread)
            ring = (4.6, 3.5, 2.3)[k]
            ax.add_patch(Circle((cx, cy), ring, fill=False, ec="#aab7b8", lw=0.5, ls=(0, (2, 1.5)), zorder=1))
            if gi == 0:
                label(ax, cx, cy - ring - 1.7, "clone %d" % (k + 1), ha="center", fontsize=5.6)
    label(ax, 29, 71.2, "one transduced cell → one clone", ha="center", fontsize=6)
    label(ax, 29, 68.6, "uneven expansion: a few clones dominate", ha="center", fontsize=6)
    # clones drawn as in the left box (dashed ring, one colour per clone), so that the reader can see which clones each half holds
    def ringed(cx, cy, c, n):
        if n == 0:
            return
        spread, ring = {1: (0.0, 1.9), 2: (1.6, 2.8), 3: (2.0, 3.3), 4: (2.4, 3.7), 5: (2.6, 4.0), 9: (3.2, 4.6)}[n]
        cluster(ax, cx, cy, c, n, rng, spread=spread)
        ax.add_patch(Circle((cx, cy), ring, fill=False, ec="#aab7b8", lw=0.5, ls=(0, (2, 1.5)), zorder=1))
    # random halves: every clone is split between the two halves (the single-cell clones fall on one side)
    label(ax, 88, 108.5, "random cell halves", ha="center", fontsize=8, color=RED, fontweight="bold"); rule(ax, 74, 102, 106.2, RED)
    halves = {62.5: ((5, 2, 1), (4, 1, 0)), 89.5: ((4, 1, 0), (5, 2, 1))}   # (guide 1 clones, guide 2 clones)
    for x0, nm in ((62.5, "half A"), (89.5, "half B")):
        card(ax, x0, 76, 24, 27, ec="#d5d8dc", z=0.5)
        for cols, ns, yy in ((G1, halves[x0][0], 96.5), (G2, halves[x0][1], 83.0)):
            for k in range(3):
                ringed(x0 + (6.5, 15.0, 20.8)[k], yy, cols[k], ns[k])
        label(ax, x0 + 12, 73.8, nm, ha="center", fontsize=6.5, color="#34495e")
    ax.plot([88, 88], [76, 103], color=RED, lw=0.8, ls=(0, (2, 1.5)), zorder=2)
    label(ax, 88, 69, "same clones in both halves", ha="center", fontsize=6.5, color=RED)
    # independent replicates: each side holds the clones of one guide only
    label(ax, 150.5, 108.5, "independent replicates", ha="center", fontsize=8, color=BLUE, fontweight="bold"); rule(ax, 133, 168, 106.2, BLUE)
    for x0, cols, nm in ((124.5, G1, "guide 1 cells"), (151.5, G2, "guide 2 cells")):
        card(ax, x0, 76, 25, 27, ec="#d5d8dc", z=0.5)
        for k, n in enumerate((9, 3, 1)):
            ringed(x0 + (7.0, 16.0, 21.6)[k], 89.5, cols[k], n)
        label(ax, x0 + 12.5, 73.8, nm, ha="center", fontsize=6.5, color="#34495e")
    label(ax, 150.5, 69, "no clone in common", ha="center", fontsize=6.5, color=BLUE)
    arrow(ax, 54.5, 89, 59.5, 89, lw=0.9)
    label(ax, 119, 89, "vs", ha="center", fontsize=8, color=GREY, fontweight="bold")
    # ---------- b: four uses ----------
    # neutral greys only: red and navy are reserved for "shares clones" / "no shared clone" in the data figures
    DK, LT = "#4d4d4d", "#a6a6a6"
    panel_letter(ax, 2, 61, "b")
    label(ax, 8, 58.5, "Uses that treat the two halves as independent", fontsize=7)
    uses = ["ceilings and\nreliability labels", "common and context-\nspecific components", "predicting the\nspecific component", "scoring transfer\nmodels"]
    for k, t in enumerate(uses):
        x0 = 2 + k * 24.0; gx, gy = x0 + 2, 24
        label(ax, x0 + 11.5, 50, t, ha="center", fontsize=5.8, fontweight="bold")
        if k in (0, 1):
            v = rng.normal(size=8); vb = v + rng.normal(scale=0.5 if k == 0 else 0.9, size=8)
            ax.plot([gx, gx + 18], [gy + 12, gy + 12], color=LT, lw=0.4)
            for j in range(8):
                ax.add_patch(Rectangle((gx + j * 1.0, gy + 12), 0.75, v[j] * 2.0, fc=DK, ec="none", zorder=3))
                ax.add_patch(Rectangle((gx + 10 + j * 1.0, gy + 12), 0.75, vb[j] * 2.0, fc=LT, ec="none", zorder=3))
            a_, b_ = ("half A", "half B") if k == 0 else ("context 1", "context 2")
            label(ax, gx + 3.5, gy + 3, a_, ha="center", fontsize=5.6); label(ax, gx + 13.5, gy + 3, b_, ha="center", fontsize=5.6)
            label(ax, gx + 9.5, gy - 3, "agreement → ceiling" if k == 0 else "mean and half-difference", ha="center", fontsize=5.6)
        else:
            l1, l2, bot = (("common\n(half A)", "specific\n(half B)", "held-out perturbations") if k == 2
                           else ("context 1\n(input)", "context 2\n(answer)", "same transduction"))
            card(ax, gx - 1.8, gy + 6, 9.8, 9, ec=DK); label(ax, gx + 3.1, gy + 10.5, l1, ha="center", fontsize=5.6)
            arrow(ax, gx + 8.2, gy + 10.5, gx + 11.0, gy + 10.5, lw=0.8, ms=6)
            card(ax, gx + 11.2, gy + 6, 9.8, 9, ec=DK); label(ax, gx + 16.1, gy + 10.5, l2, ha="center", fontsize=5.6)
            label(ax, gx + 9.5, gy - 3, bot, ha="center", fontsize=5.6)
    # ---------- c: datasets ----------
    panel_letter(ax, 100, 61, "c")
    label(ax, 106, 58.5, "Datasets (16 screens)", fontsize=7)
    x0, x1 = 101.5, 181
    ax.plot([x0, x1], [54, 54], color="k", lw=0.6)
    xs = [102, 136, 162]
    for x, t in zip(xs, ["system (screens)", "null perturbations", "direct clone test"]):
        label(ax, x, 51.3, t, fontsize=5.8, fontweight="bold")
    ax.plot([x0, x1], [48.6, 48.6], color="k", lw=0.4)
    rows = [("primary CD4 T cells, 4 donors (1)", "NT guides; unex-\npressed targets", "TRBC halves"),
            ("brain organoids (1)", "control guide", "clone barcodes"),
            ("T cells transferred into mice (3)", "control guides", "–"),
            ("other primary cells in vitro (2)", "control guides", "–"),
            ("mouse brain in vivo (2)", "control guides", "–"),
            ("cell lines (7)", "control guides", "–")]
    for i, r in enumerate(rows):
        y = 45.0 - i * 5.6
        for x, t in zip(xs, r):
            label(ax, x, y, t, fontsize=5.6)
    ax.plot([x0, x1], [y - 3.0, y - 3.0], color="k", lw=0.6)
    fig.set_size_inches(W * MM, (H - 9) * MM); ax.set_ylim(9, H)   # trim the empty strip below panels b and c
    PS_save(fig, out, "fig1_next")


def main():
    v1, res, out = sys.argv[1:4]; os.makedirs(out, exist_ok=True)
    notes = os.path.join(res, "model_pilot_summary_table.md")   # model_eval/pilot_summary.py on <results>/pilot
    fig1_next(out); fig2(v1, res, out); fig3(v1, res, out); fig4(v1, res, out); fig5(v1, res, out); fig6(v1, res, out); fig7(res, os.path.normpath(notes), out)


if __name__ == "__main__":
    main()
