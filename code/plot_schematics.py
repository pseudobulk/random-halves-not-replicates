"""Schematic panels for the clone-effect preprint. Coordinates are millimetres on the final page
(double column 183 mm). Palette fixed across the paper.
Outputs (PDF + PNG 300 dpi) in figures/preprint/: fig1_schematic, fig3a_trbc, fig5a_units, fig6_decision, and a
contact sheet schematics_overview.png for checking.
Usage: plot_schematics.py <out_dir>"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch, Rectangle, Polygon, RegularPolygon

plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"], "font.size": 7, "pdf.fonttype": 42, "svg.fonttype": "none"})
RED, BLUE, GREY, ORANGE, LIGHT = "#c0392b", "#2c3e50", "#7f8c8d", "#e67e22", "#f4f6f7"
G1 = ["#6c3483", "#9b59b6", "#c39bd3"]          # clones of guide 1 (purples)
G2 = ["#117a65", "#1abc9c", "#76d7c4"]          # clones of guide 2 (teals)
MM = 1 / 25.4


def canvas(w, h):
    fig = plt.figure(figsize=(w * MM, h * MM))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, w); ax.set_ylim(0, h); ax.axis("off"); ax.set_aspect("equal")
    return fig, ax


def card(ax, x, y, w, h, fc="white", ec="#d5d8dc", lw=0.6, r=1.6, z=0):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=%g" % r, fc=fc, ec=ec, lw=lw, zorder=z))


def cell(ax, x, y, c, r=1.15, z=3):
    ax.add_patch(Circle((x, y), r, fc=c, ec="white", lw=0.5, zorder=z))


def arrow(ax, x0, y0, x1, y1, c="#566573", lw=0.8, ms=6, rad=0.0, z=2, ls="-"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=ms, color=c, lw=lw,
                                 connectionstyle="arc3,rad=%g" % rad, zorder=z, linestyle=ls))


def label(ax, x, y, s, **k):
    k.setdefault("ha", "left"); k.setdefault("va", "center")
    if k.get("color") != "white":
        k["color"] = "k"                     # no coloured text; colour stays in the marks (borders, cells)
    ax.text(x, y, s, **k)


def panel_letter(ax, x, y, s):
    ax.text(x, y, s, fontsize=8, fontweight="bold", ha="left", va="top")


def cluster(ax, cx, cy, c, n, rng, spread=2.4):
    pts = []
    for _ in range(n):
        for _ in range(50):
            p = (cx + rng.uniform(-spread, spread), cy + rng.uniform(-spread * 0.8, spread * 0.8))
            if all((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 > 2.3 ** 2 for q in pts):
                break
        pts.append(p); cell(ax, p[0], p[1], c)
    return pts


def scatter_cells(ax, x0, y0, w, h, colors, rng, r=1.15):
    pts = []
    for c in colors:
        for _ in range(200):
            p = (rng.uniform(x0 + r, x0 + w - r), rng.uniform(y0 + r, y0 + h - r))
            if all((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 > (2 * r + 0.25) ** 2 for q in pts):
                break
        pts.append(p); cell(ax, p[0], p[1], c, r)


def fig1(out):
    W, H = 183, 122
    fig, ax = canvas(W, H); rng = np.random.default_rng(4)
    # ---------- a ----------
    panel_letter(ax, 2, H - 1.5, "a")
    label(ax, 8, H - 4, "Pooled screen: random halves and independent units", fontsize=7)
    card(ax, 4, 66, 50, 46)
    label(ax, 29, 108.5, "transduction and expansion", ha="center", fontsize=7, color="#34495e")
    for gi, (gy, cols, gname) in enumerate(((97, G1, "guide 1"), (80, G2, "guide 2"))):
        ax.add_patch(FancyBboxPatch((7, gy - 2.6), 11, 5.2, boxstyle="round,pad=0,rounding_size=1.2", fc="white", ec=cols[0], lw=0.8, zorder=2))
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
    # random halves
    card(ax, 60, 66, 56, 46, ec=RED, lw=0.9)
    label(ax, 88, 108.5, "random cell halves", ha="center", fontsize=8, color=RED, fontweight="bold")
    allc = [G1[0]] * 9 + [G1[1]] * 3 + [G1[2]] * 1 + [G2[0]] * 9 + [G2[1]] * 3 + [G2[2]] * 1
    perm = list(rng.permutation(allc)); hA, hB = perm[: len(perm) // 2], perm[len(perm) // 2:]
    for x0, cols, nm in ((62.5, hA, "half A"), (89.5, hB, "half B")):
        ax.add_patch(FancyBboxPatch((x0, 76), 24, 27, boxstyle="round,pad=0,rounding_size=1.2", fc="white", ec="#e5e8e8", lw=0.6, zorder=1))
        scatter_cells(ax, x0 + 0.5, 76.5, 23, 26, cols, rng)
        label(ax, x0 + 12, 73.8, nm, ha="center", fontsize=6.5, color="#34495e")
    ax.plot([88, 88], [76, 103], color=RED, lw=0.8, ls=(0, (2, 1.5)), zorder=2)
    label(ax, 88, 69, "same clones in both halves", ha="center", fontsize=6.5, color=RED)
    # independent units
    card(ax, 122, 66, 57, 46, ec=BLUE, lw=0.9)
    label(ax, 150.5, 108.5, "independent units", ha="center", fontsize=8, color=BLUE, fontweight="bold")
    for x0, cols, nm in ((124.5, [G1[0]] * 9 + [G1[1]] * 3 + [G1[2]] * 1, "guide 1 cells"), (151.5, [G2[0]] * 9 + [G2[1]] * 3 + [G2[2]] * 1, "guide 2 cells")):
        ax.add_patch(FancyBboxPatch((x0, 76), 25, 27, boxstyle="round,pad=0,rounding_size=1.2", fc="white", ec="#e5e8e8", lw=0.6, zorder=1))
        scatter_cells(ax, x0 + 0.5, 76.5, 24, 26, list(rng.permutation(cols)), rng)
        label(ax, x0 + 12.5, 73.8, nm, ha="center", fontsize=6.5, color="#34495e")
    label(ax, 150.5, 69, "no clone in common", ha="center", fontsize=6.5, color=BLUE)
    arrow(ax, 54.5, 89, 59.5, 89, lw=0.9)
    label(ax, 119, 89, "vs", ha="center", fontsize=8, color=GREY, fontweight="bold")
    # ---------- b ----------
    panel_letter(ax, 2, 61, "b")
    label(ax, 8, 58.5, "Uses of random halves", fontsize=7)
    uses = [("noise ceilings and\nreliability labels", ""), ("shared and context-\nspecific parts", ""), ("transfer across\ncontexts", "")]
    for k, (t, s) in enumerate(uses):
        x0 = 4 + k * 31.5
        card(ax, x0, 8, 29.5, 46)
        label(ax, x0 + 14.75, 49, t, ha="center", fontsize=6.5, fontweight="bold", color="#34495e")
        # mini graphic
        gx, gy = x0 + 4, 25
        if k == 0:
            v = rng.normal(size=10); vb = v + rng.normal(scale=0.5, size=10)
            for j in range(10):
                ax.add_patch(Rectangle((gx + j * 1.0, gy + 9), 0.75, v[j] * 2.2, fc=RED, ec="none", zorder=3))
                ax.add_patch(Rectangle((gx + 11 + j * 1.0, gy + 9), 0.75, vb[j] * 2.2, fc=RED, ec="none", alpha=0.6, zorder=3))
            label(ax, gx + 4.5, gy + 1.5, "half A", ha="center", fontsize=5.5, color=GREY); label(ax, gx + 15.5, gy + 1.5, "half B", ha="center", fontsize=5.5, color=GREY)
            label(ax, gx + 10.5, gy - 3.5, "r  →  Spearman–Brown", ha="center", fontsize=6)
        elif k == 1:
            ax.add_patch(Rectangle((gx + 2, gy + 2), 6, 13, fc=GREY, ec="none", alpha=0.5, zorder=3))
            ax.add_patch(Rectangle((gx + 2, gy + 15), 6, 5, fc=RED, ec="none", zorder=3))
            label(ax, gx + 10, gy + 8.5, "shared", fontsize=5.5, color=GREY); label(ax, gx + 10, gy + 17.5, "context-\nspecific", fontsize=5.5, color=RED)
            label(ax, gx + 10.5, gy - 3.5, "half A · half B", ha="center", fontsize=6)
        else:
            card(ax, gx, gy + 6, 8, 8, fc="white", ec=GREY, r=1)
            label(ax, gx + 4, gy + 10, "shared\n(A)", ha="center", fontsize=5.3)
            arrow(ax, gx + 8.5, gy + 10, gx + 12.5, gy + 10, lw=0.7, ms=5)
            card(ax, gx + 13, gy + 6, 9, 8, fc="white", ec=RED, r=1)
            label(ax, gx + 17.5, gy + 10, "specific\n(B)", ha="center", fontsize=5.3, color=RED)
            label(ax, gx + 10.5, gy - 3.5, "model fit on A, scored on B", ha="center", fontsize=6)
    # ---------- c ----------
    panel_letter(ax, 100, 61, "c")
    label(ax, 106, 58.5, "Datasets", fontsize=7)
    card(ax, 100, 8, 79, 46, fc="white")   # table
    hdr = ["system", "perturbation", "independent unit", "null control"]
    xs = [102, 126, 145, 164]
    for x, t in zip(xs, hdr):
        label(ax, x, 49.5, t, fontsize=6, fontweight="bold", color="#34495e")
    ax.plot([101.5, 178], [47, 47], color="#d5d8dc", lw=0.6)
    rows = [("primary CD4 T\n(4 donors)", "CRISPRi,\ngenome-wide", "guides; donors", "unexpressed\ntargets;\nNT guides"),
            ("K562 / RPE1", "CRISPRi,\nessential", "two K562\nscreens", "unexpressed\ntargets*;\nNT guides"),
            ("Jurkat / HepG2", "CRISPRi,\nessential", "none", "unexpressed\ntargets*;\nNT guides"),
            ("H1 hESC\n(VCC 2025)", "CRISPRi,\n150 targets", "none", "NT\nconstructs")]
    for i, r in enumerate(rows):
        y = 41.5 - i * 8.6
        if i % 2 == 0:
            ax.add_patch(Rectangle((101.5, y - 4.1), 76.5, 8.2, fc=LIGHT, ec="none", zorder=0.5))
        for x, t in zip(xs, r):
            col = BLUE if (x == 145 and t != "none") else ("#b03a2e" if t == "none" else "#212f3c")
            label(ax, x, y, t, fontsize=5.8, color=col)
    save(fig, out, "fig1_schematic")


def fig3a(out):
    fig, ax = canvas(89, 52); draw_fig3a(ax); save(fig, out, "fig3a_trbc")


def draw_fig3a(ax, letter=True):
    W, H = 89, 52; rng = np.random.default_rng(7)
    if letter:
        panel_letter(ax, 1.5, H - 1, "a")
    label(ax, 7, H - 3.5, "TCR β constant region as a clone marker", fontsize=7)
    # locus
    y = 38
    ax.plot([4, 59.5], [y, y], color="#aab7b8", lw=1.0, zorder=1)
    segs = [(6, "V", "#aab7b8", 7), (16, "D-J1", "#5d6d7e", 8), (26, "C1", ORANGE, 8), (40, "D-J2", "#5d6d7e", 8), (50, "C2", "#f5b041", 8)]
    for x, t, c, w in segs:
        ax.add_patch(Rectangle((x, y - 2.2), w, 4.4, fc=c, ec="white", lw=0.5, zorder=2)); label(ax, x + w / 2, y, t, ha="center", fontsize=5.8, color="white", fontweight="bold")
    label(ax, 61, y, "TCR β locus", fontsize=6, color=GREY)
    # clones
    fr = [0.9, 0.15, 0.6, 0.05, 0.75]
    for k, f in enumerate(fr):
        cx = 10 + k * 15; cy = 23
        n = 5
        for j in range(n):
            a = 2 * np.pi * j / n; c = ORANGE if rng.random() < f else "#f8c471"
            cell(ax, cx + 2.6 * np.cos(a), cy + 2.6 * np.sin(a), c, r=1.25)
        label(ax, cx, cy - 5.2, "clone %d" % (k + 1), ha="center", fontsize=5.5, color=GREY)
    ax.add_patch(Circle((8, 9), 1.2, fc=ORANGE, ec="white")); label(ax, 10, 9, "TRBC1-high cell", fontsize=5.8)
    ax.add_patch(Circle((33, 9), 1.2, fc="#f8c471", ec="white")); label(ax, 35, 9, "TRBC2-high cell", fontsize=5.8)


def fig5a(out):
    fig, ax = canvas(183, 50); draw_fig5a(ax); save(fig, out, "fig5a_units")


def draw_fig5a(ax, letter=True):
    W, H = 183, 50; rng = np.random.default_rng(11)
    if letter:
        panel_letter(ax, 1.5, H - 1, "a")
    label(ax, 7, H - 3.5, "Shared and separate units", fontsize=7)
    # shared: one clone, two states (left two thirds of the page; panel b sits on the right)
    card(ax, 4, 4, 56, 36, fc="white", ec=RED, lw=0.9)
    label(ax, 32, 36.5, "one clone, two states (CD4)", ha="center", fontsize=6.5, fontweight="bold", color=RED)
    cluster(ax, 16, 19, G1[0], 6, rng)
    arrow(ax, 20, 22, 32, 27, c=G1[0], lw=0.8); arrow(ax, 20, 16, 32, 11, c=G1[0], lw=0.8)
    for yy, nm in ((27, "rest"), (11, "stimulated")):
        for j in range(3):
            cell(ax, 35 + j * 2.6, yy, G1[0], r=1.1)
        label(ax, 44, yy, nm, fontsize=6)
    # separate screens
    card(ax, 64, 4, 55, 36, fc="white", ec=BLUE, lw=0.9)
    label(ax, 91.5, 36.5, "separate screens (K562 vs RPE1)", ha="center", fontsize=6.5, fontweight="bold", color=BLUE)
    cluster(ax, 77, 24, G1[1], 5, rng); cluster(ax, 77, 13, G2[1], 5, rng)
    label(ax, 83, 24, "line 1 clones", fontsize=6); label(ax, 83, 13, "line 2 clones", fontsize=6)


def fig6(out):
    W, H = 183, 64
    fig, ax = canvas(W, H)
    label(ax, 4, H - 3.5, "Can random halves be trusted?", fontsize=8, fontweight="bold")
    def box(x, y, w, h, head, body="", ec="#aab7b8", c="#212f3c"):
        card(ax, x, y, w, h, fc="white", ec=ec, lw=0.8, r=1.5, z=1)
        if body:
            label(ax, x + 3, y + h - 4, head, fontsize=6.6, color=c, fontweight="bold")
            label(ax, x + 3, y + h - 7, body, fontsize=6, color="#34495e", va="top", linespacing=1.45)
        else:
            label(ax, x + w / 2, y + h / 2, head, ha="center", fontsize=6.6, color=c, fontweight="bold", linespacing=1.4)
    # question
    box(4, 32, 40, 16, "Do the two halves share\nclones or wells?")
    # no -> replicates
    box(58, 42, 50, 14, "No: use them as replicates", "different guides, donors or\nreplicate plates", ec=BLUE, c=BLUE)
    # yes -> how large
    box(58, 18, 50, 18, "Yes: expect inflation", "random halves of the same cells,\nor cells from the same well", ec=RED, c=RED)
    arrow(ax, 44.5, 43, 57.5, 49, lw=0.8); arrow(ax, 44.5, 37, 57.5, 27, lw=0.8)
    label(ax, 50, 48.5, "no", fontsize=6, color=BLUE, ha="center"); label(ax, 50, 29, "yes", fontsize=6, color=RED, ha="center")
    # size of the effect
    box(122, 30, 57, 22, "Likely large", "primary cells, strong bottlenecks,\nseveral cell lines per well,\nthe same clones in several conditions", ec=RED, c=RED)
    box(122, 6, 57, 20, "Likely small", "cell lines screened separately,\nstrong effects", ec="#aab7b8")
    arrow(ax, 108.5, 31, 121.5, 40, lw=0.8); arrow(ax, 108.5, 23, 121.5, 17, lw=0.8)
    # always
    card(ax, 4, 4, 104, 10, fc="#eaf2f8", ec=BLUE, lw=0.8, r=1.5, z=1)
    label(ax, 7, 9, "Always include perturbations that should do nothing, compare them\nwith their own mean, and do not take them from growth-selected libraries.",
          fontsize=6, color=BLUE, linespacing=1.35)
    save(fig, out, "fig6_decision")


def save(fig, out, name):
    import matplotlib.text as mtext
    for t in fig.findobj(mtext.Text):          # same type floor as plot_main_figs.MIN_PT
        if t.get_text().strip() and t.get_fontsize() < 5.6:
            t.set_fontsize(5.6)
    fig.savefig("%s/%s.pdf" % (out, name)); fig.savefig("%s/%s.png" % (out, name), dpi=300); plt.close(fig)


def overview(out):
    import matplotlib.image as mpimg
    names = ["fig1_schematic", "fig3a_trbc", "fig5a_units", "fig6_decision"]
    ims = [mpimg.imread("%s/%s.png" % (out, n)) for n in names]
    fig = plt.figure(figsize=(12, 15))
    hs = [im.shape[0] / im.shape[1] for im in ims]; tot = sum(hs) + 0.15 * len(hs); y = 1.0
    for im, h in zip(ims, hs):
        hh = h / tot; y -= hh + 0.15 / tot
        a = fig.add_axes([0.02, y, 0.96 if im.shape[1] > 1500 else 0.48, hh]); a.imshow(im); a.axis("off")
    fig.savefig("%s/schematics_overview.png" % out, dpi=110); plt.close(fig)


if __name__ == "__main__":
    out = sys.argv[1]
    fig1(out); fig3a(out); fig5a(out); print("ok")
