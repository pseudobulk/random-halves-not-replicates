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
LIGHT_OR = "#f8c471"                            # TRBC2 (light orange); ORANGE is TRBC1
G1 = ["#6c3483", "#9b59b6", "#c39bd3"]          # clones of guide 1 (purples)
G2 = ["#117a65", "#1abc9c", "#76d7c4"]          # clones of guide 2 (teals)
MM = 1 / 25.4


def canvas(w, h):
    fig = plt.figure(figsize=(w * MM, h * MM))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, w); ax.set_ylim(0, h); ax.axis("off"); ax.set_aspect("equal")
    return fig, ax


def card(ax, x, y, w, h, fc="white", ec="#d5d8dc", lw=0.5, r=0, z=0):
    """a plain rectangle with a thin edge (journal style: no rounded or heavy frames)"""
    ax.add_patch(Rectangle((x, y), w, h, fc=fc, ec=ec, lw=min(lw, 0.5), zorder=z))


def rule(ax, x0, x1, y, c):
    """short coloured rule under a heading: the colour carries the meaning (red, shared cells; blue, independent replicates)"""
    ax.plot([x0, x1], [y, y], color=c, lw=1.2, solid_capstyle="butt", zorder=2)


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
    # ---------- b ----------
    panel_letter(ax, 2, 61, "b")
    label(ax, 8, 58.5, "Uses of random halves", fontsize=7)
    uses = [("noise ceilings and\nreliability labels", ""), ("shared and context-\nspecific parts", ""), ("predicting the context-\nspecific part", "")]
    for k, (t, s) in enumerate(uses):
        x0 = 4 + k * 31.5
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
            label(ax, gx + 10.5, gy - 3.5, "energy: half A · half B", ha="center", fontsize=6)
        else:
            card(ax, gx - 1, gy + 5, 10, 10, fc="white", ec=GREY, r=1)
            label(ax, gx + 4, gy + 10, "shared\n(half A)", ha="center", fontsize=5.6)
            arrow(ax, gx + 9.3, gy + 10, gx + 12.3, gy + 10, lw=0.7, ms=5)
            card(ax, gx + 12.7, gy + 5, 10.5, 10, fc="white", ec=RED, r=1)
            label(ax, gx + 17.95, gy + 10, "specific\n(half B)", ha="center", fontsize=5.6, color=RED)
            label(ax, gx + 10.5, gy - 3.5, "scored on new perturbations", ha="center", fontsize=6)
    # ---------- c ----------
    panel_letter(ax, 100, 61, "c")
    label(ax, 106, 58.5, "Datasets", fontsize=7)
    ax.plot([101.5, 178], [52.5, 52.5], color="k", lw=0.6)   # three-line table: top rule
    hdr = ["system", "perturbation", "independent replicate", "null control"]
    xs = [102, 120, 136, 163]
    for x, t in zip(xs, hdr):
        label(ax, x, 49.5, t, fontsize=6, fontweight="bold", color="#34495e")
    ax.plot([101.5, 178], [47, 47], color="k", lw=0.4)
    rows = [("primary CD4 T\n(4 donors)", "CRISPRi,\ngenome-wide", "guides; donors", "unexpressed\ntargets;\nNT guides"),
            ("K562 / RPE1", "CRISPRi,\nessential", "two K562\nscreens", "unexpressed\ntargets*;\nNT guides"),
            ("Jurkat / HepG2", "CRISPRi,\nessential", "none", "unexpressed\ntargets*;\nNT guides"),
            ("H1 hESC\n(VCC 2025)", "CRISPRi,\n150 targets", "none", "NT\nconstructs")]
    for i, r in enumerate(rows):
        y = 41.5 - i * 8.6
        for x, t in zip(xs, r):
            col = BLUE if (x == 145 and t != "none") else ("#b03a2e" if t == "none" else "#212f3c")
            label(ax, x, y, t, fontsize=5.8, color=col)
    ax.plot([101.5, 178], [y - 4.6, y - 4.6], color="k", lw=0.6)   # bottom rule
    save(fig, out, "fig1_schematic")


def fig3a(out):
    fig, ax = canvas(89, 52); draw_fig3a(ax); save(fig, out, "fig3a_trbc")


def draw_fig3a(ax, letter=True):
    W, H = 89, 52; rng = np.random.default_rng(7)
    if letter:
        panel_letter(ax, 1.5, H - 1, "a")
    label(ax, 7, H - 3.5, "TCR β constant region as a clone marker", fontsize=7)
    # locus as a process: the gene before joining, then the two ways a developing T cell joins it.
    # joining to D-J2 cuts out D-J1 and TRBC1; joining to D-J1 leaves TRBC1 as the constant region used (Maciocia 2017)
    cols = [(20, 5, "V", "#aab7b8"), (26.5, 7, "D-J1", "#5d6d7e"), (35, 8, "TRBC1", ORANGE), (44.5, 7, "D-J2", "#5d6d7e"), (53, 8, "TRBC2", LIGHT_OR)]
    def seg(i, y, alpha=1.0):
        x, w, t, c = cols[i]
        ax.add_patch(Rectangle((x, y - 1.8), w, 3.6, fc="white", ec="none", zorder=1.5))   # hides the backbone under faded segments
        ax.add_patch(Rectangle((x, y - 1.8), w, 3.6, fc=c, ec="white", lw=0.5, alpha=alpha, zorder=2))
        faded = alpha < 1
        label(ax, x + w / 2, y, t, ha="center", fontsize=5.6, alpha=0.45 if faded else 1.0,
              color="k" if faded or c == LIGHT_OR else "white", fontweight="bold")
    rows = ((42.5, "TCR β gene", None), (36.5, "joined to D-J1", "→ makes TRBC1"), (30.5, "joined to D-J2", "→ makes TRBC2"))
    for k, (y, left, right) in enumerate(rows):
        label(ax, 4, y, left, fontsize=5.8)
        if right:
            label(ax, 63, y, right, fontsize=5.8)
        if k == 0:
            ax.plot([18.5, 62.5], [y, y], color="#aab7b8", lw=0.8, zorder=1)
            for i in range(5):
                seg(i, y)
        elif k == 1:
            ax.plot([18.5, 62.5], [y, y], color="#aab7b8", lw=0.8, zorder=1)
            for i, a in zip(range(5), (1, 1, 1, 0.3, 0.3)):   # D-J2 and TRBC2 stay in the DNA but are not used
                seg(i, y, a)
        else:
            ax.plot([18.5, 26.5], [y, y], color="#aab7b8", lw=0.8, zorder=1); ax.plot([43, 62.5], [y, y], color="#aab7b8", lw=0.8, zorder=1)
            ax.add_patch(Rectangle((26.5, y - 1.8), 16.5, 3.6, fill=False, ec="#aab7b8", lw=0.5, ls=(0, (2, 1.5)), zorder=1))
            label(ax, 34.75, y, "cut out", ha="center", fontsize=5.6, color="#7f8c8d")
            for i in (0, 3, 4):
                seg(i, y)
    label(ax, 4, 25.3, "the joined DNA is inherited, so all cells of a clone make the same one", fontsize=5.8)
    ring = lambda cx, cy, r: ax.add_patch(Circle((cx, cy), r, fill=False, ec="#aab7b8", lw=0.5, ls=(0, (2, 1.5)), zorder=1))
    # many clones: one cell from each of ten clones, so the TRBC1 share is near the population average
    pat = [ORANGE, LIGHT_OR, ORANGE, ORANGE, LIGHT_OR, ORANGE, LIGHT_OR, ORANGE, ORANGE, LIGHT_OR]
    for i, c in enumerate(pat):
        cx, cy = 9 + (i % 5) * 4.4, 20.2 if i < 5 else 16.2
        cell(ax, cx, cy, c, r=1.1); ring(cx, cy, 1.8)
    label(ax, 17.8, 11.0, "many clones", ha="center", fontsize=6, fontweight="bold")
    label(ax, 17.8, 8.4, "TRBC1 share near average", ha="center", fontsize=5.8)
    # few clones: one large clone dominates, so the share is far from the average
    rng2 = np.random.default_rng(3)
    for (cx, n, c, spread, r) in ((53.5, 7, ORANGE, 3.3, 4.6), (63.5, 2, LIGHT_OR, 1.4, 2.6), (69.5, 1, ORANGE, 0.0, 1.8)):
        cluster(ax, cx, 18.2, c, n, rng2, spread=spread); ring(cx, 18.2, r)
    label(ax, 61.5, 11.0, "few clones", ha="center", fontsize=6, fontweight="bold")
    label(ax, 61.5, 8.4, "TRBC1 share far from average", ha="center", fontsize=5.8)
    # key
    cell(ax, 8, 2.6, ORANGE, r=1.1); label(ax, 10, 2.6, "makes TRBC1", fontsize=5.8)
    cell(ax, 30, 2.6, LIGHT_OR, r=1.1); label(ax, 32, 2.6, "makes TRBC2", fontsize=5.8)
    ring(53, 2.6, 1.8); label(ax, 55.6, 2.6, "one clone", fontsize=5.8)


def fig5a(out):
    fig, ax = canvas(183, 50); draw_fig5a(ax); save(fig, out, "fig5a_units")


def draw_fig5a(ax, letter=True):
    W, H = 183, 50; rng = np.random.default_rng(11)
    if letter:
        panel_letter(ax, 1.5, H - 1, "a")
    label(ax, 7, H - 3.5, "Shared and separate clones", fontsize=7)
    # shared: one clone, two states (left two thirds of the page; panel b sits on the right)
    # one clone in two states and two random halves: the clone sits in all four cell sets, so its deviation reaches
    # the shared part (computed from half A) and the context-specific part (computed from half B)
    label(ax, 32, 36.5, "one clone, two states (CD4)", ha="center", fontsize=6.5, fontweight="bold", color=RED); rule(ax, 18, 46, 34.4, RED)
    label(ax, 19, 30.5, "resting", ha="center", fontsize=6); label(ax, 31, 30.5, "stimulated", ha="center", fontsize=6)
    for yy, nm in ((22.5, "half A"), (11.5, "half B")):
        label(ax, 6.5, yy, nm, fontsize=6)
        for xx in (19, 31):
            for j in range(3):
                cell(ax, xx - 2.6 + j * 2.6, yy, G1[0], r=1.05)
    ax.plot([25, 25], [8, 28], color="#d5d8dc", lw=0.5, zorder=0); ax.plot([13, 37], [17, 17], color="#d5d8dc", lw=0.5, zorder=0)
    arrow(ax, 37.5, 22.5, 41.5, 22.5, c=G1[0], lw=0.8, ms=5); arrow(ax, 37.5, 11.5, 41.5, 11.5, c=G1[0], lw=0.8, ms=5)
    card(ax, 42.2, 18.6, 16.3, 7.8, fc="white", ec=GREY, r=1, z=1)
    label(ax, 50.35, 22.5, "shared part\n(half A)", ha="center", fontsize=5.6)
    card(ax, 42.2, 7.6, 16.3, 7.8, fc="white", ec=RED, r=1, z=1)
    label(ax, 50.35, 11.5, "specific part\n(half B)", ha="center", fontsize=5.6)
    # separate screens
    label(ax, 91.5, 36.5, "separate screens (K562 vs RPE1)", ha="center", fontsize=6.5, fontweight="bold", color=BLUE); rule(ax, 77.5, 105.5, 34.4, BLUE)
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
