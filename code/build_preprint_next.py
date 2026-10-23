"""Typeset the next version of the clone preprint with the Springer Nature LaTeX template (paper/sn_template/sn-jnl.cls,
option sn-nature, unnumbered heads), as build_preprint_sn.py does for the first version. Text from preprint_next.md and
supplementary_next.md; references from references_next_en.md, cited as [@Key] and numbered in order of first citation
(main text, then figure legends, then the Supplementary Information). Seven main figures from figures/next, each floated
at the first section that discusses it. Every table is one float set 183 mm wide (never split across pages); each table
box is measured when compiled, and the build stops if a table is wider than 183 mm or taller than the text block.
Usage: build_preprint_next.py <paper_dir> <fig_dir> <supp_fig_dir> <out.pdf> [--no-author]"""
import os, re, sys, shutil, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_preprint_sn import common, PREAMBLE as PRE_V1, fig_width, BOX_OPEN, BOX_CLOSE

TITLE = "Shared clones make random cell halves look reproducible in pooled screens"
SHORT = "Shared clones make random cell halves look reproducible"
KEYWORDS = "Perturb-seq, pooled CRISPR screens, reproducibility, pseudoreplication, clonal expansion, benchmarking"
FIGS = {n: "fig%d_next" % n for n in range(1, 8)}
FIG_AT = {1: "## Introduction", 2: "### Null perturbations look reproducible under random halves",
          3: "### The halves share clones", 4: "### The size of the shared component depends",
          5: "### Ceilings and reliability labels are inflated", 6: "### Shared clones create cross-context structure",
          7: "### Shared clones raise the scores of transfer models"}
SFIGS = {1: "figS1_model", 2: "figS3_sb", 3: "figS2_k562"}   # numbered by first citation in this version
RR = r">{\raggedright\arraybackslash}"
PREAMBLE = PRE_V1.replace(r"\usepackage{booktabs,longtable,array,calc}", r"\usepackage{booktabs,longtable,array,calc,tabularx}"
                          + "\n\\newsavebox{\\tbox}\n\\DeclareUnicodeCharacter{2295}{\\ensuremath{\\oplus}}\n\\DeclareUnicodeCharacter{2016}{\\ensuremath{\\Vert}}"
                          + "\n\\DeclareUnicodeCharacter{03A9}{\\ensuremath{\\Omega}}\n\\DeclareUnicodeCharacter{2020}{\\textdagger}")


def pandoc(md, cols=1000):
    r = subprocess.run(["pandoc", "-f", "markdown+superscript+subscript+raw_tex-implicit_figures-auto_identifiers", "-t", "latex",
                        "--shift-heading-level-by=-1", "--columns=%d" % cols], input=md, capture_output=True, text=True, encoding="utf-8")
    if r.returncode:
        sys.exit(r.stderr)
    return re.sub(r"\\(section|subsection|subsubsection)\{", r"\\\1*{", r.stdout)


def number_citations(texts, refs):
    """[@A; @B] -> <sup>n,m</sup>, numbered by first appearance over texts in order; returns new texts and the list"""
    order = []
    for t in texts:
        for grp in re.findall(r"\[(@[^\]]+)\]", t):
            for k in re.findall(r"@([A-Za-z0-9\-]+)", grp):
                if k not in order:
                    order.append(k)
    missing = [k for k in order if k not in refs]
    assert not missing, missing
    num = {k: i + 1 for i, k in enumerate(order)}

    def rng(ns):
        ns = sorted(set(ns)); out = []; i = 0
        while i < len(ns):
            j = i
            while j + 1 < len(ns) and ns[j + 1] == ns[j] + 1:
                j += 1
            out.append(str(ns[i]) if j == i else ("%d,%d" % (ns[i], ns[j]) if j == i + 1 else "%d–%d" % (ns[i], ns[j])))
            i = j + 1
        return ",".join(out)

    def rep(m):
        ks = re.findall(r"@([A-Za-z0-9\-]+)", m.group(1))
        return "<sup>%s</sup>" % rng([num[k] for k in ks])
    new = [re.sub(r"\s?\[(@[^\]]+)\]", rep, t) for t in texts]
    new = [re.sub(r"ref\. <sup>(\d+)</sup>", r"ref. \1", t) for t in new]
    unused = [k for k in refs if k not in num]
    return new, [refs[k] for k in order], unused


def bibliography(items):
    out = ["\\begin{thebibliography}{%d}" % len(items)]
    for i, t in enumerate(items, 1):
        out.append("\\bibitem{ref%d} %s" % (i, pandoc(common(t)).strip()))
    return "\n".join(out) + "\n\\end{thebibliography}\n"


def split_panels(body):
    """legend body -> chunks, each starting at a panel letter (**a**, **b**, ...); text before **a** stays with a"""
    idx = [m.start() for m in re.finditer(r"\*\*[a-h]\*\*,", body)]
    if not idx:
        return [body]
    cuts = [0] + idx[1:] + [len(body)]
    return [body[cuts[i]:cuts[i + 1]].strip() for i in range(len(cuts) - 1)]


def figure(n, path, head, body, width, keep=None):
    """figure and legend in one float; with keep set, only the first keep panel descriptions stay with the figure and
    the rest follow in a float of their own ("Figure n (continued)"), so that no float exceeds the page"""
    chunks = split_panels(body)
    k = len(chunks) if keep is None else keep
    cap = pandoc(common("**%s** %s" % (head, " ".join(chunks[:k])))).strip()
    # the float is set in a box and measured (FIG lines in the log), so that the build can move panels out of it
    t = ("\\begin{figure}[!tbp]\n\\sbox{\\tbox}{\\begin{minipage}{\\textwidth}\\makebox[\\textwidth][c]{\\includegraphics[width=%.0fmm]{%s}}\n"
         "\\caption{%s}\\label{fig:%d}\\end{minipage}}\n"
         "\\typeout{FIG %d \\the\\dimexpr\\ht\\tbox+\\dp\\tbox\\relax\\space\\the\\textheight}\\usebox{\\tbox}\n\\end{figure}\n"
         % (width, path, cap, n, n))
    if k < len(chunks):
        rest = pandoc(common(" ".join(chunks[k:]))).strip()
        t += ("\\begin{figure}[!tbp]\n{\\small\\textbf{Figure %d (continued).} %s\\par}\n\\end{figure}\n" % (n, rest))
    return t


def md_cells_to_tex(cells):
    """convert many table cells with one pandoc call"""
    sep = "\n\nQQSEPQQ\n\n"
    tex = pandoc(common(sep.join(c.strip() if c.strip() else "\u00a0" for c in cells)))
    parts = [p.strip() for p in tex.split("QQSEPQQ")]
    assert len(parts) == len(cells), (len(parts), len(cells))
    return [p.replace("\u00a0", "") for p in parts]


def table(label, caption_md, md_table):
    """one float, 183 mm wide; l/r columns unless an entry is long (X column of a tabularx); measured at compile time"""
    lines = [l for l in md_table.strip().split("\n") if l.startswith("|") and not re.match(r"^\|[-|\s:]+\|$", l)]
    rows = [[c for c in l.strip().strip("|").split("|")] for l in lines]
    head, body = rows[0], rows[1:]
    flat = md_cells_to_tex([c for r in rows for c in r])
    k = 0; tex_rows = []
    for r in rows:
        tex_rows.append(flat[k:k + len(r)]); k += len(r)
    ln = lambda c: len(re.sub(r"[*~^`]", "", c.strip()))
    num = re.compile(r"^[−-]?[\d,]*\.?\d+%?$|^–$|^$")
    spec = ""
    for j in range(len(head)):
        col = [r[j].strip() for r in body if j < len(r)]
        longest = max((ln(c) for c in col), default=0)
        kind = "X" if longest > 20 else ("r" if all(num.match(c) for c in col) else "l")
        spec += (RR + "X") if kind == "X" else kind
    hspec = "".join("X" if c == "X" else c for c in spec.replace(RR, ""))
    # headers wrap within their column
    hcells = []
    for j, h in enumerate(tex_rows[0]):
        words, out = h.split(" "), [""]
        width = max(10, max((ln(r[j]) for r in body if j < len(r)), default=10))
        for w in words:
            if hspec[j] != "X" and out[-1] and len(out[-1]) + 1 + len(w) > width:
                out.append(w)
            else:
                out[-1] = (out[-1] + " " + w).strip()
        hcells.append(out[0] if len(out) == 1 else "\\begin{tabular}[b]{@{}l@{}}%s\\end{tabular}" % " \\\\ ".join(out))
    size = "\\fontsize{8}{9.5}\\selectfont" if len(body) <= 20 else "\\fontsize{7}{8.2}\\selectfont"
    if "X" in spec:
        env = ("\\begin{tabularx}{\\linewidth}{@{}%s@{}}" % spec, "\\end{tabularx}")
    else:
        env = ("\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}%s@{}}" % spec, "\\end{tabular*}")
    cap = pandoc(common(caption_md)).strip()
    return "\n".join(["\\begin{table}[!htbp]", "\\renewcommand{\\tablename}{Supplementary Table}", "\\renewcommand{\\thetable}{%s}" % label,
                      "\\sbox{\\tbox}{\\begin{minipage}{183mm}", "\\caption{%s}" % cap,
                      size + "\\setlength{\\tabcolsep}{3pt}\\renewcommand{\\arraystretch}{1.1}", env[0], "\\toprule", " & ".join(hcells) + " \\\\",
                      "\\midrule"] + [" & ".join(r) + " \\\\" for r in tex_rows[1:]] + ["\\bottomrule", env[1], "\\end{minipage}}",
                      "\\typeout{TBL %s \\the\\wd\\tbox\\space\\the\\dimexpr\\ht\\tbox+\\dp\\tbox\\relax\\space\\the\\textheight}" % label,
                      "\\makebox[\\textwidth][c]{\\usebox{\\tbox}}", "\\end{table}"])


def check(log):
    s = open(log, encoding="latin-1").read()
    mm = lambda x: float(x.replace("pt", "")) / 72.27 * 25.4
    tabs = re.findall(r"TBL (\S+) ([\d.]+pt) ([\d.]+pt) ([\d.]+pt)", s)
    for t in tabs:
        print("table %-4s %5.1f mm wide, %5.1f mm tall (text block %5.1f mm)" % (t[0], mm(t[1]), mm(t[2]), mm(t[3])))
    bad = [t[0] for t in tabs if mm(t[1]) > 183.2 or mm(t[2]) > mm(t[3])]
    over = re.findall(r"Overfull \\hbox \(([\d.]+)pt too wide\)", s)
    print("overfull boxes:", over[:10])
    assert tabs and not bad and "Float too large" not in s, bad
    print("all %d tables fit on one page each" % len(tabs))


def supplement(s, rel_sfig, sfdir):
    s = s[s.index("Contents:"):]
    out = []
    for blk in re.split(r"\n(?=## )", s):
        if not blk.startswith("## "):
            out.append(pandoc(common(blk))); continue
        head, txt = blk.split("\n", 1)
        title = head[3:].strip()
        m = re.match(r"Supplementary Figure S(\d) \| (.*)", title)
        if m:
            n = int(m.group(1)); stem = SFIGS[n]
            import pymupdf
            wmm = min(170.0, pymupdf.open(os.path.join(sfdir, stem + ".pdf"))[0].rect.width / 72 * 25.4)
            cap = pandoc(common("**%s.** %s" % (m.group(2), txt.strip()))).strip()
            out.append("\\clearpage\n\\begin{figure}[!htbp]\n\\renewcommand{\\figurename}{Supplementary Figure}\\renewcommand{\\thefigure}{S%d}\n"
                       "\\makebox[\\textwidth][c]{\\includegraphics[width=%.0fmm]{%s/%s.pdf}}\n\\caption{%s}\n\\end{figure}\n" % (n, wmm, rel_sfig, stem, cap))
            continue
        if title.startswith("Supplementary Table"):
            parts = re.split(r"<!-- TABLE (\S+) -->\n", txt)
            pre = parts[0].strip()
            if pre:
                out.append(pandoc(common(pre)))
            for lab, body in zip(parts[1::2], parts[2::2]):
                cap_end = body.index("\n|")
                cap = body[:cap_end].strip()
                cap = re.sub(r"^\*\*Supplementary Table [^|]*\| ", "**", cap)
                rest = body[cap_end:]
                tbl = "\n".join(l for l in rest.strip().split("\n") if l.startswith("|"))
                after = "\n".join(l for l in rest.strip().split("\n") if not l.startswith("|")).strip()
                out.append("\\clearpage\n" + table(lab, cap, tbl))
                if after:
                    out.append(pandoc(common(after)))
            continue
        # a note: heading, text, and tables marked inside it
        out.append("\\clearpage\n\\subsection{%s}\n" % pandoc(common(title)).strip().replace("\\|", "|"))
        parts = re.split(r"<!-- TABLE (\S+) -->\n", txt)
        out.append(pandoc(common(parts[0])))
        for lab, body in zip(parts[1::2], parts[2::2]):
            cap_end = body.index("\n|"); cap = re.sub(r"^\*\*Supplementary Table [^|]*\| ", "**", body[:cap_end].strip())
            rest = body[cap_end:]
            tbl = "\n".join(l for l in rest.strip().split("\n") if l.startswith("|"))
            out.append(table(lab, cap, tbl))
    return "\n".join(out)


def main():
    pdir, fdir, sfdir, out = sys.argv[1:5]
    anonymous = "--no-author" in sys.argv[5:]
    bdir = os.path.join(pdir, "sn_template")
    m = open(os.path.join(pdir, "preprint_next.md"), encoding="utf-8").read()
    s = open(os.path.join(pdir, "supplementary_next.md"), encoding="utf-8").read()
    refs = dict(re.findall(r"^- ([A-Za-z0-9\-]+): (.*)$", open(os.path.join(pdir, "references_next_en.md"), encoding="utf-8").read(), flags=re.M))
    m = re.sub(r"<!--(?! BOX).*?-->", "", m, flags=re.S)
    m = m.replace("<!-- BOX START -->", BOX_OPEN).replace("<!-- BOX END -->", BOX_CLOSE)
    abstract = m[m.index("## Abstract") + len("## Abstract"):m.index("## Introduction")].strip()
    body = m[m.index("## Introduction"):m.index("## Figure legends")]
    legs = m[m.index("## Figure legends"):m.index("## Ethics")]
    end = m[m.index("## Ethics"):]
    (body, legs, s), items, unused = number_citations([body, legs, s], refs)
    print("references:", len(items), "unused:", unused)
    if anonymous:
        end = re.sub(r"## Author contributions\n.*?\n(?=## )", "", end, flags=re.S)
    leg = {}
    for n in FIGS:
        g = re.search(r"\*\*Figure %d \| (.*?)\*\* (.*?)(?=\n\n\*\*Figure %d \||\n\n## |\Z)" % (n, n + 1), legs + "\n\n## ", flags=re.S)
        leg[n] = (g.group(1).strip(), g.group(2).strip())
    for n, h in FIG_AT.items():
        i = body.index(h); j = body.index("\n", i)
        body = body[:j + 1] + "\n```{=latex}\n%%%%FIG%d%%%%\n```\n" % n + body[j + 1:]
    tex0 = pandoc(common(body)).replace("\\section*{", "\\section{").replace("\\subsection*{", "\\subsection{")
    rel = os.path.relpath(fdir, bdir).replace("\\", "/"); srel = os.path.relpath(sfdir, bdir).replace("\\", "/")
    keep = {n: None for n in FIGS}
    em = []
    for blk in re.split(r"\n(?=## )", end.strip()):
        head, txt = blk.split("\n", 1)
        em.append("\\bmhead{%s}\n%s" % (head.lstrip("# ").strip(), pandoc(common(txt.strip())).strip()))
    supp = supplement(s, srel, sfdir).replace("\\section*{", "\\section{").replace("\\subsection*{", "\\subsection{")
    pre = PREAMBLE % {"short": SHORT, "title": TITLE}
    if anonymous:
        pre = re.sub(r"\\author\*.*\n\\affil\*.*\n", "", pre)
    name = "_sn_next"
    abstract_tex = pandoc(common(abstract)).strip(); bib = bibliography(items)
    for it in range(10):   # move legend panels of figures that do not fit on one page into a continuation float
        tex = tex0; starts = {}
        for n, stem in FIGS.items():
            tex = tex.replace("%%%%FIG%d%%%%" % n, "%%%% FIGSTART %d\n" % n + figure(n, "%s/%s.pdf" % (rel, stem), *leg[n],
                              width=fig_width(os.path.join(fdir, stem + ".pdf")), keep=keep[n]))
        doc = (pre + "\\abstract{%s}\n" % abstract_tex + "\\keywords{%s}\n\\maketitle\n\n" % KEYWORDS
               + tex + "\n\\FloatBarrier\n\\backmatter\n\n" + "\n\n".join(em) + "\n\n" + bib
               + "\n\\clearpage\n\\section{Supplementary Information}\n" + supp + "\n\\end{document}\n")
        for i, l in enumerate(doc.split("\n"), 1):
            g = re.match(r"%% FIGSTART (\d)", l)
            if g:
                starts[int(g.group(1))] = i
        open(os.path.join(bdir, name + ".tex"), "w", encoding="utf-8").write(doc)
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", name + ".tex"], cwd=bdir,
                           capture_output=True, text=True, encoding="utf-8", errors="replace")
        if r.returncode:
            print(r.stdout[-3000:]); sys.exit("FAILED")
        log = open(os.path.join(bdir, name + ".log"), encoding="latin-1").read()
        moved = False
        for n, h, th in re.findall(r"FIG (\d) ([\d.]+)pt ([\d.]+)pt", log):
            n = int(n); k = len(split_panels(leg[n][1])) if keep[n] is None else keep[n]
            if float(h) > float(th) - 6 and k > 1:   # leave room for the float separation
                keep[n] = k - 1; moved = True
        print("pass %d, panels kept with the figure:" % (it + 1), {n: k for n, k in keep.items() if k is not None})
        if not moved:
            break
    r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", name + ".tex"], cwd=bdir,
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    check(os.path.join(bdir, name + ".log"))
    warn = [l for l in r.stdout.splitlines() if "Missing character" in l or "undefined" in l.lower()]
    print("\n".join(warn[:20]))
    shutil.move(os.path.join(bdir, name + ".pdf"), out)
    for ext in (".aux", ".log", ".out", ".tex"):
        p = os.path.join(bdir, name + ext)
        if os.path.exists(p):
            os.remove(p)
    print("built", out)


if __name__ == "__main__":
    main()
