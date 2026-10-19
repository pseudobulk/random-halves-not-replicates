"""Typeset the clone-effect preprint with the Springer Nature LaTeX template (paper/sn_template/sn-jnl.cls, journal
article template package, December 2024 version, from the Springer Nature LaTeX author support page; option sn-nature,
the style for Nature Portfolio journals; unnumbered heads as the template offers). The text comes unchanged from
manuscript_v2.md and supplementary_v1.md: pandoc converts the markdown to LaTeX fragments, placed into the template's
structure (title, author and affiliation commands, \abstract, figures with captions, \backmatter with \bmhead
items, a numbered reference list in the template's bibliography layout, then the Supplementary Information).
Compiles with pdflatex inside paper/sn_template/.
Usage: build_preprint_sn.py <paper_dir> <fig_dir> <out.pdf>"""
import os, re, sys, shutil, subprocess

TITLE = "Random cell halves are not independent replicates in pooled perturbation screens"
SHORT = "Random cell halves are not independent replicates"
# five main figures, each floated at the first section that discusses it; Box 1 is marked <!-- BOX START --> / <!-- BOX END -->
FIGS = {1: "fig1_schematic", 2: "fig2_null", 3: "fig3_clones", 4: "fig4_ceilings", 5: "fig5_contexts"}
SFIGS = {1: "supp/figS1_model", 2: "supp/figS2_dmso", 3: "supp/figS3_k562"}   # numbered by first citation
FIG_AT = {1: "## Introduction", 2: "### Null perturbations look reproducible under random halves",
          3: "### The shared state points to clones", 4: "### Random-half ceilings and reliability labels are inflated",
          5: "### Shared units carry cross-context structure that independent units do not reproduce"}


def common(s):
    """markdown details for pandoc: comments dropped, <sup> to pandoc superscripts, a few symbols, variable subscripts"""
    s = re.sub(r"<!--.*?-->", "", s, flags=re.S)
    s = re.sub(r"<sup>([^<]+)</sup>", lambda m: "^" + m.group(1).replace(" ", "\\ ") + "^", s)
    s = s.replace("10⁻³", "10^−3^").replace("10⁻²", "10^−2^").replace("10³", "10^3^").replace("10²", "10^2^")
    s = s.replace("W₂", "W~2~").replace("ZᵀZ", "Z^T^Z").replace("⊕", "$\\oplus$").replace("‖", "$\\Vert$")
    s = re.sub(r"\n---\n", "\n", s)
    # variable subscripts written as x_A, γ_B, K_eff, n_lines -> pandoc subscripts (not inside code names)
    s = s.replace("sd_f", "sd~f~").replace("Gσ_a", "Gσ~a~")
    s = re.sub(r"(?<![A-Za-z0-9_\\.`-])([A-Za-zΑ-Ωα-ω])_([A-Za-z0-9]{1,5}(?:,[A-Za-z0-9]{1,3})?)(?![A-Za-z0-9_])", r"\1~\2~", s)
    return s


def pandoc(md, cols=1000):
    """markdown fragment -> LaTeX fragment"""
    r = subprocess.run(["pandoc", "-f", "markdown+superscript+subscript+raw_tex-implicit_figures", "-t", "latex",
                        "--shift-heading-level-by=-1", "--columns=%d" % cols], input=md, capture_output=True, text=True, encoding="utf-8")
    if r.returncode:
        sys.exit(r.stderr)
    t = r.stdout
    t = re.sub(r"\\(section|subsection|subsubsection)\{", r"\\\1*{", t)          # unnumbered heads, as the template offers
    return t


BOX_OPEN = ("\n```{=latex}\n\\begin{center}\\fbox{\\begin{minipage}{0.94\\textwidth}\\small\\setlength{\\parskip}{0.5em}\n```\n")
BOX_CLOSE = "\n```{=latex}\n\\end{minipage}}\\end{center}\n```\n"


def legends(m):
    out = {}
    for n in FIGS:
        g = re.search(r"\*\*Figure %d \| (.*?)\*\* (.*?)(?=\n\n\*\*Figure %d \||\n---|\Z)" % (n, n + 1), m, flags=re.S)
        out[n] = (g.group(1).strip(), g.group(2).strip())
    return out

KEYWORDS = "Perturb-seq, pooled CRISPR screens, reproducibility, pseudoreplication, benchmarking"

PREAMBLE = r"""\documentclass[pdflatex,sn-nature]{sn-jnl}
\usepackage{graphicx}
\usepackage{amsmath,amssymb,amsfonts}
\usepackage{textcomp}
\usepackage{booktabs,longtable,array,calc}
\usepackage{xcolor}
\usepackage{placeins}
\usepackage{etoolbox}
\DeclareUnicodeCharacter{2212}{\ensuremath{-}}
\DeclareUnicodeCharacter{2265}{\ensuremath{\geq}}
\DeclareUnicodeCharacter{2264}{\ensuremath{\leq}}
\DeclareUnicodeCharacter{2248}{\ensuremath{\approx}}
\DeclareUnicodeCharacter{221A}{\ensuremath{\surd}}
\DeclareUnicodeCharacter{0394}{\ensuremath{\Delta}}
\DeclareUnicodeCharacter{03A3}{\ensuremath{\Sigma}}
\DeclareUnicodeCharacter{03B2}{\ensuremath{\beta}}
\DeclareUnicodeCharacter{03B3}{\ensuremath{\gamma}}
\DeclareUnicodeCharacter{03BB}{\ensuremath{\lambda}}
\DeclareUnicodeCharacter{03BC}{\ensuremath{\mu}}
\DeclareUnicodeCharacter{03C1}{\ensuremath{\rho}}
\DeclareUnicodeCharacter{03C3}{\ensuremath{\sigma}}
\DeclareUnicodeCharacter{2192}{\ensuremath{\rightarrow}}
\providecommand{\tightlist}{\setlength{\itemsep}{0pt}\setlength{\parskip}{0pt}}
\providecommand{\pandocbounded}[1]{#1}
\newcounter{none}
\raggedbottom
\unnumbered
\begin{document}
\title[%(short)s]{%(title)s}
\author*[1]{\fnm{Xuanrui} \sur{Huang}}\email{huangxuanrui@u.nus.edu}
\affil*[1]{\orgname{National University of Singapore}, \orgaddress{\country{Singapore}}}
"""


def fig_width(pdf):
    """170 mm, or less for figures taller than 160 mm, so that figure and legend fit on one page"""
    import pymupdf
    r = pymupdf.open(pdf)[0].rect
    return min(170.0, 160.0 * r.width / r.height)


def figure(n, path, head, body, width=170.0):
    cap = pandoc(common("**%s** %s" % (head, body))).strip()
    # the template's text block is narrower than the 183 mm figures; centre them at up to 170 mm so that their type stays legible
    return ("\\begin{figure}[!tbp]\n\\makebox[\\textwidth][c]{\\includegraphics[width=%.0fmm]{%s}}\n\\caption{%s}\\label{fig:%d}\n\\end{figure}\n"
            % (width, path, cap, n))


def bibliography(refs):
    items = re.findall(r"^\d+\.\s+(.*)$", refs, flags=re.M)
    out = ["\\begin{thebibliography}{%d}" % len(items)]
    for i, t in enumerate(items, 1):
        out.append("\\bibitem{ref%d} %s" % (i, pandoc(common(t)).strip()))
    return "\n".join(out) + "\n\\end{thebibliography}\n"


def main():
    pdir, fdir, out = sys.argv[1:4]
    bdir = os.path.join(pdir, "sn_template")
    m = open(os.path.join(pdir, "manuscript_v2.md"), encoding="utf-8").read()
    s = open(os.path.join(pdir, "supplementary_v1.md"), encoding="utf-8").read()
    m = m.replace("<!-- BOX START -->", BOX_OPEN).replace("<!-- BOX END -->", BOX_CLOSE)
    m = re.sub(r"<!--.*?-->", "", m, flags=re.S)
    leg = legends(m[m.index("## Figure legends"):m.index("## Ethics")])
    abstract = m[m.index("## Abstract") + len("## Abstract"):m.index("## Introduction")].replace("---", "").strip()
    body = m[m.index("## Introduction"):m.index("## Figure legends")]
    end = m[m.index("## Ethics"):m.index("## References")]
    refs = m[m.index("## References") + len("## References"):].strip()
    for n, h in FIG_AT.items():
        i = body.index(h); j = body.index("\n", i)
        body = body[:j + 1] + "\n```{=latex}\n%%%%FIG%d%%%%\n```\n" % n + body[j + 1:]
    tex = pandoc(common(body)).replace("\\section*{", "\\section{").replace("\\subsection*{", "\\subsection{")
    rel = os.path.relpath(fdir, bdir).replace("\\", "/")
    for n, stem in FIGS.items():
        tex = tex.replace("%%%%FIG%d%%%%" % n, figure(n, "%s/%s.pdf" % (rel, stem), *leg[n],
                                                       width=fig_width(os.path.join(fdir, stem + ".pdf"))))
    em = []
    for blk in re.split(r"\n(?=## )", end.strip()):
        head, txt = blk.split("\n", 1)
        head = head.lstrip("# ").strip()
        txt = re.sub(r"\n---\n", "\n", txt).strip()
        em.append("\\bmhead{%s}\n%s" % (head, pandoc(common(txt)).strip()))
    s = s[s.index("Contents:"):]
    for n, stem in SFIGS.items():
        h = re.search(r"## Supplementary Figure S%d \|[^\n]*\n" % n, s)
        s = s[:h.end()] + "\n```{=latex}\n\\begin{center}\\makebox[\\textwidth][c]{\\includegraphics[width=170mm]{%s/%s.pdf}}\\end{center}\n```\n" % (rel, stem) + s[h.end():]
    s = re.sub(r"(## Supplementary (Table|Figure) S)", lambda g: "```{=latex}\n\\clearpage\n```\n\n" + g.group(1), s)
    # a table with panels (a, b, ...) starts each later panel on a new page, and so does the cell-line part of
    # Table S1, so that no table is split across pages
    s = re.sub(r"\n(\*\*[b-z], |\*\*Cell-line screens\.\*\*)", lambda g: "\n```{=latex}\n\\clearpage\n```\n\n" + g.group(1), s)
    supp = pandoc(common(s), cols=40)
    supp = supp.replace("\\begin{longtable}", "{\\fontsize{7}{8.5}\\selectfont\\setlength{\\tabcolsep}{3pt}\\begin{longtable}").replace("\\end{longtable}", "\\end{longtable}}")
    supp = supp.replace("\\section*{", "\\section{").replace("\\subsection*{", "\\subsection{")
    doc = (PREAMBLE % {"short": SHORT, "title": TITLE}
           + "\\abstract{%s}\n" % pandoc(common(abstract)).strip()
           + "\\keywords{%s}\n\\maketitle\n\n" % KEYWORDS
           + tex + "\n\\FloatBarrier\n\\backmatter\n\n" + "\n\n".join(em) + "\n\n" + bibliography(refs)
           + "\n\\clearpage\n\\section{Supplementary Information}\n" + supp + "\n\\end{document}\n")
    name = "_sn_build"
    open(os.path.join(bdir, name + ".tex"), "w", encoding="utf-8").write(doc)
    for _ in range(2):
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", name + ".tex"], cwd=bdir,
                           capture_output=True, text=True, encoding="utf-8", errors="replace")
        if r.returncode:
            print(r.stdout[-3000:]); sys.exit("FAILED")
    shutil.move(os.path.join(bdir, name + ".pdf"), out)
    warn = [l for l in r.stdout.splitlines() if "Overfull" in l or "Missing character" in l or "undefined" in l.lower()]
    print("\n".join(warn[:20]))
    for ext in (".aux", ".log", ".out", ".tex"):
        p = os.path.join(bdir, name + ext)
        if os.path.exists(p):
            os.remove(p)
    print("built", out)


if __name__ == "__main__":
    main()
