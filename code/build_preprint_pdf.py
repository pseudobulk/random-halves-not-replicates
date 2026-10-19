"""Assemble the preprint PDF (local; pandoc 3 + xelatex): main text, Methods, references, one page per main figure with
its legend, then the Supplementary Information with Supplementary Figures placed above their legends.
Markdown details converted for LaTeX: <sup>..</sup> -> pandoc superscript; a few Unicode symbols that Arial lacks ->
inline math; the title block -> a centered LaTeX title; horizontal rules dropped.
Usage: build_preprint_pdf.py <paper_dir> <fig_dir> <out.pdf>"""
import os, re, sys, subprocess

FIGS = {1: "fig1_schematic", 2: "fig2_null", 3: "fig3_clones", 4: "fig4_ceilings", 5: "fig5_contexts", 6: "fig6_decision"}
SFIGS = {1: "supp/figS1_model", 2: "supp/figS2_dmso", 3: "supp/figS3_k562"}   # numbered by first citation (v10)

HEADER = r"""
\usepackage{fontspec}
\setmainfont{Arial}
\usepackage{etoolbox}
\AtBeginEnvironment{longtable}{\scriptsize\setstretch{1.0}}
\usepackage{setspace}
\setstretch{1.25}
\usepackage{caption}
\usepackage{float}
\usepackage{xurl}
\setlength{\parskip}{0.4em}
\setlength{\parindent}{0pt}
\usepackage{titlesec}
\titleformat*{\section}{\Large\bfseries\raggedright\hyphenpenalty=10000\exhyphenpenalty=10000}
\titleformat*{\subsection}{\large\bfseries\raggedright\hyphenpenalty=10000\exhyphenpenalty=10000}
\titleformat*{\subsubsection}{\normalsize\bfseries\raggedright\hyphenpenalty=10000\exhyphenpenalty=10000}
"""

TITLE = r"""
```{=latex}
\begin{center}
{\LARGE\bfseries Random cell halves are not independent replicates in pooled perturbation screens\par}
\vspace{1.2em}
{\large Xuanrui Huang\par}
\vspace{0.6em}
National University of Singapore, Singapore\par
Correspondence: huangxuanrui@u.nus.edu\par
\end{center}
\vspace{1em}
```
"""


def common(s):
    s = re.sub(r"<!--.*?-->", "", s, flags=re.S)
    s = re.sub(r"<sup>([^<]+)</sup>", lambda m: "^" + m.group(1).replace(" ", "\\ ") + "^", s)
    s = s.replace("10⁻³", "10^−3^").replace("10⁻²", "10^−2^").replace("10³", "10^3^").replace("10²", "10^2^")
    s = s.replace("W₂", "W~2~").replace("ZᵀZ", "Z^T^Z").replace("⊕", "$\\oplus$").replace("‖", "$\\Vert$")
    s = re.sub(r"\n---\n", "\n", s)
    # variable subscripts written as x_A, γ_B, K_eff, n_lines -> pandoc subscripts (not inside code names)
    s = s.replace("sd_f", "sd~f~").replace("Gσ_a", "Gσ~a~")
    s = re.sub(r"(?<![A-Za-z0-9_\\.`-])([A-Za-zΑ-Ωα-ω])_([A-Za-z0-9]{1,5}(?:,[A-Za-z0-9]{1,3})?)(?![A-Za-z0-9_])", r"\1~\2~", s)
    return s


def figure_block(path, legend):
    return "\n```{=latex}\n\\newpage\n```\n\n![](%s){width=100%%}\n\n%s\n" % (path, legend)


def main():
    pdir, fdir, out = sys.argv[1:4]
    m = open(os.path.join(pdir, "manuscript_v2.md"), encoding="utf-8").read()
    s = open(os.path.join(pdir, "supplementary_v1.md"), encoding="utf-8").read()
    # main: drop title block (everything before ## Abstract), split off figure legends
    m = m[m.index("## Abstract"):]
    i, j = m.index("## Figure legends"), m.index("## Ethics")
    legends, m = m[i:j], m[:i] + m[j:]
    r = m.index("## References")
    e = m.index("## Ethics")
    back = m[:e] + "\n```{=latex}\n\\small\\setstretch{1.05}\n```\n\n" + m[e:]   # end matter and references in small type
    figs = []
    for n, stem in FIGS.items():
        mm = re.search(r"(\*\*Figure %d \|.*?)(?=\n\n\*\*Figure %d \||\Z)" % (n, n + 1), legends, flags=re.S)
        figs.append(figure_block(os.path.join(fdir, stem + ".pdf").replace("\\", "/"), mm.group(1).strip()))
    # supplementary: drop its own title lines, add figure images under their headings
    s = s[s.index("Contents:"):]
    for n, stem in SFIGS.items():
        h = re.search(r"## Supplementary Figure S%d \|[^\n]*\n" % n, s)
        img = "\n![](%s){width=100%%}\n" % os.path.join(fdir, stem + ".pdf").replace("\\", "/")
        s = s[:h.end()] + img + s[h.end():]
        s = s.replace(h.group(0), "```{=latex}\n\\newpage\n```\n\n" + h.group(0), 1)
    s = re.sub(r"## Supplementary Table (S[2-5]) \|", r"```{=latex}\n\\newpage\n```\n\n## Supplementary Table \1 |", s)
    supp = "```{=latex}\n\\newpage\n```\n\n## Supplementary Information\n\n" + s
    wide = "\n```{=latex}\n\\newgeometry{margin=1.2cm}\n```\n"; narrow = "\n```{=latex}\n\\restoregeometry\n\\normalsize\\setstretch{1.25}\n```\n"
    doc = TITLE + common(back) + wide + "".join(common(f) for f in figs) + narrow + "\n" + common(supp)
    tmp = os.path.join(pdir, "_build_preprint.md"); hdr = os.path.join(pdir, "_build_header.tex")
    open(tmp, "w", encoding="utf-8").write(doc); open(hdr, "w", encoding="utf-8").write(HEADER)
    cmd = ["pandoc", tmp, "-o", out, "--pdf-engine=xelatex", "-H", hdr, "--shift-heading-level-by=-1",
           "-V", "geometry:margin=2.2cm", "-V", "papersize=a4", "-V", "fontsize=10pt", "-V", "colorlinks=true",
           "-f", "markdown+superscript+subscript+raw_tex-implicit_figures", "--columns=40"]
    res = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace")
    print(res.stdout[-2000:], res.stderr[-4000:])
    os.remove(hdr); os.remove(tmp)
    print("built" if res.returncode == 0 else "FAILED", out)


if __name__ == "__main__":
    main()
