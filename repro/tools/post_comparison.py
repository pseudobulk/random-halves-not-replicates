"""Post-process repro/comparison.tsv: earlier-check sections, Tahoe notes, text-side notes.
Usage: post_comparison.py <printed_numbers.md> <comparison.tsv> <figure note>"""
import sys

t = open(sys.argv[1], encoding="utf-8").read().split("\n## ")[1:]
L = open(sys.argv[2], encoding="utf-8").read().rstrip("\n").split("\n")
fig_note = sys.argv[3]
for i, l in enumerate(L):
    c = l.split("\t")
    if l.startswith("figure (reproduced data)\tfig3_clones.png"):
        c[6] = "6 of 4.3 million pixels near x 1185-1192, y 1500-1502 (panel e), each by 1/255 (anti-aliasing; the plotted values differ by less than 1e-5)"
    if "cl_ntc" in l and l.startswith("results"):
        c[3] = "max rel " + c[5].split("max rel ")[1].split(";")[0].rstrip(")")
    if l.startswith("results\ttahoe") or l.startswith("extraction\textract/tahoe"):
        c[6] = (c[6] + "; " if c[6] else "") + "Tahoe-100M: not in the paper since v22 (same-well paper)"
    if l.startswith("figure") or l.startswith("supplement"):
        c[6] = (c[6] + "; " if c[6] else "") + fig_note
    L[i] = "\t".join(c)
for sec in t:
    title = sec.split("\n")[0]
    for l in sec.split("\n"):
        if l.startswith("| ") and not l.startswith("| Item"):
            c = [x.strip() for x in l.strip().strip("|").split("|")]
            L.append("\t".join(["earlier check (" + title + ")", c[0], c[2], c[3], "identical at printed precision",
                                "match" if c[2] == c[3] else "MISMATCH", "printed earlier: %s; %s" % (c[1], c[4])]))
open(sys.argv[2], "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
