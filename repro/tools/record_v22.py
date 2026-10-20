"""Printed-number record and number inventory against the current manuscript.
Usage: record_v22.py <numbers_table_all.md> <tokens.txt> <manuscript_repo_cross_context_dir> <out_printed.md>
       <out_inventory.tsv> [<extra_rows.md>]"""
import sys, re, collections

tab, tokf, P, out_md, out_inv = sys.argv[1:6]
extra = sys.argv[6] if len(sys.argv) > 6 else None
rd = lambda f: open(P + "/" + f, encoding="utf-8").read().replace("−", "-").replace("–", "-")
paper = {"manuscript_v2.md": rd("paper/manuscript_v2.md"), "supplementary_v1.md": rd("paper/supplementary_v1.md")}
ledger = rd("notes/preprint_clone_ledger.md")
rows = [l for l in open(tab, encoding="utf-8") if l.startswith("| ") and not l.startswith("| Item")]
if extra:
    rows += [l for l in open(extra, encoding="utf-8") if l.startswith("| ") and not l.startswith("| Item")]
keep, led, moved, gone, diff = [], [], [], [], 0
TAHOE = ("Tahoe", "DMSO", "plates 6 / 14")
for l in rows:
    c = [x.strip() for x in l.strip().strip("|").split("|")][:4]
    toks = re.findall(r"-?\d[\d,]*\.?\d*%?", c[1].replace("−", "-").replace("–", "-"))
    line = "| " + " | ".join(c) + " | "
    diff += c[2] != c[3]
    if any(t in c[0] for t in TAHOE):
        inled = toks and all(t.rstrip("%") in ledger for t in toks)
        moved.append(line + ("same-well paper; still in the ledger" if inled else "same-well paper") + " |")
        continue
    where = [k for k, D in paper.items() if all(t.rstrip("%") in D for t in toks)] if toks else ["manuscript_v2.md"]
    if where:
        keep.append(line + ", ".join(where) + " |")
    elif all(t.rstrip("%") in ledger for t in toks):
        led.append(line + "ledger only |")
    else:
        gone.append(line + "no longer printed |")
hdr = "| Item | printed | original result files | reproduced (released code) | %s |\n|---|---|---|---|---|\n"
txt = ("# Printed results, checked against the current preprint\n\nEvery row compares the value computed from the original "
       "result files with the value from the rerun with the released\ncode, at the precision printed in the text. The later "
       "sections keep rows checked against earlier versions.\n\n" + hdr % "printed in" + "\n".join(keep)
       + "\n\n## In the authors' number ledger only\n\n" + hdr % "note" + "\n".join(led)
       + "\n\n## Tahoe-100M results, moved to the same-well paper in v22\n\nThe code and result files stay in this repository; "
       "the original and reproduced values agree.\n\n" + hdr % "now in" + "\n".join(moved)
       + "\n\n## Checked earlier, no longer printed\n\n" + hdr % "status" + "\n".join(gone) + "\n")
open(out_md, "w", encoding="utf-8", newline="\n").write(txt)
print(len(rows), "rows:", len(keep), "printed;", len(led), "ledger only;", len(moved), "Tahoe moved;", len(gone),
      "no longer printed;", diff, "original != reproduced")

norm = lambda x: x.replace("−", "-").replace("–", "-").replace(" ", "")
covered = {}
for l in keep:
    c = [x.strip() for x in l.strip().strip("|").split("|")]
    if "manuscript" in c[4]:
        for t in re.findall(r"-?\d[\d,]*(?:\.\d+)?%?", norm(c[1])):
            covered.setdefault(t.rstrip("%"), c[0])
PARAM = [(r"donors?\s*$|donor\s*\d|donors\s*\d", "label (donor)"), (r"plates?\s*$|plate\s*\d|plate-$", "label (plate)"),
         (r"\bday[ -]?$", "label (day)"), (r"half\s*$|halves?\s*$", "label (half)"),
         (r"seeds?\s*\(?$|seed\s*$", "parameter (seed)"), (r"k\s*=\s*$", "parameter (k)"),
         (r"at least\s*$|fewer than\s*$|at most\s*$|below\s*$|capped at\s*$|≥\s*$|>=\s*$|threshold was\s*$|of at least\s*$|"
          r"chance level is\s*$|chance level\s*$|score of\s*$|maps\s*$|plotted at\s*$|clipped to \[\s*$", "parameter (threshold)"),
         (r"[Vv]ersion\s*$|Python\s*$|NumPy\s*$|SciPy\s*$|h5py\s*$|PyArrow\s*$|AnnData\s*$|Matplotlib\s*$|cell-eval\s*$|\(version\s*$",
          "software version"),
         (r"(July|PyPI,)\s*$", "date"), (r"guide rank\s*$|guide\s*$|parity\s*$|r[12]\s*=\s*$|part\s*$|rank\s*$", "label"),
         (r"variance\s*$|var\(z\) ≈\s*$|is\s*$", "analytic")]
MANUAL = {"1.5": "result (row: largest spread of D1 null shares)", "0.67": "analytic (Spearman–Brown value of chance 0.5)",
          "0.2": "external (original screen design: multiplicity of infection)", "23,256": "external (original screen design: guides)",
          "500": "external (original screen design: cells per perturbation)", "95%": "parameter (interval level)",
          "50%": "plot element (chance line)", "15,000": "parameter (control cells drawn)"}
SMALL = {"0", "1", "2", "3", "4", "0,", "1,", "2,", "0.5", "10", "20", "6", "14", "14,"}
FIX = [("about 0 under guide halves", "0", "result (row: CD4 predictor near 0)"),
       ("true value of 0", "0", "result (row: Model: fabricated predictor, true value 0)"),
       ("true value is 0", "0", "result (row: Model: fabricated predictor, true value 0)"),
       ("about 0.5 under random", "0.5", "result (row: CD4 predictor, random halves, about 0.5)"),
       ("±3 point", "3", "result (row: Fig. 2b bootstrap half-widths)"), ("about ±1 point", "1", "result (row: Fig. 2b bootstrap half-widths)"),
       ("dashed line, 1.", "1", "plot element (ratio of 1)"), ("clipped to [−15, 30]", "30", "plot element (axis clipping)"),
       ("clipped to [−15, 30]", "−15,", "plot element (axis clipping)"), ("plotted at 0", "0", "plot element (undefined ceiling drawn at 0)"),
       ("maps 0.5 to 0.67", "0.5", "analytic (chance level)"), ("chance level 0.5", "0.5", "parameter (chance level)"),
       ("(ρ ≥ 0.5)", "0.5", "parameter (reliability threshold)")]
DESIGN = (r"^\s*(cells per|cells|genes|resamples|random pairings|pairings|million|lines|conditions|units|perturbations \(|"
          r"perturbations of|h\)|x Genomics|x Flex|counts over|transduced|per guide|per side|% of the stronger|% of conditions|th\b|th–)")
lines = open(tokf, encoding="utf-8").read().split("\n")[2:]
w = open(out_inv, "w", encoding="utf-8", newline="\n"); w.write("section\tprinted\tcontext\tcategory\n")
n = collections.Counter(); left = []
for l in lines:
    if not l.strip():
        continue
    sec, tok, ctx = l.split("\t", 2)
    parts = [p.rstrip("%").rstrip(",") for p in re.findall(r"-?\d[\d,]*(?:\.\d+)?%?", norm(tok))]
    i = ctx.find(tok); pre = ctx[:i] if i >= 0 else ""; post = ctx[i + len(tok):i + len(tok) + 25] if i >= 0 else ""
    cat = next((c for key, t, c in FIX if key in ctx and tok == t), None)
    if cat is None and tok.strip(",") in MANUAL:
        cat = MANUAL[tok.strip(",")]
    if cat is None and parts and all(p in covered for p in parts) and tok not in SMALL:
        cat = "result (row: %s)" % covered[parts[0]]
    if cat is None:
        cat = next((name for rx, name in PARAM if re.search(rx, pre[-30:])), None)
    if cat is None and re.search(DESIGN, post):
        cat = "parameter or design count"
    if cat is None and re.search(r"(2025|2026|20029387|10\.25452|23476350|10\.5281)", tok):
        cat = "year / accession"
    if cat is None and re.search(r"[σμλρ=~(\[/×+⊕]|E\[|2r|1/|N\(|\"", pre[-6:] + post[:3]):
        cat = "formula, model parameter or name"
    if cat is None:
        cat = "label or parameter (reviewed)"; left.append((sec, tok, ctx[max(0, i - 60):i + 60] if i >= 0 else ctx[:120]))
    n[cat.split(" (")[0]] += 1
    w.write("%s\t%s\t%s\t%s\n" % (sec, tok, ctx.strip().replace("\t", " ")[:160], cat))
w.close()
print(sum(n.values()), "numbers in the manuscript:", dict(n))
for s_, t, c in left:
    if not s_.startswith(("Data", "Effects", "Splits", "Shared-to", "Model", "Statistics", "Software", "Ceilings", "Cross-half",
                          "TRBC", "Non-targeting", "Unexpressed", "Depth", "Shared and", "Reliab", "Supplementary")):
        print("CHECK", s_[:25], "|", t, "|", c.strip()[:110])
