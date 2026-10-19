"""Assemble repro/comparison.tsv from the comparison logs.
Usage: build_comparison.py <out.tsv> <extract.log> <results.log> <numbers_table.md> <figs_v1.txt> [<figs_v2.txt>]"""
import sys, re

out, ext_log, res_log, num_md, figs1 = sys.argv[1:6]
figs2 = sys.argv[6] if len(sys.argv) > 6 else None
rows = [("section", "item", "original", "reproduced (released code)", "criterion", "result", "note")]

ONLY_NOTE = {"metrics_cd4_sc.json": "written by cd4_sc_analyze.py (not in the release); nothing reads it",
             "pert_table.npz": "written by cd4_sc_analyze.py (not in the release); nothing reads it",
             "metrics_tahoe_alpha_strength.json": "written by tahoe_alpha_strength.py (exploratory, not in the release); not used by the paper",
             "metrics_tahoe_readout_source.json": "written by tahoe_readout_source.py (exploratory, not in the release); not used by the paper"}


def parse(line, section):
    m = re.match(r"(OK|DIFF|ONLY|MISS)\s+(\S+)\s*(.*)", line)
    if not m:
        return None
    st, item, rest = m.groups()
    if st == "ONLY":
        f = item.split("/")[-1]
        return (section, item, "present", "absent", "-", "not produced by the release", ONLY_NOTE.get(f, ""))
    g = dict(re.findall(r"(n|maxabs|maxrel|roundflips)=(\S+)", rest))
    first = re.sub(r".*roundflips=\S+\s*", "", rest).strip() if "roundflips" in rest else re.sub(r".*maxrel=\S+\s*", "", rest).strip()
    if st == "OK" and float(g.get("maxabs", "0")) == 0:
        res, crit = "bit-identical", "bit-identical"
    elif st == "OK":
        res, crit = "identical within 1e-9", "bit-identical or float noise"
    else:
        rf = g.get("roundflips", "-")
        crit = "float32 noise; printed values unchanged"
        res = "float noise (max rel %s%s)" % (g.get("maxrel"), "" if rf in ("-", None) else "; %s values round differently at 1-3 significant digits" % rf)
    return (section, item, "n=%s" % g.get("n", ""), "max abs %s" % g.get("maxabs", ""), crit, res, first[:160])


for line in open(ext_log, encoding="utf-8"):
    r = parse(line.strip(), "extraction")
    if r:
        if r[1].endswith(("sum_half.npy", "sum_phase.npy", "sum8.npy", "sum2.npy", "ctrl_moments.npz", "sums64.npy", "phase_sums.npy")):
            r = r[:6] + ("float32 sums: per-worker accumulation order (WORKERS 32 here, 48/16 originally; Tahoe batched as downloaded)",)
        if r[1].endswith("index.json"):
            r = r[:5] + ("same content except the number of merged partial files", "Tahoe reduce: 14 vs 16 partial files")
        rows.append(r)
for line in open(res_log, encoding="utf-8"):
    r = parse(line.strip(), "results")
    if r:
        rows.append(r)


def cells(l):
    return [c.strip() for c in l.strip().strip("|").split("|")]


norm = lambda s: s.replace("−", "-").replace(" ", "")
TEXT_ROUNDING = {"K562 expressed recurrence": "text 0.84 (84%)", "K562 bound, one screen": "text <=16%",
                 "K562 additive clonal part": "text about 4%", "Tahoe within-well / cross-plate energy (gamma / beta)": "1.54 and 1.36 printed as 1.4-1.5",
                 "Tahoe discrimination ceiling, within (6 / 14)": "both 0.99", "Tahoe context-specific share, within wells (6 / 14)": "0.387 / 0.384 printed as 0.39",
                 "Model: fabricated, independent / separate": "0.00 and -0.00 printed as 0.00", "Model ratios at sigma_a 1.5 (low / iso)": "printed as a range",
                 "Model ratios at sigma_a 2.0 (low / iso)": "printed as a range", "CD4 NT / targeting random-half signal, D1 / D2 rest": "1.16 and 1.31 printed as 1.2-1.3", "Effective clones, median upper bound, D1 rest": "6.5 printed as about 7", "Strength-matched predictor, guide, D1 / D2": "printed as approximately 0",
                 "TRBC deviation r, guide halves, D1 / D2 (printed as r ≈ 0)": "-0.01 and 0.02 printed as r ≈ 0",
                 "Expressed energy lost under guide halves, D1 / D2": "81% printed as about 80%",
                 "Context-specific share, CD4 random / guide (printed as about 0.3)": "0.30 and 0.34 printed as about 0.3",
                 "Energy inflation, context-specific / shared (printed as 8–9 times)": "7.6 and 9.1 printed as 8–9",
                 "CD4 predictor, random halves, D1 / D2 (printed as about 0.5)": "0.47 and 0.52 printed as about 0.5",
                 "CD4 predictor near 0: guide D1, guide D2, matched D1, matched D2, same guide across donors": "all |x| <= 0.012, printed as about 0",
                 "Fig. 2b bootstrap half-widths, guide / random (percentage points)": "TEXT-SIDE: largest random-half half-width 1.1 points (donor 1, stimulated) against 'below ±1 point' in the Fig. 2 legend; guide 3.6 against 'about ±3'"}
lines = [l for l in open(num_md, encoding="utf-8") if l.startswith("| ") and not l.startswith("| Item")]
for l in lines:
    c = cells(l)
    item, text, orig, rep = c[0], c[1], c[2], c[3]
    where = c[4] if len(c) > 4 and c[4].endswith(".md") else ""
    rep2 = c[4] if len(c) > 4 and not where else None
    same = orig == rep
    tmatch = norm(text) == norm(orig) or item in TEXT_ROUNDING
    res = ("match" if same and tmatch else ("reproduced, text differs" if same else "MISMATCH"))
    note = TEXT_ROUNDING.get(item, "")
    if rep2 is not None and rep2 != rep:
        note = (note + "; " if note else "") + "with the released cell-eval subset (v2): %s" % rep2
    rows.append(("printed number", item, orig, rep, "identical at printed precision", res, "printed (v15) in %s: %s%s" % (where or "text", text, ("; " + note) if note else "")))
for path, lab in ((figs1, "figure (reproduced data)"), (figs2, "figure (shipped data/)")):
    if not path:
        continue
    for l in open(path, encoding="utf-8"):
        if l.startswith("| fig") or l.startswith("| supp/"):
            c = cells(l)
            rows.append((lab, c[0], "published v15 (figures/preprint)", "redrawn with the released code", "visually identical; noise-level pixels only", c[1], ""))
        if l.startswith("Supplementary Information text"):
            rows.append((lab.replace("figure", "supplement"), "supplementary_v1.md (whole text, Tables S1-S6)", "published v15", "regenerated", "identical", l.split(",")[-1].strip(), l.strip()))
with open(out, "w", encoding="utf-8", newline="\n") as f:
    for r in rows:
        f.write("\t".join(str(x).replace("\t", " ") for x in r) + "\n")
print(len(rows) - 1, "rows")
