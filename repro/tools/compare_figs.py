"""Redraw figures and the Supplementary Information from a figdata folder with the release code and compare with the
published versions. Usage: compare_figs.py <release_repo> <figdata> <out_dir> <published_fig_dir> <published supplementary_v1.md>"""
import sys, os, glob, subprocess, difflib
import matplotlib.image as mi
import numpy as np

repo, figdata, out, pub, supmd = sys.argv[1:6]
os.makedirs(out + "/supp", exist_ok=True)
for cmd in (["plot_schematics.py", out], ["plot_main_figs.py", figdata, out],
            ["make_supplement.py", figdata, os.path.join(repo, "paper", "supplementary_v1.template.md"), out + "/supplementary_v1.md", out + "/supp"]):
    r = subprocess.run([sys.executable, os.path.join(repo, "code", cmd[0])] + cmd[1:], capture_output=True, text=True)
    print(cmd[0], "->", (r.stdout.strip().splitlines() or ["(no output)"])[-1], r.stderr.strip()[-300:])


def cmp(a, b):
    if not os.path.exists(b):
        return "MISSING"
    A, B = mi.imread(a), mi.imread(b)
    if A.shape != B.shape:
        return "shape %s vs %s" % (A.shape, B.shape)
    d = np.abs(A.astype(float) - B.astype(float)).max(-1)
    return "identical" if not d.any() else "differs: %.4f%% of pixels, max channel difference %.3f" % (100 * (d > 0).mean(), d.max())


print("\n| Figure | vs published (figures/preprint) |\n|---|---|")
for f in sorted(glob.glob(pub + "/fig[1-9]*.png")) + sorted(glob.glob(pub + "/supp/figS*.png")):
    rel = os.path.relpath(f, pub)
    print("| %s | %s |" % (rel.replace("\\", "/"), cmp(f, os.path.join(out, rel))))
A = open(supmd, encoding="utf-8").read().replace("\r", "").splitlines()
B = open(out + "/supplementary_v1.md", encoding="utf-8").read().replace("\r", "").splitlines()
diff = list(difflib.unified_diff(A, B, "published", "regenerated", lineterm="", n=0))
print("\nSupplementary Information text: %d lines published, %d regenerated, %s" % (len(A), len(B), "identical" if not diff else "%d differing lines" % sum(1 for l in diff if l[:1] in "+-" and l[:3] not in ("+++", "---"))))
for l in diff[:40]:
    print("   ", l[:200])
