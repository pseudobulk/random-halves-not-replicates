"""Methods consistency numbers that live in intermediate results (A100). Usage: methods_numbers.py <tag> <map.json>"""
import sys, json, os
import numpy as np
tag, mp = sys.argv[1], json.load(open(sys.argv[2]))
J = lambda k: json.load(open(mp[k]))
out = {}
# null shares over separate extractions and halvings, donor 1
z = np.load(mp["npz_D1"]); m = z["inert"]
rest = [float((z["cross_random_rest"][m] > 0).mean())]; stim = [float((z["cross_random_stim"][m] > 0).mean())]
ci = J("ci_D1")
for k in ("random_r1", "random_r2", "random_xor"):
    rest.append(ci[k]["inert_frac_pos_rest"][0]); stim.append(ci[k]["inert_frac_pos_stim"][0])
rb = J("robust_D1")["all_common|all_expressed|cell"]; rest.append(rb["inert_frac_pos_rest"]); stim.append(rb["inert_frac_pos_stim"])
c2 = J("inert2_cell_D1")
out["inert2_cell_keys"] = sorted(c2.keys())[:20]
for path in (("A_original", "inert"), ("B_template_removed", "inert"), ("C_null_baseline", "inert")):
    x = c2
    try:
        for p in path: x = x[p]
        if "rest_frac_pos" in x: rest.append(x["rest_frac_pos"]); stim.append(x["stim_frac_pos"])
    except KeyError: pass
v = J("verify_D1"); out["verify_keys"] = sorted(v.keys())
for k, val in v.items():
    if isinstance(val, dict):
        for kk, vv in val.items():
            if "frac_pos" in kk and isinstance(vv, (int, float)):
                (rest if "rest" in kk else stim).append(vv)
out["null_share_rest_min_max"] = [min(rest), max(rest)]; out["null_share_stim_min_max"] = [min(stim), max(stim)]
out["null_share_rest_all"] = rest; out["null_share_stim_all"] = stim
# second implementation readout D1 (cd4_inert2 cell split, template removed, k=20) and main (ci)
def find(x, pre=""):
    if isinstance(x, dict):
        for k, vv in x.items(): yield from find(vv, pre + "/" + k)
    else: yield pre, x
out["inert2_cell_readout"] = {k: v for k, v in find(c2) if "readout" in k and not isinstance(v, list)}
out["inert2_guide_readout"] = {k: v for k, v in find(J("inert2_guide_D1")) if "readout" in k and not isinstance(v, list)}
out["ci_D1_random_r1_readout"] = ci["random_r1"]["readout_k20"][0]
# Tahoe: robust (second impl) and ci/verify (main)
rt = J("robust_tahoe"); out["robust_tahoe_readout"] = {k: v for k, v in find(rt) if "readout" in k and not isinstance(v, list) and "shuffled" not in k}
out["verify_tahoe"] = {k: v for k, v in find(J("verify_tahoe")) if not isinstance(v, list)}
json.dump(out, open("methods_numbers_%s.json" % tag, "w"), indent=1)
print("ok", tag)
