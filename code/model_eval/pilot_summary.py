"""Model pilot: one table and one figure of every model's score against same-clone truth (A) and independent-replicate
truth (B), from the score files written by the pilot scripts. Models whose fit check showed that they did not learn
per-guide effects (CPA default and adapted, biolord tutorial settings) are listed in the table but not drawn.
Usage: pilot_summary.py <out_dir> -> <out_dir>/pilot_summary.md, pilot_summary.png"""
import os, sys, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RED, BLUE, GREY = "#c0392b", "#2c3e50", "#95a5a6"


def get(o, *ks):
    for k in ks:
        o = o[k]
    return o


def main():
    out = sys.argv[1]
    J = lambda p: json.load(open(os.path.join(out, p)))
    b = J("baselines.json")["models"]
    rows = [  # (label, group, scores dict with inert/strong, drawn, note)
        ("Context mean", "baseline", b["mean"], True, "no perturbation information"),
        ("Copy resting effect", "baseline", b["copy"], True, ""),
        ("Additive", "baseline", b["additive"], True, ""),
        ("Linear (ridge)", "baseline", b["linear"], True, ""),
        ("Arc perturb_mean", "Arc", J("arc_perturb_mean_scores.json")["scores"], True, ""),
        ("Arc pseudobulk", "Arc", J("arc_pseudobulk_scores.json")["scores"], True, ""),
        ("STATE Transition (HVG)", "Arc", J("state_st_scores.json")["scores"], True, ""),
    ]
    if os.path.exists(os.path.join(out, "state_se_scores.json")):
        rows.append(("State (SE-600M + ST)", "Arc", J("state_se_scores.json")["scores"], True, "SE loaded from safetensors"))
    rows += [
        ("PIE (expression shift)", "Arc", get(J("pie_scores.json"), "scores", "delta_p_pred"), True, "without ncbi_text source"),
        ("PIE (log2 fold change)", "Arc", get(J("pie_scores.json"), "scores", "lfc_pred"), False, ""),
        ("Stack (zero-shot)", "Arc", J("stack/stack_scores.json")["scores"], True, "1,672 of 2,152 genes"),
        ("biolord (adapted)", "other", J("biolord_adapted/biolord_scores.json")["scores"], True, "fit r = %.2f" % J("biolord_adapted/biolord_scores.json")["fit_check"]["resting_fit_pearson"]),
        ("biolord (tutorial)", "other", J("biolord/biolord_scores.json")["scores"], False, "fit r = %.3f; did not learn guides" % J("biolord/biolord_scores.json")["fit_check"]["resting_fit_pearson"]),
        ("CPA (default)", "other", J("cpa/cpa_scores.json")["scores"], False, "fit r = %.3f; did not learn guides" % J("cpa/cpa_scores.json")["fit_check"]["resting_fit_pearson"]),
        ("CPA (adapted)", "other", J("cpa_adapted/cpa_scores.json")["scores"], False, "fit r = %.3f; did not learn guides" % J("cpa_adapted/cpa_scores.json")["fit_check"]["resting_fit_pearson"]),
    ]
    ref = J("baselines.json")["truth_A_vs_B"]
    md = ["| Model | Unexpressed A | Unexpressed B | A − B (SE) | Knockdown A | Knockdown B | A − B (SE) | Note |", "|---|---|---|---|---|---|---|---|"]
    for lab, g, s, drawn, note in rows:
        i, k = s["inert"], s["strong"]
        md.append("| %s | %.3f | %.3f | %.3f (%.3f) | %.3f | %.3f | %.3f (%.3f) | %s |" % (
            lab, i["pearson_A"], i["pearson_B"], i["leak"], i["leak_se"], k["pearson_A"], k["pearson_B"], k["leak"], k["leak_se"], note))
    md.append("")
    md.append("Reference: the two guides' stimulated effects agree at r = %.3f (unexpressed) and %.3f (knockdown)." % (
        ref["inert"]["pearson_B"], ref["strong"]["pearson_B"]))
    open(os.path.join(out, "pilot_summary.md"), "w").write("\n".join(md) + "\n")
    drawn = [r for r in rows if r[3]]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 0.32 * len(drawn) + 1.2), sharey=True)
    for ax, kind, title in zip(axes, ("inert", "strong"), ("Unexpressed targets", "Strength-matched knockdowns")):
        for y, (lab, g, s, _, _) in enumerate(drawn[::-1]):
            a, bb = s[kind]["pearson_A"], s[kind]["pearson_B"]
            ax.plot([bb, a], [y, y], color="#d5dbdb", lw=2, zorder=1)
            ax.scatter([a], [y], color=RED, s=18, zorder=3); ax.scatter([bb], [y], color=BLUE, s=18, zorder=3)
        ax.axvline(drawn[0][2][kind]["pearson_B"], color=GREY, lw=0.6, ls=":")
        ax.set_title(title, loc="left", fontsize=8); ax.set_xlabel("Pearson r, predicted vs observed effect", fontsize=7)
        ax.tick_params(labelsize=7); ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_yticks(range(len(drawn))); axes[0].set_yticklabels([r[0] for r in drawn[::-1]], fontsize=7)
    axes[1].scatter([], [], color=RED, s=18, label="same clones (A)"); axes[1].scatter([], [], color=BLUE, s=18, label="independent replicate (B)")
    axes[1].legend(fontsize=6.5, frameon=False, loc="lower right")
    fig.tight_layout(); fig.savefig(os.path.join(out, "pilot_summary.png"), dpi=200)
    print("\n".join(md))


if __name__ == "__main__":
    main()
