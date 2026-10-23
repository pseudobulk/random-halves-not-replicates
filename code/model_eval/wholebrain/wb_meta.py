"""Whole-brain in vivo Perturb-seq (Shi et al. 2026; HF perturbai/wholebrain_crispr_atlas): per-nucleus metadata for the
clone null test, as a numpy file for the A100 (no pandas there). Kept: passes_qc, exactly one guide (num_guides == 1),
neuron_type Glut or GABA. Per nucleus: cell_id, shard (scp_name), mouse (source), guide (guide_call), target
(gene_target; 'Non_target' for the 399 non-targeting guides), class (0 Glut, 1 GABA), and a random half bit (seed 11),
fixed here so that every downstream step uses the same halves.
Usage: wb_meta.py <all_obs.parquet> <out.npz>"""
import sys
import numpy as np, pandas as pd

o = pd.read_parquet(sys.argv[1], columns=["cell_id", "scp_name", "source", "guide_call", "gene_target", "num_guides",
                                          "neuron_type", "passes_qc"])
o = o[o["passes_qc"].astype(bool) & (o["num_guides"] == 1) & o["neuron_type"].isin(["Glut", "GABA"])].reset_index(drop=True)
half = np.random.default_rng(11).integers(0, 2, len(o)).astype(np.int8)
np.savez(sys.argv[2], cell_id=o["cell_id"].astype(str).to_numpy(), shard=o["scp_name"].astype(str).to_numpy(),
         mouse=o["source"].astype(str).to_numpy(), guide=o["guide_call"].astype(str).to_numpy(),
         target=o["gene_target"].astype(str).to_numpy(), cls=(o["neuron_type"] == "GABA").astype(np.int8).to_numpy(), half=half)
print("nuclei", len(o), "Glut", int((o["neuron_type"] == "Glut").sum()), "GABA", int((o["neuron_type"] == "GABA").sum()),
      "guides", o["guide_call"].nunique(), "mice", o["source"].nunique(),
      "non-targeting nuclei", int((o["gene_target"] == "Non_target").sum()))
