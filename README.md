# random-halves-not-replicates

Analysis code for the preprint *Random cell halves are not independent replicates in pooled perturbation screens*
(Xuanrui Huang, National University of Singapore).

Archived at Zenodo: version 1.0.0, https://doi.org/10.5281/zenodo.23476350; all versions,
https://doi.org/10.5281/zenodo.23476349.

The scripts are the ones that produced the numbers in the paper. For release, only paths and imports were changed:
input and output locations are command-line arguments or come from one configuration file, and the scripts import each
other from the `code/` folder. No computation was changed.

## Quick start: figures from the included result files

`data/` holds the small result files (2.3 MB) that the figure and table scripts read, so the figures can be redrawn
without the raw data:

```bash
pip install -r requirements.txt
mkdir -p figures/supp
python code/plot_schematics.py figures
python code/plot_main_figs.py data figures
python code/make_supplement.py data paper/supplementary_v1.template.md figures/supplementary_v1.md figures/supp
```

This writes Figures 1–5 (`fig1_schematic` … `fig5_contexts`, PDF and PNG), Supplementary Figures S1–S3 and the
Supplementary Information text with Tables S1–S6 (`paper/supplementary_v1.template.md` filled in). Figures 1–5, Supplementary Figures S1–S3 and the Supplementary Information text
drawn this way are identical to the published ones (pixel for pixel and line for line; checked with Python 3.14,
NumPy 2.5, Matplotlib 3.11 on Windows). Arial is used if installed, otherwise DejaVu Sans, which changes the glyphs but
not the content. `paper/` holds the manuscript text and the Supplementary Information template the tables are filled
into; `code/build_preprint_sn.py` typesets the preprint PDF from them with pandoc and pdflatex, using the Springer
Nature LaTeX template (journal article template package, December 2024), which is not redistributed here: place its
files in `paper/sn_template/` and the generated `supplementary_v1.md` in `paper/` first.
`data/README.md` lists, for every file, the script and stage that produce it.

## Full analysis

### Configuration

Copy `config.example.env` to `config.env` and set the locations. One file serves both the shell pipelines and the
Python scripts (`code/config.py`); environment variables of the same name take precedence.

| Key | Meaning |
|---|---|
| `DATA_ROOT` | raw public data, laid out as below |
| `TAHOE_ROOT` | Tahoe-100M and the PIE_tahoe100m files (default `${DATA_ROOT}/tahoe`) |
| `WORK` | intermediate arrays and per-analysis outputs (`${WORK}/results/...`) |
| `FIGDATA` | where the `collect` stage gathers the small result files (default `${WORK}/figdata`) |
| `WORKERS` | processes for the extraction scripts (the random halves do not depend on it) |
| `PYTHON`, `CELLEVAL_PYTHON` | interpreters of the analysis environment and of the cell-eval environment |

### Data

All data are public. `bash pipelines/download.sh` builds the manifest (`code/make_manifests.py`) and downloads it with
size and checksum checks: 3,401 files, about 1.04 TB.

| Dataset | Source | Files under `DATA_ROOT` | Used for |
|---|---|---|---|
| Genome-scale CRISPRi Perturb-seq in primary CD4+ T cells (Zhu et al. 2026) | https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data | `cd4_sc/{D1,D2}_{Rest,Stim48hr}.assigned_guide.h5ad`, `cd4/GWCD4i.pseudobulk_merged.h5ad` | null perturbations, clone marker, depth, ceilings, predictor; donor halves |
| Tahoe-100M (Zhang et al. 2026) | https://huggingface.co/datasets/tahoebio/Tahoe-100M | `tahoe/Tahoe-100M/data/train-*.parquet` (3,388 shards), `tahoe/Tahoe-100M/metadata/{sample,cell_line,gene}_metadata.parquet` | well effects, replicate plates, ceilings, predictor |
| Arc Virtual Cell Challenge preprocessing of Tahoe-100M | https://huggingface.co/datasets/arcinstitute/PIE_tahoe100m | `tahoe/PIE_tahoe100m/preprocessed/{meta.json,ctrl_means.npy}` | gene axis and DMSO means for the Tahoe extraction |
| Replogle et al. 2022: K562 essential (day 6), K562 genome-wide part 1 (day 8), RPE1 essential | figshare+ 20029387 | `replogle/{K562_essential,K562_gwps,rpe1}_raw_singlecell_01.h5ad` | independent screens, template check, Supplementary Table S4; non-targeting null (Fig. 2d, Table S1) |
| Nadig et al. 2025, Jurkat and HepG2 essential | GEO GSE264667 (by hand) | `nadig/jurkat.h5ad`, `nadig/hepg2.h5ad` (raw single-cell files, renamed) | pass test (Supplementary Table S4); non-targeting null (Fig. 2d, Table S1) |
| VCC 2025 H1 hESC training data (Arc Virtual Cell Challenge) | Arc Virtual Cell Atlas, Google Cloud public dataset https://console.cloud.google.com/marketplace/product/bigquery-public-data/arc-institute (by hand) | `vcc2025/adata_Training.h5ad` | non-targeting constructs as null perturbations (Fig. 2d, Supplementary Table S1) |
| Human Protein Atlas, cell-line RNA | https://www.proteinatlas.org (by hand) | `hpa/rna_celline.tsv` | Tahoe line matching in `tahoe_plate_recompute.py` |

### Environments

- **Analysis:** Python ≥ 3.10 with NumPy, SciPy, h5py and PyArrow; Matplotlib for the figures (`requirements.txt`). The
  paper used Python 3.11, NumPy 2.4, SciPy 1.17, h5py 3.16, PyArrow 20 and Matplotlib 3.10.
- **cell-eval:** version 0.8.2 in its own environment (`requirements-celleval.txt`), used only by
  `run_celleval_ceiling.py`.

### Order of runs

`bash pipelines/run_all.sh` runs every stage in order; `bash pipelines/run_all.sh <stage> ...` runs selected stages.
Each stage reads only the outputs of earlier stages.

| Stage | What it does | Scripts |
|---|---|---|
| `extract_cd4` | per-cell aggregates of the four CD4 single-cell files: random halves, guide halves, 8-way per-cell split with TRBC1 / TRBC2 / TRAC counts, non-targeting guides | `cd4_sc_extract.py`, `cd4_sc_extract_guide.py`, `cd4_sc_extract_cells.py`, `cd4_ntc_extract.py` |
| `extract_tahoe` | per line × sample × random half aggregates over all Tahoe shards (map, then reduce) | `tahoe_extract.py` |
| `cd4` | single-cell analyses for donors 1 and 2; non-targeting-guide null; TRBC across states and effective-clone bound | `cd4_inert2.py`, `robust_cd4_tahoe.py cd4`, `cd4_kd_matched.py`, `readout_inert_line.py cd4`, `ceiling_compare.py cd4`, `clone_tests_cd4.py`, `ci_bootstrap.py cd4`, `export_fig_data.py`, `verify_headline.py cd4`, `cd4_ntc_null.py`, `cd4_ntc_ratio.py`, `cd4_trbc_state_keff.py` |
| `cd4_pseudobulk` | released four-donor pseudobulk: donor halves; same vs different guide across donors | `cd4_donor_split.py`, `cd4_guide_donor.py` |
| `celleval` | cell-eval 0.8.2 itself on a donor-1 resting subset | `cd4_celleval_subset.py` → `run_celleval_ceiling.py` |
| `tahoe` | plates sharing conditions, replicate plates 6 / 14, DMSO wells, ceilings, predictor | `tahoe_plate_overlap.py`, `tahoe_plate_recompute.py`, `tahoe_dmso_null.py`, `ceiling_compare.py tahoe`, `robust_cd4_tahoe.py tahoe`, `ci_bootstrap.py tahoe`, `verify_headline.py tahoe` |
| `replogle` | K562 / RPE1 independent screens, template check, pass test for K562–RPE1 and Jurkat–HepG2, non-targeting constructs as null perturbations in the five cell-line screens and the H1 hESC data | `replogle_indep_readout.py` → `k562_inert_screens.py`, `k562_inert.py`, `readout_inert_line.py h5ad`, `cl_ntc_null.py` |
| `sim` | three-term unit model (simulation only, about 4 minutes) | `clone_model_sim.py` |
| `collect` | copy the small result files into `FIGDATA` under the names in `data/` | — |
| `figures` | Figures 1–5, Supplementary Figures S1–S3, Supplementary Information with Tables S1–S6 | `plot_schematics.py`, `plot_main_figs.py`, `make_supplement.py` |

Shared modules, imported by the scripts above: `cd4_m3.py`, `cd4_readout.py`, `tahoe_arms.py`, `tahoe_readout.py`,
`tahoe_hpa_check.py`, `tahoe_tone_check.py` (gene sets, ridge readout, line matching) and `explore_V_day6_day8.py`
(readers for the Replogle files, per-GEM-group control correction). Their own command-line entry points belong to
earlier exploratory analyses and are not part of the paper.

Determinism: random halves, cell assignments and all counts are fixed by seeds and do not depend on `WORKERS`. The
extraction scripts accumulate expression sums in float32 per worker, so a different number of workers (or, for
Tahoe, a different batching of shards) changes the sums in the last bits (relative differences up to about 1e-6).
Everything downstream agrees with the published values at their printed precision; the full reproduction is described
in `REPRODUCTION.md`.

Compute: the extraction stages stream each raw file once (the CD4 files are 140–160 GB each, the Tahoe shards 0.34 TB
in total) and were run with 16–48 processes on a large-memory server. The later stages work on the extracted
aggregates.

## From data to numbers

Result file names refer to `data/` (see `data/README.md` for the producing output).

### Figures

| Manuscript item | Script(s) | Result files |
|---|---|---|
| Fig. 1 (schematic), Fig. 3a, Fig. 5a | `plot_schematics.py` | — |
| Fig. 2a: unexpressed targets positive under random vs guide halves, CD4 donors 1–2 | `export_fig_data.py` → `plot_main_figs.py` | `cd4_D1.npz`, `cd4_D2.npz` |
| Fig. 2b: guide and donor halves agree only by chance | `export_fig_data.py`, `cd4_donor_split.py` (null baseline, variant C) | `cd4_D*.npz`, `donor_split.json` |
| Fig. 2c: Tahoe DMSO well difference | `tahoe_dmso_null.py` | `dmso.json` |
| Fig. 2d: non-targeting random-half signal relative to targeting perturbations, five cell-line screens, H1 hESC and CD4 | `cl_ntc_null.py`, `cd4_ntc_ratio.py` | `cl_ntc_replogle.json`, `cl_ntc_nadig.json`, `cl_ntc_h1.json`, `cd4_ntc_ratio.json` |
| Fig. 3b–d: TRBC overdispersion, TRBC deviation between halves, clone skew vs reproducibility | `cd4_sc_extract_cells.py` → `export_fig_data.py` | `cd4_D*.npz` |
| Fig. 3e: depth (quarter vs half of the cells) | `clone_tests_cd4.py` | `clone_D1.json`, `clone_D2.json` |
| Fig. 3f: composition of the random-half signal | `clone_tests_cd4.py`, `cd4_guide_donor.py` | `clone_D*.json`, `guide_donor.json` |
| Fig. 4a: cell-eval 0.8.2 ceilings | `cd4_celleval_subset.py` → `run_celleval_ceiling.py` | `celleval.json` |
| Fig. 4b: reliability labels (split-half rule of Wang et al. 2026) | `ci_bootstrap.py cd4` | `ci_cd4.json` |
| Fig. 4c: Spearman–Brown at chance | analytic, in `plot_main_figs.py` | — |
| Fig. 4d: Tahoe ceilings, same well vs replicate plate | `ceiling_compare.py tahoe` | `ceiling_tahoe.json` |
| Fig. 5b: context-specific share and absolute energies, CD4 | `ci_bootstrap.py cd4` | `ci_cd4.json` |
| Fig. 5c: Tahoe predictor, same well vs other plate | `ci_bootstrap.py tahoe` (uses `verify_headline.cv_readout`) | `ci_tahoe.json` |
| Fig. 5d: CD4 predictor (donors 1–2, strength-matched knockdowns) and K562 with an independent screen | `ci_bootstrap.py cd4`, `cd4_kd_matched.py`, `replogle_indep_readout.py` | `ci_cd4_D2.json`, `kd_D1.json`, `kd_D2.json`, `replogle_indep.json` |
| Fig. 5e: three-term model, fabrication and attenuation | `clone_model_sim.py` | `sim_v6.json` |

### Supplementary material

| Item | Script(s) | Result files |
|---|---|---|
| Table S1: non-targeting guides as null perturbations, CD4, five cell-line screens and H1 hESC | `cd4_ntc_extract.py` → `cd4_ntc_null.py`; `cl_ntc_null.py` → `make_supplement.py` | `ntc_null.json`, `cl_ntc_replogle.json`, `cl_ntc_nadig.json`, `cl_ntc_h1.json` |
| Table S2: same guide vs different guides across donors | `cd4_guide_donor.py` → `make_supplement.py` | `guide_donor.json` |
| Table S3: all cell-eval metrics, corrected and uncorrected | `run_celleval_ceiling.py` → `make_supplement.py` | `celleval.json` |
| Table S4: template-direction pass test (CD4 donors 1–2, K562, RPE1, Jurkat, HepG2) | `readout_inert_line.py` → `make_supplement.py` | `pass_*.json` |
| Table S5: robustness (target set, gene set, halves, k; Tahoe target plate) | `robust_cd4_tahoe.py` → `make_supplement.py` | `robust_D1.json`, `robust_D2.json`, `robust_tahoe.json` |
| Table S6: donor-2 expectations written before analysis | text of the template, not computed | — |
| Fig. S1: three-term model scenarios | `clone_model_sim.py` | `sim_v5.json` (the same values are in `sim_v6.json`) |
| Fig. S2: Tahoe DMSO wells per plate | `tahoe_dmso_null.py` | `dmso.json` |
| Fig. S3: K562 unexpressed targets across two screens; template artefact | `k562_inert_screens.py`, `k562_inert.py` | `k562_inert_screens.json`, `k562_template.json` |

### Numbers in the text

| Number | Script | Result file |
|---|---|---|
| Unexpressed targets positive: random halves 96–99%, guide halves 49–54% (both donors, both states) | `export_fig_data.py`; cross-checks `cd4_inert2.py`, `verify_headline.py cd4`, `ci_bootstrap.py cd4` | `cd4_D*.npz`, `ci_cd4*.json` |
| Donor halves 55% / 53% (null baseline) | `cd4_donor_split.py` | `donor_split.json` |
| Non-targeting guides: random halves 96.6–98.8%, different guides 48.8–49.8% | `cd4_ntc_null.py` | `ntc_null.json` |
| Non-targeting constructs in five cell-line screens and H1 hESC, relative to targeting perturbations; CD4 1.2–1.3 times | `cl_ntc_null.py`, `cd4_ntc_ratio.py` | `cl_ntc_*.json`, `cd4_ntc_ratio.json` |
| TRBC1 deviation correlated between resting and stimulated cells (0.89, 0.94); effective clones per perturbation at most about 7 and 4 | `cd4_trbc_state_keff.py` | `trbc_state_keff.json` |
| Tahoe DMSO wells: same-well difference reproducible, size relative to drugs, share common to lines, cross-plate correlation | `tahoe_dmso_null.py` | `dmso.json` |
| Tahoe: only plates 6 and 14 are full replicates | `tahoe_plate_overlap.py` | printed |
| Tahoe within-well vs cross-plate energy, signal correlation, context-specific share | `tahoe_plate_recompute.py`, `ci_bootstrap.py tahoe` | `tahoe_recompute.json`, `ci_tahoe.json` |
| TRBC overdispersion, correlation between halves, clone-skew Spearman | `clone_tests_cd4.py`, `export_fig_data.py` | `clone_D*.json`, `cd4_D*.npz` |
| Expressed-target energy random vs guide halves (8.9 vs 1.7; 16.3 vs 1.7) | `clone_tests_cd4.py` | `clone_D*.json` |
| Same vs different guide across donors (3.09 vs 1.71; 2.85 vs 1.92), guide-specific ratios 0.80 / 0.48 | `cd4_guide_donor.py` | `guide_donor.json` |
| K562 unexpressed structure recurring across screens (85%) | `replogle_indep_readout.py` → `k562_inert_screens.py` | `k562_inert_screens.json` |
| K562 template artefact (83% → 89%) | `k562_inert.py` | `k562_template.json` |
| CD4 ceilings, random vs guide halves (all, unexpressed, strong knockdowns) | `ceiling_compare.py cd4` | `ceiling_D1_rest.json`, `ceiling_D2_rest.json` |
| Tahoe ceilings | `ceiling_compare.py tahoe` | `ceiling_tahoe.json` |
| Context-specific share (CD4 random / guide / donor halves; Tahoe), absolute energies, reliable shares, predictor with intervals | `ci_bootstrap.py`, `cd4_donor_split.py` | `ci_*.json`, `donor_split.json` |
| Predictor, second implementation and robustness | `robust_cd4_tahoe.py` | `robust_*.json` |
| Predictor, strength-matched knockdowns | `cd4_kd_matched.py` | `kd_D*.json` |
| Predictor, K562 side with an independent K562 screen | `replogle_indep_readout.py` | `replogle_indep.json` |
| Model: fabricated predictor 0.10; attenuation at the CD4 unit-to-signal ratios | `clone_model_sim.py` | `sim_v6.json` |

### Derived quantities

A few numbers in the text are simple ratios of the outputs above:

- **Unit-to-signal energy ratio (CD4 4.0 / 7.5):** unexpressed-target cross-half signal under random halves divided by
  the expressed-target signal under guide halves, both from `clone_D*.json` (`depth`), the same definition as in
  `clone_model_sim.py`.
- **Composition of the random-half signal (Fig. 3f; 66% / 77% / 81% shared within one transduction):**
  rep = guide-half / random-half expressed signal (`clone_D*.json`), guide-specific = rep × the same-guide excess of
  `guide_donor.json`, transduction-level = 1 − rep − guide-specific. Intervals: the same transfer applied to the
  bootstrap of the ratios (`bootstrap_over_targets` in `guide_donor.json`).

## Notes

- Comments and `rules` strings sometimes carry labels such as "Z1 (am)", "findings ca" or "review 2.5". They refer to
  the authors' analysis log and are kept so that each script can be traced to the record it was written for.
- The `collect` stage fills `sim_v5.json` (Supplementary Fig. S1) from the current simulation; the values used by
  Fig. S1 are unchanged between the two versions (largest difference 2 × 10⁻¹⁶).

## License

The code is released under the MIT License (see `LICENSE`). If you use the code or the result files in `data/`, please cite the preprint.
