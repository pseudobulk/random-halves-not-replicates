# random-halves-not-replicates

Analysis code for the preprint *Shared clones make random cell halves look reproducible in pooled screens*
(Xuanrui Huang, Department of Biomedical Informatics, Yong Loo Lin School of Medicine, National University of Singapore).

Archived at Zenodo, all versions: https://doi.org/10.5281/zenodo.23476349 (each release has its own DOI there;
version 1.1.0 is the one the first version of the preprint uses; the files of the second version are described in "Second version of the preprint").

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
Supplementary Information text with Tables S1–S5 (`paper/supplementary_v1.template.md` filled in). Figures 1–5, Supplementary Figures S1–S3 and the Supplementary Information text
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
| `figures` | Figures 1–5, Supplementary Figures S1–S3, Supplementary Information with Tables S1–S5 | `plot_schematics.py`, `plot_main_figs.py`, `make_supplement.py` |

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
| Fig. 2c: non-targeting random-half signal relative to targeting perturbations, five cell-line screens, H1 hESC and CD4 | `cl_ntc_null.py`, `cd4_ntc_ratio.py` | `cl_ntc_replogle.json`, `cl_ntc_nadig.json`, `cl_ntc_h1.json`, `cd4_ntc_ratio.json` |
| Fig. 3b–d: TRBC overdispersion, TRBC deviation between halves, clone skew vs reproducibility | `cd4_sc_extract_cells.py` → `export_fig_data.py` | `cd4_D*.npz` |
| Fig. 3e: depth (quarter vs half of the cells) | `clone_tests_cd4.py` | `clone_D1.json`, `clone_D2.json` |
| Fig. 3f: composition of the random-half signal | `clone_tests_cd4.py`, `cd4_guide_donor.py` | `clone_D*.json`, `guide_donor.json` |
| Fig. 4a: cell-eval 0.8.2 ceilings | `cd4_celleval_subset.py` → `run_celleval_ceiling.py` | `celleval.json` |
| Fig. 4b: reliability labels (split-half rule of Wang et al. 2026), donors 1–2 | `ci_bootstrap.py cd4` | `ci_cd4.json`, `ci_cd4_D2.json` |
| Fig. 5b: TRBC1 deviation in resting against stimulated cells | `cd4_sc_extract_cells.py` → `export_fig_data.py` | `cd4_D1.npz` (z / √n per state) |
| Fig. 5c: context-specific share and absolute energies, CD4 | `ci_bootstrap.py cd4` | `ci_cd4.json` |
| Fig. 5d: CD4 predictor (donors 1–2, strength-matched knockdowns) and K562 with an independent screen | `ci_bootstrap.py cd4`, `cd4_kd_matched.py`, `replogle_indep_readout.py` | `ci_cd4_D2.json`, `kd_D1.json`, `kd_D2.json`, `replogle_indep.json` |
| Fig. 5e: three-term model, fabrication and attenuation | `clone_model_sim.py` | `sim_v6.json` |

### Supplementary material

| Item | Script(s) | Result files |
|---|---|---|
| Table S1: non-targeting guides as null perturbations, CD4, five cell-line screens and H1 hESC | `cd4_ntc_extract.py` → `cd4_ntc_null.py`; `cl_ntc_null.py` → `make_supplement.py` | `ntc_null.json`, `cl_ntc_replogle.json`, `cl_ntc_nadig.json`, `cl_ntc_h1.json` |
| Table S2: same guide vs different guides across donors | `cd4_guide_donor.py` → `make_supplement.py` | `guide_donor.json` |
| Table S3: all cell-eval metrics, corrected and uncorrected | `run_celleval_ceiling.py` → `make_supplement.py` | `celleval.json` |
| Table S4: template-direction pass test (CD4 donors 1–2, K562, RPE1, Jurkat, HepG2) | `readout_inert_line.py` → `make_supplement.py` | `pass_*.json` |
| Table S5: robustness (target set, gene set, halves, k) | `robust_cd4_tahoe.py` → `make_supplement.py` | `robust_D1.json`, `robust_D2.json` |
| Fig. S1: three-term model scenarios | `clone_model_sim.py` | `sim_v5.json` (the same values are in `sim_v6.json`) |
| Fig. S2: K562 unexpressed targets across two screens; template artefact | `k562_inert_screens.py`, `k562_inert.py` | `k562_inert_screens.json`, `k562_template.json` |
| Fig. S3: Spearman–Brown at chance | analytic, in `make_supplement.py` | — |

Tahoe-100M: the pipeline stage `tahoe`, its scripts and its result files (`dmso.json`, `tahoe_recompute.json`,
`ceiling_tahoe.json`, `ci_tahoe.json`, `robust_tahoe.json`) remain in the repository and in the reproduction record, but
since preprint v22 the paper reports the pooled CRISPR screens only.

### Numbers in the text

| Number | Script | Result file |
|---|---|---|
| Unexpressed targets positive: random halves 96–99%, guide halves 49–54% (both donors, both states) | `export_fig_data.py`; cross-checks `cd4_inert2.py`, `verify_headline.py cd4`, `ci_bootstrap.py cd4` | `cd4_D*.npz`, `ci_cd4*.json` |
| Donor halves 55% / 53% (null baseline) | `cd4_donor_split.py` | `donor_split.json` |
| Non-targeting guides: random halves 96.6–98.8%, different guides 48.8–49.8% | `cd4_ntc_null.py` | `ntc_null.json` |
| Non-targeting constructs in five cell-line screens and H1 hESC, relative to targeting perturbations; CD4 1.2–1.3 times | `cl_ntc_null.py`, `cd4_ntc_ratio.py` | `cl_ntc_*.json`, `cd4_ntc_ratio.json` |
| TRBC1 deviation correlated between resting and stimulated cells (0.89, 0.94); effective clones per perturbation at most about 7 and 4 | `cd4_trbc_state_keff.py` | `trbc_state_keff.json` |
| TRBC overdispersion, correlation between halves, clone-skew Spearman | `clone_tests_cd4.py`, `export_fig_data.py` | `clone_D*.json`, `cd4_D*.npz` |
| Expressed-target energy random vs guide halves (8.9 vs 1.7; 16.3 vs 1.7) | `clone_tests_cd4.py` | `clone_D*.json` |
| Same vs different guide across donors (3.09 vs 1.71; 2.85 vs 1.92), guide-specific ratios 0.80 / 0.48 | `cd4_guide_donor.py` | `guide_donor.json` |
| K562 unexpressed structure recurring across screens (85%) | `replogle_indep_readout.py` → `k562_inert_screens.py` | `k562_inert_screens.json` |
| K562 template artefact (83% → 89%) | `k562_inert.py` | `k562_template.json` |
| CD4 ceilings, random vs guide halves (all, unexpressed, strong knockdowns) | `ceiling_compare.py cd4` | `ceiling_D1_rest.json`, `ceiling_D2_rest.json` |
| Context-specific share (CD4 random / guide / donor halves), absolute energies, reliable shares, predictor with intervals | `ci_bootstrap.py`, `cd4_donor_split.py` | `ci_*.json`, `donor_split.json` |
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

## Second version of the preprint

The second version adds constant-region (TRBC) halves, clone-barcode halves in brain organoids, one clone-strength
estimator across 16 screens, a depth grid and an evaluation of transfer models. Its code, result files and figures sit
next to those of the first version, which stay unchanged:

| Path | Content |
|---|---|
| `code/plot_next_figs.py` | Figures 1–7 of the second version (imports `plot_main_figs.py` and `plot_schematics.py`) |
| `code/build_preprint_next.py` | typesets the second version (imports `build_preprint_sn.py`) |
| `code/cd4_sc_extract_grid.py`, `depth_grid_cd4.py` | depth grid (Fig. 2d) |
| `code/cd4_trbc_split.py`, `trbc_split_test.py` | constant-region halves (Fig. 3d) |
| `code/cd4_pilot_extract.py`, `code/model_eval/` | transfer models (Fig. 7, Supplementary Table S6) and clonality readouts (Supplementary Table S2) |
| `code/model_eval/generality/`, `jin2020/`, `wholebrain/` | clone strength in 16 screens (Fig. 4c, Supplementary Table S3) |
| `code/model_eval/choose/` | organoid clone barcodes, clone halves and formal tests (Fig. 3f, Fig. 4a, Supplementary Table N1) |
| `code/stats/clone_theory_check.py` | Supplementary Note 1: symbolic checks (sympy) and simulation recovery |
| `pipelines/run_next.sh` | all new stages in order |
| `data/next/` | the result files the new figures and tables read (6.5 MB); `data/next/pilot/` holds the model score files of Fig. 7a, `data/next/inputs/` the run lists of the organoid and whole-brain downloads |
| `paper/preprint_next.md`, `supplementary_next.md`, `references_next_en.md` | text of the second version |
| `paper/figures_next/` | Figures 1–7 as published, for comparison |

### Figures from the included result files

```bash
python code/plot_next_figs.py data data/next figures/next
```

`data/` holds the first-version files that several panels reuse, `data/next/` the new ones. Drawn with Python 3.10,
NumPy 2.2.6 and Matplotlib 3.10.8 (Windows, Arial), all seven figures are identical pixel for pixel to
`paper/figures_next/`; with Matplotlib 3.11 the content is the same but the rendering differs at the pixel level.
`code/build_preprint_next.py paper <fig_dir> <supp_fig_dir> <out.pdf> [--no-author]` typesets the preprint as
`build_preprint_sn.py` does for the first version (Springer Nature template in `paper/sn_template/`, not redistributed;
Supplementary Figures S1–S3 from `make_supplement.py`). Built this way from the files in this repository, the PDF has
the same 48 pages, with the same text on every page, as the authors' build. Supplementary Tables N1 and S1–S7 are
written in `paper/supplementary_next.md`; their values come from the files listed below.

| Manuscript item | Script(s) | Result files (`data/next/` unless noted) |
|---|---|---|
| Fig. 1 (schematics) | `plot_next_figs.py` | — |
| Fig. 2a–c | as Fig. 2a–c of the first version | `data/`: `cd4_D*.npz`, `donor_split.json`, `ntc_null.json`, `cl_ntc_*.json` |
| Fig. 2d: depth grid | `cd4_sc_extract_grid.py` → `depth_grid_cd4.py` | `depth_grid/depth_grid_D1.json`, `depth_grid_D2.json` |
| Fig. 3a–c | as Fig. 3 of the first version | `data/cd4_D1.npz` |
| Fig. 3d: constant-region halves; Table N1, row 1 | `cd4_trbc_split.py` → `trbc_split_test.py` | `trbc_split/*.json` |
| Fig. 3f, Fig. 4a: organoid clone halves; Table N1, rows 2–4 | `choose/export_object.R`, `cbc_parse.py` → `choose_clone.py`, `choose_stats.py` | `choose_clone.json`, `choose_stats.json` |
| Fig. 4b, Fig. 6b, Fig. 7b (TRBC readout), Table S2 | `me_clonality.py`, `generality/unified_null.py` | `clonality.json`, `unified/cd4_*.json` |
| Fig. 4c, Table S3: clone strength in 16 screens | `generality/unified_null.py`, `generality/run_generality.py`, `jin2020/jin_clone.py`, `wholebrain/wb_*.py`, `choose_clone.py` | `unified/*.json`, `generality/*.json`, `jin2020_clone.json`, `wholebrain_wb_results.json`, `wholebrain_wb_split.json`, `choose_clone.json` |
| Fig. 5, Fig. 6 (other panels) | as Fig. 4–5 of the first version | `data/`: `celleval.json`, `clone_D*.json`, `guide_donor.json`, `ci_cd4*.json`, `kd_D*.json`, `replogle_indep.json`, `sim_v6.json` |
| Fig. 7a, Table S6: transfer models, donor 1 | `cd4_pilot_extract.py`, `run_transfer.sh`, `run_state_model.sh`, `pilot_se_embed.py` + `run_state_se.sh`, `pilot_pie_prep.py` + `run_pie.sh` + `pilot_pie_score.py`, `pilot_stack.py`, `pilot_biolord.py`, `pilot_cpa.py`, `pilot_summary.py` | `pilot/*.json`, `model_pilot_summary_table.md` |
| Fig. 7b: four donors | `me_select_targets.py`, `cd4_pilot_extract.py`, `run_transfer.sh` (`me_score_detail.py`) | `detail_D*.json`, `clonality.json` |
| Fig. 7c: answer from another donor | `me_export_preds.py` → `me_cross_donor.py` | `cross_donor.json` |
| Fig. 7d: equal cells, CD4 against melanoma | `me_downsample.py` + `pilot_baselines.py`; `me_frangieh_extract.py` + `run_transfer.sh` | `baselines_ds48_75_D*.json`, `detail_fr_IFNg_f*.json` |
| Expression clonality, fit checks, K562 → RPE1 control (text; `summary_table.md`) | `me_vc_expr.py`, `run_fitcheck.sh` (`me_fitcheck.py`), `me_replogle_extract.py`, `me_summary.py` | `vc_expr_*.json`, `fitcheck_*.json`, `detail_K562_RPE1.json`, `baselines_*.json`, `clone_tests_D*.json`, `pilot_design_*.json`, `summary_table.md` |
| Supplementary Note 1 | `stats/clone_theory_check.py` | `stats/clone_theory_check.json` |

The files in `data/next/` are the original outputs of the runs behind the manuscript (a sample of 10 was checked byte for
byte against the server copies), except that the `file` fields of `unified/*.json` hold only the input file name.

### Full analysis of the second version

`bash pipelines/run_next.sh` runs the new stages after `run_all.sh` (it reuses the per-cell extractions of `extract_cd4`
and the cell-eval subset of `celleval`). New configuration keys are in `config.example.env` (`GPU`, `STATE_BIN`,
`STATE_PYTHON`, `SE_DIR`, `PIE_DIR`, `PIE_BIN`, `RSCRIPT`); `code/model_eval/mkenv_models.sh` creates the model
environments (arc-state 0.11.1, biolord 0.0.3, cpa-tools 0.8.8, arc-pie 1.0.0, arc-stack). The original runs were split
between a CPU server (extraction, clone analyses) and a GPU server (models), with files copied between them; the pipeline
runs the same commands with the same arguments on one machine.

| Stage | What it does |
|---|---|
| `extract_next` | per-cell aggregates of CD4 donors 3 and 4; 16-code aggregates and TRBC classes of donors 1 and 2 |
| `cd4_next` | depth grid, constant-region halves, TRBC clone tests of donors 3 and 4 |
| `unified`, `generality`, `jin2020`, `wholebrain` | clone strength in the 16 screens |
| `organoid` | organoid clone barcodes (streamed from ENA), clone halves, formal tests |
| `stats` | Supplementary Note 1 checks |
| `model_inputs` | cell subsets of CD4 donors 1–4, K562 → RPE1 and the Frangieh folds |
| `model_gpu` | baselines, STATE, scoring, prediction export and fit check per run |
| `model_pilot` | donor 1: Arc perturb_mean and pseudobulk, State with SE-600M, PIE, Stack, biolord, CPA |
| `model_post` | clonality, cross-donor control, expression clonality, matched-cell baselines |
| `collect_next`, `figures_next` | gather the result files under the names in `data/next/`; draw Figures 1–7 |

Additional data (all public; Supplementary Table S7), under `DATA_ROOT`:

| Dataset | Source | Files under `DATA_ROOT` |
|---|---|---|
| Zhu et al., donors 3 and 4 | as donors 1 and 2 | `cd4_sc/{D3,D4}_{Rest,Stim48hr}.assigned_guide.h5ad` |
| Frangieh et al. 2021 | scPerturb, Zenodo 13350497 | `frangieh/FrangiehIzar2021_RNA.h5ad` |
| Li et al. 2023 (CHOOSE organoids) | Zenodo 7083558 (`CHOOSE_ASD_Modules.rds`); clone-barcode reads of ArrayExpress E-MTAB-13148, streamed by `cbc_parse.py` from the ENA runs in `data/next/inputs/choose_cbc_runs.txt` | `choose/CHOOSE_ASD_Modules.rds` |
| Jin et al. 2020 | GEO GSE157977 (RAW tar, unpacked); barcode table from the paper's Table S5 (`data/next/jin2020_barcodes.csv`) | `jin2020/` |
| Shi et al. (whole brain) | Hugging Face `perturbai/wholebrain_crispr_atlas`: 14 h5ad files (sizes and SHA-256 in `data/next/inputs/wholebrain_h5ads.txt`) and `all_obs.parquet` | `wholebrain/h5ads/`, `wholebrain/all_obs.parquet` |
| Zhou et al. 2023 | GEO GSE216800 (GSM6697660–GSM6697675, CellRanger tar files) | `zhou/` |
| Knudsen et al. | GEO GSE302504 (`GSE302504_adata_106_21_GEX_final.h5ad`) | `knudsen/` |
| Fagerberg et al. 2025; Schmidt et al. 2022 | FITdb, Zenodo 22212134 (`Fagerberg_2025.h5ad`, `Schmidt_2022_resting.h5ad`, `Schmidt_2022_restim.h5ad`) | `fitdb/` |
| Cao et al. | GEO GSE221578 (CellRanger tar files of the in vitro and in vivo samples) | `cao/` |

What has not been rerun, and what the pipeline cannot recreate exactly:

- The pipeline has not been run end to end from raw data. Checked: the figures (identical from `data/next/`), the
  preprint PDF, and three clone-strength runs recomputed with the released `unified_null.py` and the arguments in
  `run_next.sh` (HepG2, Frangieh control, H1), which equal the result files exactly. The other commands in `run_next.sh`
  were taken from the run scripts and logs on the two servers.
- Model training on GPUs is not bit-reproducible, so a rerun of `model_gpu` and `model_pilot` gives scores that differ
  in the last digits.
- In the original run of the donor-1 depth grid, `cd4_sc_extract_grid.py` was pointed at the per-cell extraction made
  before the guide-parity fix (`REPRODUCTION.md`); the script takes from it only the gene list, control means and cell
  counts, which the fix did not change, and computes guide parity itself. The pipeline uses the current extraction.
- Fig. 7a and Supplementary Table S6 (donor-1 model pilot): the PIE knowledge sources (the published sources without
  `ncbi_text`), PIE's preprocessing and its `pie eval` step, and the SE-600M weights were prepared by hand following the
  packages' documentation; the TOML of the SE run is written by `run_next.sh` as the donor-1 TOML pointed at the SE
  inputs. The epoch setting of the adapted CPA run is not recorded. CPA and biolord (tutorial settings) did not learn
  per-sgRNA effects (fit check r ≤ 0.025) and appear in the table only.
- `cd4_pilot_extract.py` and `unified_null.py` are the latest versions. The versions that ran for donor 1 and for the CD4
  and Replogle screens lacked, respectively, the `PILOT_DONOR` option (default donor 1) and the reader for
  column-compressed matrices, which those runs do not use.

## License

The code is released under the MIT License (see `LICENSE`). If you use the code or the result files in `data/`, please cite the preprint.
