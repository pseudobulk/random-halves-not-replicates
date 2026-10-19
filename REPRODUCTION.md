# Reproduction from raw data (2026-10-19)

The whole analysis was rerun from the raw public data with the code in this repository, following
`pipelines/run_all.sh`. Every result was compared with the original outputs behind the preprint, and every printed
number, figure and the Supplementary Information was checked against preprint v16. Preprint v17 differs from v16 only in the layout of the supplementary tables (ranges
written with an en dash, as in the main text, and no table split across pages); every value is unchanged. `repro/comparison.tsv` lists each
comparison; `repro/number_inventory.tsv` classifies every number in the manuscript text (422), and
`repro/printed_numbers.md` checks the 182 printed results one by one (`repro/intermediate/` holds the few values
taken from intermediate outputs); `repro/tools/` holds
the comparison scripts.

## Acceptance criteria

- Deterministic stages must match bit for bit.
- Stages whose float32 sums depend on accumulation order pass only if every number printed in the manuscript, the
  Supplementary Information and the authors' number ledger is identical at its printed precision, and every figure is
  visually identical, with any pixel differences at noise level and the largest difference stated.
- Any printed number that changes is a mismatch and is reported with both values and its cause.

## How it was run

- Code: this repository; configuration in `config.env` as documented. Raw data were read in place, read-only, through a
  folder of symbolic links laid out as `DATA_ROOT` in the README; all outputs went to a new folder.
- Machines: a 256-core server (Python 3.11.11, NumPy 2.4.6, SciPy 1.17.1, h5py 3.16, PyArrow 20) for all stages except
  two: the analyses of the Nadig and H1 files ran on a second server where those files are stored (Python 3.11.16, NumPy 2.4.6,
  SciPy 1.16.3, h5py 3.16), and cell-eval 0.8.2 ran there in its own environment (Python 3.11.16, AnnData 0.12.19).
  Figures were drawn on Windows (Python 3.14, NumPy 2.5, Matplotlib 3.11, Arial).
- Wall time with `WORKERS=32` (16 for Tahoe): `extract_cd4` 69 min, `extract_tahoe` 38 min, `cd4` 125 min,
  `cd4_pseudobulk` 33 min, `replogle` 36 min, `tahoe` 26 min, `sim` 2 min, cell-eval 90 min per run.

## Results

| Stage | Outputs compared | Result |
|---|---|---|
| CD4 extraction (4 files × random halves, guide halves, per-cell, non-targeting guides) | all arrays | cell-level content (target, guide parity, random bits, TRBC counts, cell counts) bit-identical; expression sums differ by float32 accumulation order only (relative ≤ 7 × 10⁻⁷) |
| Tahoe extraction (3,388 shards, 95,624,334 cells) | aggregates | counts bit-identical; sums relative ≤ 1.5 × 10⁻⁷ (accumulation order; the original run was batched as shards arrived) |
| CD4 analyses, donors 1 and 2; non-targeting null; TRBC across states and effective-clone bound | 26 result files | float noise only (largest relative difference 1.6 × 10⁻⁴, in near-zero medians); no value rounds differently at 1–3 significant digits; per-target arrays without a single sign change |
| Four-donor pseudobulk (donor halves, same vs different guide) | 3 | bit-identical |
| Tahoe analyses (replicate plates, DMSO wells, ceilings, robustness, intervals, independent check) | 7 | bit-identical, except the DMSO well test (relative 4 × 10⁻⁸) |
| Replogle (independent screens, K562 across screens, template check, K562–RPE1) | 5 | bit-identical |
| Non-targeting constructs in five cell-line screens (Replogle, Nadig) and the VCC 2025 H1 hESC data | 3 | float noise only (relative ≤ 3 × 10⁻⁶; the worker partials are merged in completion order); no value rounds differently; run times and file names differ by design |
| Nadig (Jurkat–HepG2 pass test) | 1 | bit-identical |
| Simulation | 1 | identical up to 4 × 10⁻¹⁶ |
| cell-eval 0.8.2 | 1 | see below |
| Printed numbers (manuscript, Supplementary Information, ledger) | 182 | all identical between the original and the rerun at printed precision; every value found in the v16 text, supplement or ledger |
| Figures 1–5, Supplementary Figures S1–S3 | 10 images | identical from `data/`; from the rerun, identical except Fig. 3 (6 pixels by 1/255, anti-aliasing) |
| Supplementary Information text (Tables S1–S6) | 201 lines | identical |

Files present only in the original output folders were written by scripts outside this repository
(`cd4_sc_analyze.py` and two exploratory Tahoe scripts); no number in the paper depends on them.

## cell-eval: the one input built before an audit fix

The reproduction found that the original cell-eval run (Fig. 4a, Supplementary Table S3 up to preprint v11) used a cell
subset built from a per-cell extraction made before the 10-18 code audit, in which guide parity followed the guide's
rank instead of its name. Rebuilding that subset from the pre-audit extraction gives a byte-identical file, and running
cell-eval on it reproduces the original output exactly (1,476 values). The released pipeline builds the subset from
the corrected extraction, as the Methods describe: 60,238 instead of 60,310 cells, median 62 instead of 66 cells per
side and a different random draw of the 150 unexpressed targets (the 150 strong knockdowns are the same). Eight
printed cell-eval values changed in the second decimal and one guide-half ceiling became undefined (split-half
correlation −0.01; cell-eval reports no ceiling for r ≤ 0). Preprints v12 to v16 use the released result, `data/celleval.json`,
and all its conclusions hold. All other inputs were produced from post-audit extractions or by scripts the audit did
not change, and they reproduce as listed above.

## Changes to the pipeline made during the reproduction

- `run_all.sh` and `download.sh`: variables already exported now take precedence over `config.env`, as documented.
- `run_all.sh`: creates the output folder of `cd4_ntc_null.py`.
- Neither change affects a result.

## Not covered

- Supplementary Table S6 is text written before the donor-2 analysis; it is not computed.
- `code/build_preprint_sn.py` was not rerun here; it needs the Springer Nature LaTeX template, which is not part of this
  repository.
