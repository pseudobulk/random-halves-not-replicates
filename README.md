# random-halves-not-replicates

Analysis code for the preprint

**Random cell halves are not independent replicates in pooled perturbation screens**
Xuanrui Huang, National University of Singapore
Preprint: bioRxiv (link will be added on posting)

## Status

The code is being prepared for release and will be published in this repository. Until then, the Methods and Supplementary Information of the preprint describe every analysis in enough detail to reimplement it.

## What will be released

- Extraction of per-cell and pseudobulk aggregates from the public datasets: CD4+ T-cell CRISPRi Perturb-seq (Zhu et al. 2026), Tahoe-100M, Replogle et al. 2022, Nadig et al. 2025.
- Null-perturbation, clone-marker (TRBC1/TRBC2) and depth analyses.
- Same-guide against different-guide comparison across donors.
- Ceilings: cell-eval 0.8.2 with random and guide halves, and the reliability-label comparison.
- Shared-to-specific predictor and robustness checks.
- The three-term simulation model.
- Scripts that regenerate all main and supplementary figures and tables.
- Environment specifications and a data manifest.

## Data

All datasets are public; the preprint's Data availability section gives the sources.
