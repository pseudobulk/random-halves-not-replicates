# Supplementary Information

Contents: Supplementary Note 1; Supplementary Tables S1–S7; Supplementary Figures S1–S3.

## Supplementary Note 1 | Statistical model of clone strength: propositions and checks

**Model.** Consider one perturbation measured on G genes. Its cells come from K clones; clone k has n~k~ cells, a share p~k~ = n~k~/N of all N cells. Cell i belongs to clone k(i), and its expression vector, after subtracting a fixed reference, is x~i~ = μ + c~k(i)~ + e~i~. Here μ is the perturbation's true effect; c~k~ is the deviation of clone k, independent between clones, with mean zero and variance σ~c~^2^ per gene; e~i~ is the cell's own variation and sampling noise, independent between cells, with mean zero and variance σ~e~^2^ per gene. The intraclass correlation is ICC = σ~c~^2^/(σ~c~^2^ + σ~e~^2^), the share of the differences between single cells that is set by clones. A null perturbation has μ = 0. Clone deviations and cell noise are assumed independent of each other and of the perturbation, and the reference is estimated from many cells, so that its error can be neglected (the handling of reference error is described in Methods, Clone strength).

**Proposition 1 (expected cross-half signal).** Let A and B be two sets of cells of a perturbation with no cell in common, in which clone k has shares a~k~ and b~k~. Then E[$\bar{x}_A$ · $\bar{x}_B$] = |μ|^2^ + Gσ~c~^2^ Σ~k~ a~k~b~k~. If A and B share no clone (independent replicates), the second term is zero and only the energy of the true effect, |μ|^2^, remains. *Reason.* $\bar{x}_A$ = μ + Σ~k~ a~k~c~k~ + $\bar{e}_A$. In the expanded inner product only the terms of the same clone, E[c~k~ · c~k~] = Gσ~c~^2^, are non-zero; A and B have no cell in common, so $\bar{e}_A$ and $\bar{e}_B$ are uncorrelated and do not enter. *Check.* With two clones and symbolic shares, the expansion by sympy equals the right-hand side.

**Proposition 2 (random halves and clone strength).** When the N cells are split into two halves by a fair coin, E[Σ~k~ a~k~b~k~] = Σ~k~ p~k~^2^ (1 + O(1/N)). For a null perturbation, the difference between the cross-half signals of random halves and independent replicates, divided by the per-cell variance G(σ~c~^2^ + σ~e~^2^), is therefore the clone strength, ICC × Σ~k~ p~k~^2^ = ICC/K~eff~ with K~eff~ = 1/Σ~k~ p~k~^2^, and it does not change with the number of cells. *Reason.* The number of cells of clone k in the first half is binomial, Bin(n~k~, 1/2), and each half has about N/2 cells, so the leading term is (n~k~/2)(n~k~/2)/((N/2)(N/2)) = p~k~^2^. *Check.* The leading term was verified with sympy. The estimator was tested in simulated screens with known truth: ICC of 0.01, 0.05 and 0.2, 5–200 clones, coefficients of variation of clone size of 0.5–2, 100 or 400 cells per perturbation and 60 null perturbations. Computed as on the data (random halves minus two different null perturbations, divided by the per-cell variance), clone strength correlated with the true ICC × Σ p~k~^2^ at 0.999, with a ratio of 0.98.

**Proposition 3 (effective number of clones, transduction and expansion).** If N~0~ cells were transduced per perturbation and the clone sizes after expansion have mean m and coefficient of variation CV, then to leading order Σ~k~ p~k~^2^ = (1 + CV^2^)/N~0~, that is, K~eff~ = N~0~/(1 + CV^2^). *Reason.* Σ~k~ p~k~^2^ = Σ~k~ n~k~^2^/(Σ~k~ n~k~)^2^ ≈ N~0~E[n^2^]/(N~0~m)^2^, and E[n^2^] = m^2^(1 + CV^2^). *Check.* Verified with sympy. The relation is approximate: it treats the sequenced cells as a proportional reduction of the clone sizes, and when only a small fraction of cells is sequenced, fewer than N~0~ clones are observed.

**Proposition 4 (limit of the Spearman–Brown extrapolation).** For a null perturbation with a fixed number of clones and σ~c~ > 0, the correlation over genes between the two half-means tends to 1 as the number of cells grows, and so does the Spearman–Brown corrected reliability; yet its effect in an independent experiment is zero. The Spearman–Brown formula extrapolates the agreement within the same clones. *Reason.* With more cells the noise terms of both halves vanish, and both halves tend to the same Σ~k~ p~k~c~k~. *Check.* In a simulated null perturbation from 4 clones with an ICC of 0.05, the split-half correlation rose from 0.23 to 0.95 and the corrected reliability from 0.38 to 0.98 as the number of cells rose from 50 to 3,200.

**Proposition 5 (leakage across contexts).** Let two contexts R and S contain the same clones, with clone k deviating by c~k~ in R and by qc~k~ in S. Let the common part be s = (x~R~ + x~S~)/2 and the specific part d = (x~R~ − x~S~)/2, taken from two sets of cells A and B with no cell in common. Then E[s~A~ · d~B~] = Gσ~c~^2^ (1 + q)(1 − q)/4 × Σ~k~ a~k~b~k~. When q ≠ 1 and the two sets share clones this is non-zero, so a predictor of d~B~ from s~A~ can score positively even without any true cross-context relation; when they share no clone it is zero. *Check.* Verified with sympy.

The formal tests reported in the main text are summarized in Supplementary Table N1.

<!-- TABLE N1 -->
**Supplementary Table N1 | Formal tests of the clone mechanism and of the clone-strength relation.**

| Test | Statistic | Null distribution | Result |
|---|---|---|---|
| CD4, constant-region halves | random minus constant-region halves, per target | paired sign test; Wilcoxon signed-rank test | p < 10^−150^ in each donor and state |
| Organoids, clone halves | random minus clone halves, per perturbation | sign test; clone labels shuffled within each perturbation (200 times) | 32 of 32, p = 2 × 10^−10^ (1,000 splits per perturbation); shuffled clone halves 4.98, random halves 4.97, true clone halves 3.43; permutation p = 0.005 |
| Relation: slope equals ICC | slope of clone strength on Σ p~k~^2^ | bootstrap over perturbations (2,000) | slope 0.0096 (0.0085–0.0109); median ICC 0.0104 (0.0034–0.0171); interval of their ratio 0.57–2.70 includes 1 |
| Relation: out-of-sample prediction | ICC of one half of the perturbations predicts the clone strength of the other | random splits of perturbations (500) | median Spearman 0.80, Pearson 0.83; median measured/predicted 0.94 (quartiles 0.71–1.39) |

## Supplementary Table S1 | Non-targeting sgRNAs as null perturbations

<!-- TABLE S1a -->
**Supplementary Table S1 | Non-targeting sgRNAs as null perturbations. a, CD4.** Each non-targeting sgRNA of the CD4 screen is a separate transduction with no target. Effects are referenced to the mean over the kept non-targeting sgRNAs in each half; genes with a mean of at least 0.1 over non-targeting cells, TCR genes excluded; sgRNAs with at least 10 cells per half. Different sgRNAs: the first half of one sgRNA against the second half of another, at the same depth, averaged over 50 random pairings (range over pairings in parentheses). Random halves: bootstrap 95% interval over sgRNAs. TRBC: variance of the z-score of each sgRNA's TRBC1 share (sgRNAs with at least 30 cells carrying TRBC counts; in parentheses after permuting sgRNA labels) and the Pearson correlation of the TRBC1 deviation between random halves (at least 15 cells per half).

| Donor, state | sgRNAs | Cells per half (median) | Positive, random halves | Positive, different sgRNAs | Mean signal, random / different | TRBC var(z) (permuted) | TRBC deviation, between random halves |
|---|---|---|---|---|---|---|---|
| Donor 1, resting | 739 | 39 | 98.2% (97.2–99.1) | 49.8% (42.8–54.7) | 13.9 / −0.04 | 34.4 (0.99) | 0.95 |
| Donor 1, stimulated | 785 | 36 | 96.6% (95.3–97.8) | 49.8% (45.2–55.1) | 13.5 / −0.01 | 29.1 (1.01) | 0.94 |
| Donor 2, resting | 578 | 41 | 98.8% (97.9–99.5) | 48.8% (43.6–56.7) | 20.2 / −0.13 | 79.9 (0.98) | 0.97 |
| Donor 2, stimulated | 622 | 40 | 98.7% (97.7–99.5) | 49.8% (43.4–56.6) | 21.3 / 0.01 | 76.0 (1.00) | 0.97 |

<!-- TABLE S1b -->
**Supplementary Table S1 (continued) | b, Cell-line screens.** Raw single cells of the K562 essential-gene, K562 genome-wide (part 1) and RPE1 screens of Replogle et al., the Jurkat and HepG2 screens of Nadig et al. and the H1 hESC training data of the 2025 Virtual Cell Challenge; one pseudo-perturbation per non-targeting construct. The last two columns compare the mean random-half signal of non-targeting constructs with that of targeting perturbations of the same screen (CD4: resting cells, TCR genes excluded).

| Screen | Constructs | Cells per half (median) | Positive, random halves | Positive, different constructs | Mean signal, non-targeting / targeting | Ratio |
|---|---|---|---|---|---|---|
| CD4 donor 1, resting | 739 | – | see a | see a | 13.9 / 10.5 | 1.31 |
| CD4 donor 2, resting | 578 | – | see a | see a | 20.2 / 17.4 | 1.16 |
| K562 essential (day 6) | 97 | 52 | 74.2% (64.9–82.5) | 49.1% (33.3–72.9) | 0.36 / 15.7 | 0.023 |
| K562 genome-wide (day 8, part 1) | 507 | 68 | 71.4% (67.6–75.0) | 49.9% (42.7–56.9) | 0.28 / 4.4 | 0.064 |
| RPE1 essential | 113 | 49 | 76.1% (68.1–84.1) | 50.9% (37.5–69.6) | 1.06 / 42.1 | 0.025 |
| Jurkat essential | 128 | 44 | 73.4% (65.6–81.2) | 51.4% (34.4–64.1) | 1.39 / 20.0 | 0.070 |
| HepG2 essential | 112 | 19 | 59.8% (50.0–68.8) | 50.0% (33.9–71.4) | 1.25 / 29.6 | 0.042 |
| H1 hESC, VCC 2025 training data | 31 | 574 | 71.0% (54.8–87.1) | 49.2% (13.3–80.0) | 0.09 / 16.7 | 0.006 |

## Supplementary Table S2 | Two clone readouts in the CD4 screen

<!-- TABLE S2 -->
**Supplementary Table S2 | Two clone readouts in the CD4 screen.** TRBC readout: the variance of each sgRNA's mean TRBC1 share beyond its sampling variance, divided by the per-cell variance of the share and averaged over sgRNAs. Expression readout: the clone strength of non-targeting sgRNAs (random-half minus different-sgRNA cross-half signal, divided by the per-cell variance), with effects referenced to the null cells or to all assigned cells of the same sample, and their geometric mean. Both readouts are expected to equal their own ICC divided by the same effective number of clones, so their ratio should not depend on the donor. Donors are ordered by the TRBC readout.

| Donor, state | TRBC readout | Clone strength, null reference | Clone strength, all-cell reference | Clone strength, geometric mean | Ratio | Null sgRNAs |
|---|---|---|---|---|---|---|
| Donor 4, resting | 0.156 | 0.0044 | 0.0044 | 0.0044 | 0.028 | 821 |
| Donor 4, stimulated | 0.166 | 0.0049 | 0.0048 | 0.0049 | 0.029 | 875 |
| Donor 1, resting | 0.277 | 0.0071 | 0.0071 | 0.0071 | 0.026 | 688 |
| Donor 1, stimulated | 0.275 | 0.0080 | 0.0080 | 0.0080 | 0.029 | 723 |
| Donor 2, resting | 0.429 | 0.0103 | 0.0103 | 0.0103 | 0.024 | 541 |
| Donor 2, stimulated | 0.425 | 0.0105 | 0.0105 | 0.0105 | 0.025 | 587 |
| Donor 3, resting | 0.535 | 0.0120 | 0.0120 | 0.0120 | 0.022 | 507 |
| Donor 3, stimulated | 0.519 | 0.0107 | 0.0107 | 0.0107 | 0.021 | 515 |

The TRBC readout differs 3.4-fold and the expression readout 2.7-fold between donors, while the ratio varies between 0.021 and 0.029 and falls slightly as the readouts rise.

## Supplementary Table S3 | Clone strength in 16 screens

<!-- TABLE S3 -->
**Supplementary Table S3 | Clone strength in 16 screens.** One estimator for all screens with several null sgRNAs: the random-half cross-half signal of null sgRNAs (mean of 20 splits) minus the cross-half signal between two different null sgRNAs (50 pairings), divided by the per-cell variance of the null cells; at least 30 cells per sgRNA. Ranges are over the two references (null cells or all assigned cells of the same sample) and over states, treatments or donors. Null / targeting: the random-half signal of null sgRNAs relative to that of a typical targeting perturbation of the same screen. Screens with a single null sgRNA were estimated with halves that share no clone (Methods). Screens are ordered as in Fig. 4c.

| Screen | System | Null sgRNAs | Clone strength | Null / targeting |
|---|---|---|---|---|
| Knudsen et al. | human CAR-T cells, 21 days in mice | 4 | 0.045–0.10 | 1.0–3.3 |
| Fagerberg et al. | mouse P14 CD8 T cells, 8 days in LCMV-infected mice | 5 | 0.010–0.015 | 2.0–2.4 |
| Zhu et al. | primary human CD4 T cells, genome-wide CRISPRi, 4 donors, resting and stimulated | 507–875 | 0.0044–0.012 | 0.83–1.05 |
| Schmidt et al. | primary human T cells, CRISPRa, resting and restimulated | 8 | 0.0013–0.0029 | 0.10–0.29 |
| Zhou et al. | mouse OT-I CD8 T cells, 7 days in tumour-bearing mice | 37 | 0.0013–0.0017 | 0.44–0.48 |
| Cao et al. | primary human AML cells, in vitro | 10 | 0.0008–0.0018 | 0.07–0.12 |
| Jin et al. | embryonic mouse cortex, in utero lentivirus | 1 | 0.001–0.003 (library halves, upper bound) | – |
| Li et al. | human brain organoids, clone barcodes | 1 | 0.0009 (clone halves; 0.0046 with cell-type composition) | – |
| Nadig et al. | HepG2, essential genes, CRISPRi | 87 | 0.0011–0.0015 | 0.04–0.20 |
| Nadig et al. | Jurkat, essential genes, CRISPRi | 127 | 0.0010–0.0012 | 0.06–0.16 |
| Replogle et al. | RPE1, essential genes, CRISPRi | 112 | 0.0011–0.0012 | 0.02–0.18 |
| Frangieh et al. | melanoma cell line, CRISPR knockout, three treatments | 71–74 | 0.0003–0.0009 | 0.35–0.82 |
| Replogle et al. | K562, essential genes, CRISPRi | 95 | 0.0004–0.0006 | 0.03–0.16 |
| Replogle et al. | K562, genome-wide, CRISPRi | 503 | 0.0003 | 0.09 |
| VCC 2025 | H1 hESC, CRISPRi, public training data | 31 | 0.0002–0.0003 | 0.01–0.02 |
| Shi et al. | mouse whole brain, AAV at postnatal day 16, non-dividing neurons; reference within mouse | 399 | 0.0001–0.0002 | 0.35–0.56 |

## Supplementary Table S4 | cell-eval 0.8.2 ceilings for all metrics

<!-- TABLE S4a -->
**Supplementary Table S4 | cell-eval 0.8.2 ceilings for all metrics. a, Spearman–Brown corrected ("ceiling").** Resting CD4+ T cells of donor 1: all 300 perturbations, the 150 unexpressed targets and the 150 strength-matched knockdowns, with equal numbers of cells per sgRNA (median 62 per side). Random: the tool's compute_ceiling (disjoint random halves of each perturbation and of the controls), range over seeds 0–2. sgRNA: the tool's compute with one sgRNA's cells as "real" and the other's as "prediction", range over the two directions. Asterisks, metrics whose chance level is 0.5, which the correction maps to 0.67. Daggers, metrics that cell-eval computes over the whole data set, so they are identical across subsets. Dashes, metrics not returned for a subset. n.d., not defined: the split-half correlation was at most 0; the half-depth correlation is given in parentheses.

| Metric | All: random | All: sgRNA | Unexpressed: random | Unexpressed: sgRNA | Knockdowns: random | Knockdowns: sgRNA |
|---|---|---|---|---|---|---|
| DE overlap (top N) | 0.30–0.31 | 0.11 | 0.20–0.22 | 0.05 | 0.37–0.40 | 0.16–0.17 |
| DE overlap (top 50) | 0.30–0.32 | 0.11 | 0.20–0.22 | 0.05 | 0.38–0.41 | 0.17 |
| DE overlap (top 100) | 0.30–0.31 | 0.11 | 0.20–0.22 | 0.05 | 0.38–0.40 | 0.17 |
| DE overlap (top 200) | 0.30–0.31 | 0.11 | 0.20–0.22 | 0.05 | 0.37–0.40 | 0.16–0.17 |
| DE overlap (top 500) | 0.30–0.31 | 0.11 | 0.20–0.22 | 0.05 | 0.37–0.40 | 0.16–0.17 |
| DE precision (top N) | 0.30–0.31 | 0.11 | 0.19–0.22 | 0.05 | 0.39 | 0.16–0.17 |
| DE precision (top 50) | 0.30–0.31 | 0.11 | 0.20–0.22 | 0.05 | 0.39–0.40 | 0.17 |
| DE precision (top 100) | 0.30–0.31 | 0.11 | 0.19–0.22 | 0.05 | 0.39 | 0.17 |
| DE precision (top 200) | 0.30–0.31 | 0.11 | 0.19–0.22 | 0.05 | 0.39 | 0.16–0.17 |
| DE precision (top 500) | 0.30–0.31 | 0.11 | 0.19–0.22 | 0.05 | 0.39 | 0.16–0.17 |
| DE recall | 0.36–0.40 | 0.15–0.18 | 0.28–0.30 | 0.09–0.11 | 0.43–0.49 | 0.21–0.26 |
| DE direction match* | 0.84–0.85 | 0.62–0.64 | 0.81–0.82 | 0.57–0.61 | 0.86–0.88 | 0.67 |
| DE Spearman (significant genes)† | 0.80–0.82 | 0.46–0.47 | 0.80–0.82 | 0.46–0.47 | 0.80–0.82 | 0.46–0.47 |
| DE Spearman of log fold changes† | – | – | – | – | 0.42 | n.d. (r = −0.05) |
| PR-AUC | – | – | – | – | 0.39–0.43 | 0.13–0.16 |
| ROC-AUC* | – | – | – | – | 0.84–0.86 | 0.76–0.77 |
| Pearson Δ | 0.32 | 0.06 | 0.26 | n.d. (r = −0.01) | 0.38–0.39 | 0.13 |
| Discrimination (L1)* | 0.83 | 0.69 | 0.81–0.82 | 0.67–0.69 | 0.84 | 0.69–0.70 |
| Discrimination (L2)* | 0.86–0.87 | 0.70 | 0.85–0.86 | 0.68–0.70 | 0.88 | 0.70–0.72 |
| Discrimination (cosine)* | 0.97–0.98 | 0.69–0.70 | 0.97 | 0.63–0.65 | 0.98 | 0.74 |

<!-- TABLE S4b -->
**Supplementary Table S4 (continued) | b, Uncorrected half-depth means.** As in a, before the Spearman–Brown correction.

| Metric | All: random | All: sgRNA | Unexpressed: random | Unexpressed: sgRNA | Knockdowns: random | Knockdowns: sgRNA |
|---|---|---|---|---|---|---|
| DE overlap (top N) | 0.18 | 0.06 | 0.11–0.12 | 0.03 | 0.23–0.25 | 0.09 |
| DE overlap (top 50) | 0.18–0.19 | 0.06 | 0.11–0.13 | 0.03 | 0.23–0.25 | 0.09 |
| DE overlap (top 100) | 0.18 | 0.06 | 0.11–0.12 | 0.03 | 0.23–0.25 | 0.09 |
| DE overlap (top 200) | 0.18 | 0.06 | 0.11–0.12 | 0.03 | 0.23–0.25 | 0.09 |
| DE overlap (top 500) | 0.18 | 0.06 | 0.11–0.12 | 0.03 | 0.23–0.25 | 0.09 |
| DE precision (top N) | 0.18 | 0.06 | 0.11–0.12 | 0.03 | 0.24–0.25 | 0.09 |
| DE precision (top 50) | 0.18–0.19 | 0.06 | 0.11–0.13 | 0.03 | 0.25 | 0.09 |
| DE precision (top 100) | 0.18 | 0.06 | 0.11–0.12 | 0.03 | 0.24–0.25 | 0.09 |
| DE precision (top 200) | 0.17–0.18 | 0.06 | 0.11–0.12 | 0.03 | 0.24 | 0.09 |
| DE precision (top 500) | 0.18 | 0.06 | 0.11–0.12 | 0.03 | 0.24 | 0.09 |
| DE recall | 0.22–0.25 | 0.08–0.10 | 0.16–0.17 | 0.05–0.06 | 0.27–0.33 | 0.12–0.15 |
| DE direction match* | 0.72–0.74 | 0.45–0.47 | 0.68–0.69 | 0.40–0.44 | 0.75–0.79 | 0.50 |
| DE Spearman (significant genes)† | 0.67–0.69 | 0.30 | 0.67–0.69 | 0.30 | 0.67–0.69 | 0.30 |
| DE Spearman of log fold changes† | – | – | – | – | 0.26 | −0.05 |
| PR-AUC | – | – | – | – | 0.24–0.27 | 0.07–0.09 |
| ROC-AUC* | – | – | – | – | 0.72–0.76 | 0.61–0.62 |
| Pearson Δ | 0.19 | 0.03 | 0.15 | −0.01 | 0.23–0.24 | 0.07 |
| Discrimination (L1)* | 0.70–0.71 | 0.52–0.53 | 0.68–0.70 | 0.51–0.52 | 0.72–0.73 | 0.53–0.54 |
| Discrimination (L2)* | 0.76–0.77 | 0.54 | 0.74–0.75 | 0.51–0.54 | 0.78 | 0.54–0.56 |
| Discrimination (cosine)* | 0.95 | 0.53 | 0.94 | 0.46–0.48 | 0.96–0.97 | 0.59 |

## Supplementary Table S5 | Structure of unexpressed targets carried by the common-to-specific predictor

<!-- TABLE S5 -->
**Supplementary Table S5 | Structure of unexpressed targets carried by the common-to-specific predictor.** For a pair of contexts X and Y, the predictor maps the common part of half A, (X~A~ + Y~A~)/2, to the specific part of half B, (X~B~ − Y~B~)/2, and gives a held-out prediction F for every perturbation. Because the specific part contains +X/2, a predictor that reads a perturbation's half-A deviation in context X passes that deviation on to F. For a set of targets S, the pass share is 2 Σ~S~ F · X~B~/Σ~S~ X~A~ · X~B~ (with the opposite sign for context Y): the share of the reproducible energy of the targets' own responses in X that reappears in the prediction. Responses were referenced to the mean of the unexpressed targets, and the direction of the mean of all perturbations (the template) was then projected out. The predictor used the top k principal components of the half-A common part, ridge regression with λ chosen by inner two-fold cross-validation, and two outer folds over perturbations (seed 11). CD4 halves are random halves or sgRNA halves; in the cell-line pairs only random halves exist. Intervals are 2.5th–97.5th percentiles over 500 bootstrap resamples of targets. The shuffled control gave values of at most 0.0079 in absolute value.

| Context | Halves | Unexpressed, n | Pass, k = 20 (95% CI) | k = 10 | Expressed, n | Pass, k = 20 (95% CI) | k = 10 |
|---|---|---|---|---|---|---|---|
| CD4 donor 1, resting | random | 1,671 | 0.27 (0.26–0.28) | 0.21 | 2,120 | 0.25 (0.24–0.26) | 0.19 |
| CD4 donor 1, stimulated | random | 1,671 | 0.23 (0.22–0.24) | 0.18 | 2,239 | 0.21 (0.20–0.22) | 0.16 |
| CD4 donor 1, resting | sgRNA | 1,241 | 0.00 (−0.01 to 0.01) | 0.00 | 1,424 | 0.00 (0.00–0.00) | 0.00 |
| CD4 donor 1, stimulated | sgRNA | 1,241 | 0.00 (0.00–0.00) | 0.00 | 1,511 | 0.00 (0.00–0.00) | 0.00 |
| CD4 donor 2, resting | random | 1,477 | 0.29 (0.28–0.30) | 0.23 | 1,974 | 0.27 (0.27–0.28) | 0.22 |
| CD4 donor 2, stimulated | random | 1,477 | 0.26 (0.25–0.28) | 0.22 | 1,951 | 0.26 (0.25–0.28) | 0.22 |
| CD4 donor 2, resting | sgRNA | 874 | 0.00 (−0.01 to 0.01) | 0.00 | 1,124 | 0.00 (0.00–0.00) | 0.00 |
| CD4 donor 2, stimulated | sgRNA | 874 | 0.00 (−0.03 to 0.02) | −0.01 | 1,068 | 0.00 (0.00–0.00) | 0.00 |
| K562 (K562–RPE1) | random | 181 | 0.18 (0.06–0.30) | 0.15 | 732 | 0.17 (0.13–0.20) | 0.14 |
| RPE1 (K562–RPE1) | random | 222 | 0.29 (0.25–0.33) | 0.24 | 608 | 0.38 (0.34–0.41) | 0.32 |
| Jurkat (Jurkat–HepG2) | random | 209 | 0.14 (0.10–0.17) | 0.10 | 773 | 0.23 (0.21–0.25) | 0.19 |
| HepG2 (Jurkat–HepG2) | random | 199 | 0.23 (0.15–0.31) | 0.21 | 656 | 0.32 (0.30–0.35) | 0.26 |

Under random halves in CD4, the predictor carried the structure of unexpressed targets at the same rate as that of expressed targets; under sgRNA halves it carried neither. In the four cell lines, the pass share of unexpressed targets was 0.6–1.1 times that of expressed targets. In these libraries, which were selected on a growth phenotype, the structure of unexpressed targets largely comes from the sgRNAs or target sites (Supplementary Fig. S3), so carrying it does not by itself indicate clones.

## Supplementary Table S6 | Transfer models: scores, settings and fit checks

<!-- TABLE S6 -->
**Supplementary Table S6 | Transfer models in CD4 donor 1: scores against the two answers.** A, answer from the same sgRNA (shares clones with the input); B, answer from the other sgRNA against the same target (another transduction). Scores are Pearson correlations over the evaluation genes between predicted and answer effects, averaged over test sgRNAs; standard errors of A − B are clustered by target. For reference, the stimulated effects of the two sgRNAs of a target agree at r = 0.041 (unexpressed targets) and 0.059 (knockdowns).

| Model | Unexpressed: A | Unexpressed: B | Unexpressed: A − B (SE) | Knockdowns: A | Knockdowns: B | Knockdowns: A − B (SE) | Note |
|---|---|---|---|---|---|---|---|
| Context mean | 0.081 | 0.081 | 0.000 (0.012) | 0.106 | 0.106 | 0.000 (0.012) | ignores the perturbation |
| Copy resting effect | 0.220 | −0.001 | 0.221 (0.008) | 0.228 | 0.023 | 0.205 (0.009) | |
| Additive | 0.226 | 0.007 | 0.220 (0.008) | 0.239 | 0.034 | 0.204 (0.009) | |
| Linear (ridge) | 0.463 | 0.020 | 0.443 (0.020) | 0.457 | 0.042 | 0.415 (0.020) | |
| Arc perturb_mean | 0.165 | 0.007 | 0.158 (0.009) | 0.179 | 0.026 | 0.153 (0.010) | |
| Arc pseudobulk | 0.222 | 0.006 | 0.216 (0.012) | 0.238 | 0.030 | 0.208 (0.014) | |
| STATE | 0.260 | 0.005 | 0.254 (0.016) | 0.260 | 0.020 | 0.241 (0.016) | |
| State (SE-600M + STATE) | 0.280 | 0.014 | 0.266 (0.016) | 0.297 | 0.038 | 0.259 (0.018) | |
| PIE (expression shift) | 0.228 | 0.108 | 0.120 (0.015) | 0.233 | 0.122 | 0.111 (0.018) | |
| PIE (log2 fold change) | 0.026 | 0.011 | 0.015 (0.003) | 0.029 | 0.014 | 0.015 (0.003) | |
| Stack (zero-shot) | 0.038 | 0.003 | 0.036 (0.014) | 0.023 | −0.005 | 0.028 (0.015) | 1,672 of 2,152 genes |
| biolord (adapted) | 0.030 | 0.004 | 0.027 (0.026) | 0.020 | −0.033 | 0.053 (0.028) | fit r = 0.35 |
| biolord (tutorial) | 0.005 | 0.010 | −0.006 (0.022) | −0.035 | −0.031 | −0.005 (0.023) | did not learn sgRNAs |
| CPA (default) | −0.001 | 0.003 | −0.003 (0.017) | −0.031 | 0.000 | −0.031 (0.020) | did not learn sgRNAs |
| CPA (adapted) | 0.001 | −0.002 | 0.003 (0.017) | −0.021 | −0.003 | −0.018 (0.020) | did not learn sgRNAs |

**Settings.** Context mean, copy, additive and linear are our baselines: the mean stimulated effect of the training sgRNAs; the sgRNA's resting effect; the resting effect plus the mean change between states of the training sgRNAs; and a ridge regression from 50 principal components of the resting effect to the stimulated effect, with λ chosen by five-fold cross-validation over training sgRNAs. STATE, Arc pseudobulk and Arc perturb_mean were trained with arc-state 0.11.1 (cell-load 0.10.4; `state tx train`) with the settings of the package's README example (40,000 steps, batch size 8; for STATE, a transformer backbone of 8 layers with hidden size 768 and the energy loss), on log1p CP10K expression of about 2,150 evaluation genes rather than the authors' highly variable genes, with the sgRNA as the perturbation and the state as the cell type, the few-shot split described in Methods, and the stimulated cells of 40 training sgRNAs for validation; predictions used `state tx predict --predict-only`, and effects were taken relative to the predicted non-targeting stimulated cells. State used the SE-600M foundation model to embed cells before STATE. PIE [@Verma2026] (arc-pie 1.0.0) used the authors' setting for the Replogle data (5,000 steps, effective batch size 128), with the released knowledge sources except ncbi_text, and the checkpoint with the best validation AUPRC. Stack [@Dong2026] (Stack-Large-Aligned) was run zero-shot as in its tutorial, with the resting cells of each test sgRNA as context and 1,000 stimulated control cells as target. CPA [@Lotfollahi2023] and biolord [@Piran2024] were run with their default or tutorial settings and with adapted settings. **Fit checks.** Re-predicting 200 training sgRNAs as if they were test sgRNAs, STATE correlated with their observed stimulated effects at r = 0.43–0.52 across the four donors, against −0.001 to 0.012 with other training sgRNAs, so it learned per-sgRNA effects; CPA and biolord with tutorial settings had fit correlations of at most 0.025 and did not learn per-sgRNA effects, so their A − B cannot be judged. STATE, copy, additive and linear were repeated in donors 2–4; the other models were run in donor 1 only.

## Supplementary Table S7 | Data sets

<!-- TABLE S7 -->
**Supplementary Table S7 | Data sets, their uses and accessions.**

| Data | System | Use | Source |
|---|---|---|---|
| Zhu et al. | primary human CD4+ T cells, CRISPRi, genome-wide, four donors | main analyses | authors' data release (S3, marson2025_data) |
| Replogle et al. | K562 (essential, genome-wide) and RPE1, CRISPRi | null controls, clone strength, predictor | figshare+ 20029387 |
| Nadig et al. | HepG2 and Jurkat, CRISPRi, essential genes | null controls, clone strength | GEO GSE264667 |
| Virtual Cell Challenge 2025 | H1 hESC, CRISPRi | null controls, clone strength | public training data, Arc Virtual Cell Atlas |
| Li et al. | human brain organoids, clone barcodes | clone halves, test of the relation | Zenodo 7083558; ArrayExpress E-MTAB-13148 |
| Shi et al. | mouse whole brain, AAV after birth, single nuclei | clone strength, reference within mouse | Hugging Face perturbai/wholebrain_crispr_atlas |
| Jin et al. | embryonic mouse cortex, lentivirus in utero | clone strength (library halves) | GEO GSE157977 |
| Knudsen et al. | human CAR-T cells in NSG mice | clone strength | GEO GSE302504 |
| Fagerberg et al. | mouse P14 CD8 T cells in LCMV-infected mice | clone strength | FITdb, Zenodo 22212134 |
| Zhou et al. | mouse OT-I CD8 T cells in tumour-bearing mice | clone strength | GEO GSE216800 |
| Schmidt et al. | primary human T cells, CRISPRa, two donors | clone strength | FITdb, Zenodo 22212134 |
| Cao et al. | primary human AML cells, in vitro | clone strength | GEO GSE221578 |
| Frangieh et al. | melanoma cell line, three treatments | clone strength, transfer models | scPerturb, Zenodo 13350497 |

## Supplementary Figure S1 | Three-term model: expectations and simulation scenarios

Simulated screens with two contexts, under random halves (red; both halves from one transduction) and independent replicates (blue; each half from its own transduction), for six scenarios: no clone effect; additive clone states only, with clones separate between contexts or shared by them; multiplicative clone gains only, separate or shared; and both terms, shared. **a**, Percentage of null perturbations with a positive cross-half signal; dotted line, 50%. **b**, **c**, Spearman–Brown corrected Pearson ceiling for all perturbations and for null perturbations. **d**, Context-specific share; dashed line, true value. **e**, Explained share of the predictor of the specific part from the common part (labelled shared-to-specific in the figure); dashed line, true value. **f**, Model and parameters. Additive clone states make null perturbations reproducible and raise both ceilings, under random halves only. Gains leave null perturbations unaffected. The context-specific share changes little in any scenario. Clone effects that are separate between contexts lower the explained share under random halves, and shared gains lower it across independent replicates. The scenarios without and with true structure are shown in Fig. 6e,f. Means over three seeds.

## Supplementary Figure S2 | The Spearman–Brown correction at chance level

The correction 2r/(1 + r) extrapolates a score computed between two halves to the full number of cells. For metrics whose chance level is 0.5, such as direction match and discrimination, it maps chance to 0.67 (dashed lines), so a corrected ceiling near 0.67 can mean no reproducible signal. Applying the correction to the score above chance would avoid this. Dotted line, identity.

## Supplementary Figure S3 | Unexpressed targets in K562 screens

**a**, Share of the within-screen structure of unexpressed targets (n = 178) and expressed targets (n = 651) of the K562 essential-gene library that recurs in an independent K562 screen (genome-wide screen, part 1): summed cross-screen signal divided by the geometric mean of the summed within-screen signals. Error bars, bootstrap 95% intervals over targets (1,000 resamples). Unexpressed targets recur as much as expressed ones, so their structure comes mainly from the sgRNAs or target sites, not from clones. **b**, Cross-half signal per target within the day-6 essential-gene screen, within the day-8 genome-wide screen and across the two screens. The non-recurring part of the null signal, (1 − 0.85) × 4.6, is at most about 4% of the expressed signal of 18.4 per target; because it also contains differences between the two screens, it is an upper bound on the clonal part. **c**, K562 genome-wide screen (part 1): share of unexpressed targets with a positive cross-half signal when responses are referenced to the mean of all perturbations (template removed) or to the mean of the unexpressed targets (null baseline), for random halves and for halves formed by sequencing lane. Template removal turns each null perturbation into a reproducible negative template and raises the share. Dashed line, 50%.
