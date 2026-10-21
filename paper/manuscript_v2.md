# Shared clones make random cell halves look reproducible in pooled screens

Xuanrui Huang

Department of Biomedical Informatics, Yong Loo Lin School of Medicine, National University of Singapore, Singapore

Correspondence: huangxuanrui@u.nus.edu

<!-- Draft v2, 2026-10-18. -->

---

## Abstract

Single-cell perturbation screens are used to train and evaluate models of perturbation effects. Whether a measured effect would reappear in a repeated experiment is usually judged by splitting each perturbation's cells into two random halves and asking whether the two measured effects agree. This agreement sets the score ceilings of benchmarks and selects which perturbations are kept. It assumes that the halves are independent replicates. In pooled screens, however, the cells carrying a guide descend from the few hundred cells it transduced, and expansion is uneven, so both halves contain the same clones. We tested this in two ways: with perturbations that should have no effect, whose agreement cannot come from an effect, and with the two guides against each gene, which come from separate transductions and share no clone. In a genome-scale CRISPR interference screen of primary T cells, perturbations that should have no effect agreed between random halves almost always, but between the two guides only half the time, as noise would. A T-cell clone marker shows that the halves share clones, and more cells do not remove this. Benchmark ceilings and reliability labels are therefore inflated, and inferring how an effect differs between cell states from its general form appears to work where no such relation exists. The problem is severe in primary T cells and small in cell lines. Judging reproducibility requires independent replicates.

---

## Introduction

Large perturbation atlases now measure the transcriptional responses to thousands of genetic perturbations across cell types and states<sup>1–3</sup>. They are used to train and evaluate models that predict the effects of perturbations or cellular contexts that have not been measured<sup>4,5</sup>. A central question follows: how much of a measured effect would reappear in a repeated experiment, since only that part can be predicted. Recent work has therefore argued that benchmarks should report the reproducibility ceiling of the data alongside model scores<sup>6</sup>.

The usual way to estimate reproducibility is to split each perturbation's cells into two random halves and compare the effects measured in each. At least three practices rest on this. cell-eval, the evaluation tool of the Virtual Cell Challenge, uses the agreement between random halves, corrected by the Spearman–Brown formula, as the ceiling of the data<sup>4,7–9</sup>. Per-perturbation reliability labels, proposed to clean training and evaluation sets, are based on the same agreement<sup>6</sup>. And methods that split perturbation effects into a part shared across cellular contexts and a part specific to each context use random halves to separate reproducible variation from noise, as a basis for transferring effects to new contexts<sup>10</sup>. All three assume that the cells of a perturbation are independent samples, so that the two halves differ only by sampling noise.

In pooled screens this assumption need not hold. Lentiviral libraries are delivered at low multiplicity of infection, and the cells that carry a given guide RNA (guide) at sequencing descend from a limited number of transduced cells, whose expansion is often very uneven<sup>11–13</sup>. The descendants of one transduced cell form a clone, and a clone has its own expression state. We call a group of cells that shares such a state, such as a clone, a unit, and halves of a perturbation's cells that share no unit independent units. Random halves therefore contain the same clones, and their agreement may reflect the state of those clones rather than the effect of the perturbation.

A similar statistical problem is well known in single-cell differential expression, where treating cells from the same individual as independent inflates false discoveries<sup>14,15</sup>. It has also been noted that random halves measure sampling precision rather than experimental reproducibility<sup>16</sup>, and clone barcodes have been used to follow clonal differences in pooled screens with phenotypic readouts<sup>11</sup>. To our knowledge, however, how large the effect is in screens with a transcriptomic readout, where it comes from and what it does to the three practices above have not been measured.

Here we measure it in screens that contain independent units for the same perturbation and context (Fig. 1): the two guides against each gene, which are independent transductions; different donors; and independent screens of the same cell line. We use perturbations that should have no effect as negative controls, and in T cells a natural clone marker to observe the clones directly. We find that the agreement between random halves largely reflects shared clones; that the size of this component differs between screens and between donors; and that it inflates benchmark ceilings and reliability labels and produces structure across cell states that independent units do not reproduce.

---

## Results

### Null perturbations look reproducible under random halves

We first asked whether a perturbation with no expected effect looks reproducible when its cells are split at random. The genome-scale CRISPR interference (CRISPRi) screen of primary human CD4+ T cells<sup>3</sup> profiles cells of four donors in a resting and a stimulated state, with two guides against each gene. Its probe-based readout also measures the target genes, except those outside its probe panel, and some of them are not expressed in T cells (1,239 targets in donor 1 and 872 in donor 2; Methods). Knocking down a gene that is not expressed cannot change the transcriptome through that gene; we call such perturbations null perturbations.

For each target we compared its cells with control cells in each of two halves, which gives two effects, and computed their cross-half signal: the gene-wise product of the two effects, summed over genes. It estimates the energy, the sum of squares over genes, of the part of the effect that the two halves share. When two halves differ only by sampling noise, their effects are unrelated, and the cross-half signal of a null perturbation is positive about half the time; our pipeline returns this in simulated screens (Supplementary Fig. S1a).

Under random halves, 96–99% of the unexpressed targets had a positive cross-half signal, in both donors and both states (Fig. 2a,b). With one guide in each half (guide halves) the share was 49–54%, and with donors 1–2 against donors 3–4, from the released data of all four donors, it was 53–55%, close to chance (Fig. 2b). Random halves therefore report a reproducible response for perturbations that should have none; independent units do not.

A guide against an unexpressed target can still act elsewhere, for example at a neighboring promoter. Non-targeting guides, which have no target at all, exclude this, and they behaved in the same way: 97–99% were positive under random halves, against 49–50% for halves taken from two different non-targeting guides (Fig. 2c; Supplementary Table S1). Non-targeting cells are also the controls that benchmark tools split at random when they compute ceilings.

### The shared state points to clones

What do random halves share that guide halves do not? In T cells a natural marker answers this directly. During T-cell development the T-cell receptor (TCR) β gene is rearranged once, and depending on the rearrangement the cell uses one of two constant regions, TRBC1 or TRBC2; it passes this choice to all its descendants. The cells of one clone therefore share a TRBC1 fraction (Fig. 3a), and the readout of the CD4 screen measures both constant regions. If a perturbation's cells came from many clones, their TRBC1 fraction would be close to that of the whole population; if they came from few, it would deviate from it. Scaled as a z-score against random assignment of cells, the deviation would have variance 1 under the first possibility; it was 34 in donor 1 and 81 in donor 2 (Fig. 3b). It was nearly identical between the random halves of a perturbation (r = 0.93 and 0.96) and unrelated between its guide halves (r ≈ 0; Fig. 3c): random halves contain the same clones, and the two guides contain different ones. The size of the deviation implies that a perturbation's cells behave as if they came from at most about 7 (donor 1) and 4 (donor 2) equally sized clones (Methods), although of the order of 500 cells per perturbation were transduced<sup>3</sup>. A few clones expanded far more than the rest.

The marker also tracks the expression signal. Among unexpressed targets, the more a perturbation's TRBC1 fraction deviated, the larger its random-half signal (Spearman correlation 0.21 and 0.24, TCR genes excluded; Fig. 3d). A marker with two states captures only part of the clonal composition, so a modest correlation is what clones would produce. The guide sequence, by contrast, does not explain the null signal. If it did, the same guide would agree with itself in independent donors. Pairing donors 1–2 with donors 3–4, unexpressed targets were positive in 57% (resting) and 54% (stimulated) of cases with the same guide, against 51% with different guides (Supplementary Table S2). The shared component is therefore a property of the transduced cells, not of the guide.

This argument needs null perturbations whose agreement cannot come from the guide, and not every library provides them. Libraries selected on a growth phenotype, such as the K562 essential-gene library<sup>1</sup>, include unexpressed targets whose guides affected growth, most likely by acting at a neighboring promoter<sup>17</sup>. Their effects are then real and recur across independent screens: in two separate K562 screens, 85% of the structure of unexpressed targets recurred (Supplementary Fig. S2). In such libraries, unexpressed targets cannot separate clones from guide effects; non-targeting guides and independent transductions of the same guide can.

### More cells do not remove the shared component

If the shared component were sampling noise, measuring more cells would shrink it. To test this, we also split each perturbation's cells into quarters and compared the two depths on the same targets (Fig. 3e). The null cross-half signal under random splits was the same at a quarter and at half of the cells, and under guide splits it was near zero at both. The between-half Pearson correlation of these null perturbations nevertheless rose with depth, from 0.14 to 0.24 in donor 1. This is what a clonal component predicts: with more cells, sampling noise falls while the clonal component stays, so the correlation between random halves approaches the agreement within the same clones. The Spearman–Brown correction, which extrapolates a split-half correlation to the full number of cells, predicted this rise correctly (0.25 against an observed 0.24). What it extrapolates, however, is the agreement within one transduced population, not across experiments.

### How large the shared component is, and where

Most of what random halves call reproducible in the CD4 screen is shared only within one transduction. For expressed targets, whose knockdown has real effects, the reproducible energy fell by about 80% in donor 1 and 90% in donor 2 when guide halves replaced random halves. Part of this loss could come from differences between the two guides rather than from clones. To separate the two, we compared pairings across the two donor groups, which share no clones (Supplementary Table S2; Methods). The same guide in donors 1–2 and in donors 3–4 shares the effect of the knockdown and the guide's own effects; two different guides share only the effect of the knockdown. The difference between the two pairings measures the guide's own effects. Transferred to single cells, the random-half signal of expressed targets splits into three parts: 10–19% is reproduced by the other guide, 8–15% is specific to the guide, and 66–81% is shared only by cells of one transduction (Fig. 3f). Donor 2, whose perturbations were more clonal by the TRBC marker, had the larger share.

Across screens, the size of the shared component differs widely. Non-targeting constructs allow a comparison, because they have no effect in any screen. In five cell-line CRISPRi screens<sup>1,2</sup> and in the H1 embryonic stem cell data of the 2025 Virtual Cell Challenge<sup>4</sup>, they were positive somewhat more often than chance under random halves (60–76%), but their random-half signal was only 0.6–7% of that of a typical targeting perturbation. In the CD4 screen it exceeded that of a typical targeting perturbation (1.2–1.3 times; Fig. 2c). The inflation is therefore large in the primary T-cell screen and small in the cell lines and in the H1 benchmark data. It also differs between donors of the same screen, so no single correction factor applies.

### Ceilings and reliability labels are inflated

Benchmarks use random halves to set the score that a perfect model could reach. cell-eval (version 0.8.2)<sup>7</sup>, the evaluation tool of the Virtual Cell Challenge, computes this ceiling by comparing two random halves of each perturbation on each of its metrics, and corrects the result from half to full depth with the Spearman–Brown formula. We ran the tool on 300 perturbations from resting cells of donor 1: 150 unexpressed targets and 150 strength-matched knockdowns, in which both guides knocked down the target to a similar extent, with equal numbers of cells per guide. We compared its random-half ceiling with the same metrics computed between the two guides of each perturbation.

The random-half ceiling exceeded the guide ceiling on every metric (Fig. 4a; Supplementary Table S3). For the Pearson correlation of effects it was 0.32 against 0.06 over all 300 perturbations, and 0.38 against 0.13 for the strength-matched knockdowns, so the gap is not explained by one guide being weaker. Perturbations with no expected effect received high ceilings, for example 0.26 for the Pearson correlation and 0.97 for discrimination between perturbations, while their guide ceilings stayed near the corrected chance level. For metrics whose chance level is 0.5, the Spearman–Brown correction itself turns chance into 0.67 (Supplementary Fig. S3); this is a separate problem, which applying the correction to the score above chance would avoid.

Reliability labels behave the same way. Under a recently proposed split-half rule<sup>6</sup>, a perturbation is labelled reliable if the corrected correlation between its halves is at least 0.5. With random halves, 18% of unexpressed and 24% of expressed targets in donor 1 were labelled reliable, against 0.6% and 2.2% with guide halves. In donor 2, whose perturbations were more clonal, the shares were 42% and 52% with random halves and 0.9% and 2.3% with guide halves (Fig. 4b).

### Cross-context structure comes from shared units

Perturbation effects measured in several cellular contexts, such as cell types or states, are often split into a part shared across contexts and a part specific to each context, with random halves used to separate reproducible variation from noise<sup>10</sup>. For two contexts, such as resting and stimulated cells, in which a perturbation has effects x_R and x_S, the shared part is s = (x_R + x_S)/2 and the context-specific part is d = (x_R − x_S)/2. Models that transfer perturbation effects to unmeasured contexts rest on such splits<sup>10</sup>, and current benchmarks increasingly test this task<sup>4</sup>. If contexts share units, a unit's deviation enters both parts (Fig. 5a). This is the case in the CD4 screen: resting and stimulated cells were split from one culture after expansion<sup>3</sup>, and a perturbation's TRBC1 deviation correlated between the two states nearly as strongly as between random halves (r = 0.90 and 0.94; Fig. 5b), so the same clones are present in both states.

The fraction of the reproducible energy that is context-specific was about 0.3 under both random and guide halves, but the energies behind it were not the same: both parts were 8–9 times larger under random halves (Fig. 5c). A stable fraction therefore does not show that an estimate is free of shared-unit effects.

We then asked whether the shared part predicts the context-specific part, as transfer models assume. A linear predictor f was trained on many perturbations to map the shared part from one half of the cells to the context-specific part from the other half, d_B ≈ f(s_A), and was scored on perturbations it had not seen. Taking the two parts from different halves is meant to keep sampling noise from linking them, so that only reproducible structure can be predicted. Random halves keep the noise apart but not the clones: the halves contain the same clones (Figs. 2 and 3), the clones are present in both states (Fig. 5b), and a clone's deviation is related but not identical in the two states. A clone that shifts expression by c in resting and by qc in stimulated cells, with q ≠ 1, adds (1 + q)c/2 to s and (1 − q)c/2 to d; when halves A and B contain the same clone, f can predict d_B from s_A through c alone (Fig. 5a). In CD4 the predictor explained about half of the reproducible context-specific energy under random halves and about none under guide halves, including strength-matched knockdowns and the same guide in independent donors (Fig. 5d; Supplementary Table S2). Under random halves it also carried the deviations of unexpressed targets into its predictions as much as those of expressed targets (Supplementary Table S4): what it learned was the state of the units, not the response to the perturbation. For K562 and RPE1, which were screened separately and share no units<sup>1</sup>, taking the K562 part of the shared response from an independent K562 screen left the predictor unchanged (Fig. 5d).

A model with three terms reproduces these patterns (Fig. 5e, Supplementary Fig. S1; Methods). Each cell carries an additive state of its unit, a multiplicative gain of its unit on the perturbation's response, and sampling noise. The additive term makes null perturbations reproducible and inflates ceilings; neither unit term shrinks with more cells; and when contexts share units, a predictor appears where none exists (an explained share of 0.10 when the true value is 0). The model also shows the reverse risk: large unit deviations along the directions of the response can hide a real cross-context structure across independent units. At the ratios of unit to signal energy found in the CD4 screen (4.0 and 7.5; Methods), a true explained share of 0.70 falls to 0.10 and 0.03. The absence of a predictor across guides therefore shows only that it cannot be used across independent units, not that no cross-context structure exists.

Box 1 summarizes how to check whether random halves can be trusted.

<!-- BOX START -->
**Box 1 | Checking whether random halves can be trusted**

*Do the two halves share clones?* Halves formed by different guides or donors do not. They can be used as replicates, and a null control should still be reported.

*If they do* (random halves of the same cells), reliability labels, ceilings and cross-context structure will be inflated. The inflation is likely to be large in primary cells, after strong bottlenecks and when the same clones are measured in several conditions. It was small in five cell-line screens and in the H1 data of the 2025 Virtual Cell Challenge, where the random-half signal of non-targeting constructs was 0.6–7% of that of a typical targeting perturbation.

*In every case*, include perturbations that should do nothing, such as non-targeting guides or unexpressed targets. Compute their effects relative to their own mean rather than to the mean of all perturbations, and do not take them from libraries selected on a growth phenotype.
<!-- BOX END -->

---

## Discussion

Random halves of a perturbation's cells answer a narrower question than the one they are usually asked. They estimate how precisely the response of a given set of clones has been sampled, not how well an independent experiment would reproduce it. Random halves remain the right tool when the question is sampling precision within one experiment, for example whether more cells of the same population would sharpen an estimate; they are the wrong tool when the question is what an independent experiment would show. The same applies wherever a group of cells shares a state that random halves cannot separate, for example the cell lines pooled in one well of a multiplexed drug screen. Independent transductions of the same guide answer that question. Different guides answer a third: how reproducible the effect of knocking down the gene is, a question in which differences between guides are part of the answer. The difference is not a small correction. In a primary T-cell screen, most of what random halves call reproducible was shared only within one transduction, most likely by clones, and perturbations with no expected effect received high ceilings and reliability labels.

The statistical problem is known from single-cell differential expression<sup>14,15</sup>, and that random halves measure sampling precision has been noted before<sup>16</sup>. Our results add three points. First, we compared random halves with independent units within the same data, together with null perturbations and a direct clone marker; this points to clones and measures the size of the shared component. Second, the component does not shrink when more cells of the same population are sequenced, so more cells, the remedy implied by sampling models of reliability<sup>6</sup>, do not remove it; with more cells, the agreement between random halves approaches the agreement within clones. Third, shared units affect not only reproducibility but also structure across contexts. When the same clone is present in several contexts, its deviation enters the shared and the context-specific responses together, and a predictor that links them can hold within the same units and fail across independent units. This matters for models that transfer perturbation effects across cellular contexts, a task that current benchmarks increasingly target<sup>4</sup>.

The original study used appropriate replicate units for its own analyses, modeling donors as biological replicates<sup>3</sup>. The issue arises when random halves are used downstream, in benchmark ceilings, reliability labels, denoising and cross-context modeling. In the H1 data of the 2025 Virtual Cell Challenge the shared component was small (Fig. 2c), so ceilings computed on those data are little affected; the inflation matters most in primary T cells, where each perturbation's cells come from few clones.

Our analysis has limitations. The same-data comparisons in CD4 used two donors. Guide halves differ in the guide as well as in clones. We estimated the guide-specific part by pairing the same guide across donors and matched knockdown strength, but that comparison crosses donors and uses pseudobulk; independent transductions of the same guide within one donor would be a cleaner unit. Targets outside the probe panel were excluded from the null set. Finally, the model reproduces the observed patterns but was not fitted to the data.

Future screens can build in independent units. Two or more guides per gene are already standard, and clone barcodes add little cost<sup>11</sup>; separate transductions and donors cost more but give the most direct test. For existing data, a set of null perturbations is the minimum check.

---

## Methods

### Data

**Primary human CD4+ T cells.** The genome-scale CRISPRi Perturb-seq screen of Zhu et al.<sup>3</sup> profiles cells of four donors, resting and stimulated (48 h), with a probe-based readout and two guides per target. We used the single-cell files of donors 1 and 2, without cells flagged as low quality or lacking an assigned guide, and the released pseudobulk counts of all four donors.

**Other screens.** We used the raw single cells of the K562 essential-gene (day 6), K562 genome-wide (day 8, part 1) and RPE1 essential-gene screens of Replogle et al.<sup>1</sup>, normalized to 10,000 counts over the genes common to the three screens, log1p transformed and corrected additively for the control mean of each GEM group (one 10x Genomics channel); the Jurkat and HepG2 essential-gene screens of Nadig et al.<sup>2</sup> (Fig. 2c and Supplementary Table S4); and the H1 hESC training data of the 2025 Virtual Cell Challenge<sup>4</sup> (221,273 cells, 150 targets and 31 non-targeting constructs, 10x Flex; Fig. 2c only).

### Effects and gene filters

The effect of a perturbation in a set of its cells is their mean log1p CP10K expression minus that of control cells, the non-targeting guides. Controls were split like the perturbed cells in CD4; in the Replogle data, both halves used all controls, which cancels after centering. The target gene was set to zero in its own effect, except in the ceiling analyses and the template check, which keep it as the benchmark tools do. Unless stated otherwise, effects were centered over perturbations within each half and context.

Genes were kept if their control mean log1p CP10K was at least 0.1: in CD4 single cells, in both states for analyses that combine states and in the analyzed state otherwise; in the CD4 pseudobulk, in every donor and state; and in the Replogle data, in all three screens (the 2,000 genes of highest control mean for the template check). cell-eval used genes with a resting control mean of at least 0.01, together with the target genes. TCR genes (TRAC, TRBC, TRAV, TRBV, TRGC and TRDC) were excluded from the CD4 null-share, depth, clonal-energy and TRBC-link analyses.

### Splits

In CD4, each cell carries two random bits, r1 and r2 (seed 21). Random halves split on r1; the bootstrap and the reliability labels also used r2 and r1 ⊕ r2. Random quarters are (r1, r2) = (0, 0) against (1, 0). Guide halves split on the parity of the guide number in the guide name, so that each guide falls on the same side in every file, and controls were split likewise by non-targeting guide; guide quarters add r2 = 0. Donor halves sum the pseudobulk counts over donors 1–2 or 3–4, both guides and both sequencing runs (log1p CP10K of the summed counts minus that of the summed controls; at least 40 cells per half and state). In the Replogle data, the cells of each perturbation were permuted (seeds 0, 1 and 2 for the three screens) and cut into equal halves.

### Unexpressed and expressed targets

Unexpressed targets are probed targets with control mean log1p CP10K below 0.01; expressed targets have at least 0.5. They are an operational null: their knockdown cannot act through the target gene, although a guide can act elsewhere, for example at a neighboring promoter. Analyses that combine states (Fig. 2a,b, the context-specific share, the predictor and the reliability labels) apply the threshold in both states; analyses within one state (depth, clonal energy and ceilings) apply it in that state. In the Replogle data, which release only expressed genes, targets absent from the gene list counted as unexpressed; the analysis across K562 screens used the 178 unexpressed in both.

In CD4, perturbations needed at least 10 cells on each side of both the random and the guide halves in both states for Fig. 2a,b (1,239 unexpressed targets in donor 1 and 872 in donor 2), at least 8 on each side of every quarter design for the depth and clone tests, and at least 10 on each side of all four halvings in both states for the predictor, the context-specific share and the reliability labels.

### Non-targeting guides

Each non-targeting guide of CD4 donors 1 and 2 was treated as a perturbation, with the same random bit per cell and effects referenced to the mean over all kept non-targeting guides in the half. Genes had a mean of at least 0.1 over non-targeting cells (TCR genes excluded), and guides needed at least 10 cells per half. For halves from different guides, guides were paired at random and the first half of one was compared with the second half of the other, at the same depth; the share of positive pairs is the mean over 50 pairings. The TRBC marker was computed as for targeting guides. The same analysis was run in the Replogle, Nadig and H1 data, with one pseudo-perturbation per non-targeting construct (sgID_AB, or guide_id in the H1 data). The ratio in Fig. 2c is the mean random-half cross-half signal of non-targeting constructs divided by that of targeting perturbations, which were grouped by target gene, referenced to all non-targeting cells of the half and centered, with the same genes and at least 10 cells per half (CD4: resting cells, TCR genes excluded).

### Cross-half signal and null perturbations

The cross-half signal of a perturbation is Σ_g x_A,g x_B,g, the gene-wise product of its two half-effects summed over genes; with independent noise in the two halves, it is an unbiased estimate of the energy of the reproducible effect. For null perturbations, every perturbation was referenced to the null baseline, the mean effect of the unexpressed targets in each half. Referencing to the mean of all perturbations instead (template removal) makes a null perturbation equal to the negative template in both halves; in the K562 genome-wide screen (random bit per cell, seed 9; at least 15 cells per half) it raised the share of reproducible null targets from 83% to 89%.

In CD4, the share of the random-half signal of expressed targets that the other guide does not reproduce is 1 minus the ratio of the guide-half to the random-half cross-half signal, for the same state and targets. To split this share into guide-specific and transduction-level parts, we formed four pseudobulk effects per target and state from all four donors: guide rank 0 or 1 in donors 1–2 or 3–4 (both sequencing runs, minus the non-targeting controls of the same donors and state; at least 20 cells each). The same-guide design pairs each guide across the two donor groups, and the different-guide design pairs guide 0 with guide 1, on the same cells. Their difference, divided by the different-guide signal, estimates guide-specific energy relative to the signal shared across guides (0.80 resting, 0.48 stimulated). Multiplied by the single-cell guide-half signal, it gives the guide-specific part of the random-half signal, and the remainder above the guide-half signal is the transduction-level part. The single-cell guide-half signal also contains donor-specific responses, so this transfer tends to overestimate the guide-specific part.

For K562, the share of the unexpressed-target structure that recurs across the two independent screens is the summed cross-screen signal divided by the geometric mean of the summed within-screen signals (0.85). The unit-to-signal energy ratio is the null cross-half signal under random halves divided by the expressed cross-half signal under guide halves, as in the simulation.

### TRBC clone marker

We used cells with at least two TRBC1 plus TRBC2 counts (1.75 million resting and 1.58 million stimulated cells of donor 1) and their TRBC1 fraction f = TRBC1 / (TRBC1 + TRBC2). For a perturbation with n such cells, z = (mean f − f0) / (sd_f / √n), where f0 and sd_f are the mean and standard deviation over all cells, so that var(z) ≈ 1 under random assignment. var(z) used perturbations with at least 30 cells, and its permutation value is the mean over five shuffles of perturbation labels. Reproducibility of the deviation (mean f − f0) is its between-half Pearson correlation (at least 15 cells per half). Clonal skew is |z|; its link to the cross-half signal (Fig. 3d) used the unexpressed targets of Fig. 2a,b with at least 30 cells carrying TRBC counts (1,234 in donor 1 and 864 in donor 2), with quartiles over these targets. With the per-state definition of unexpressed targets, the Spearman correlations were 0.229 and 0.227 (resting) and 0.159 and 0.160 (stimulated) in donors 1 and 2.

If the n cells of a perturbation come from K equally sized clones with intra-clone correlation ICC of the TRBC1 fraction, E[z²] ≈ 1 + ICC(n/K − 1), so K ≤ n/z² for ICC ≤ 1. We report the median of n/z² over resting perturbations with z² > 1 (85–91% of them), with all guides of a perturbation pooled. The screen transduced 30 million cells per donor at a multiplicity of infection of about 0.2 with 23,256 guides<sup>3</sup>, about 250 transduced cells per guide and 500 per perturbation. Between states, the TRBC1 deviation in units of the cell-level standard deviation, (mean f − f0)/sd_f = z/√n, was correlated over the perturbations of the main analysis with at least 30 cells carrying TRBC counts in both states (6,582 in donor 1 and 4,963 in donor 2).

### Depth

Half and quarter designs were compared on the same targets. The between-half Pearson correlation is the mean over targets of the per-target correlation between half-effects, with the Spearman–Brown value 2r/(1 + r) applied to that mean. In the model below, the unit term scales with 1/K_eff, the inverse effective number of clones, not with the number of cells sampled.

### Shared and context-specific responses

For two contexts X and Y, the shared response is (X + Y)/2 and the context-specific response is (X − Y)/2, after centering each context over perturbations. The context-specific share is the cross-half signal of the context-specific part divided by the sum of the shared and context-specific cross-half signals.

### Shared-to-specific predictor

In the main implementation (CD4 donors 1 and 2), the features of a perturbation are its half-A shared response projected on the top k right singular vectors of the training perturbations' shared responses (k = 20). A ridge regression maps the features to the context-specific response in half B, with λ chosen from {10⁻³, 10⁻², …, 10²} × tr(ZᵀZ)/k by five-fold inner cross-validation, and perturbations were split into five outer folds. The explained share is (Σ‖γ_B‖² − Σ‖γ_B − F‖²) / Σ γ_A · γ_B, where F is the held-out prediction and the denominator is the reproducible energy of the target. A shuffled control, permuting the features across perturbations, gave values near zero throughout. Centering preceded cross-validation and was the same for every design; principal directions and λ were estimated within the training folds.

A second implementation, used for the CD4 robustness checks, the strength-matched subsets and K562–RPE1, uses two outer and two inner folds, λ up to 10³ and k = 10 and 20 (5–50 in the robustness checks; Supplementary Table S5). For K562–RPE1, perturbations needed at least 20 cells in all three screens, p53-pathway and cell-cycle genes and TP53 were excluded, and the independent shared response was (other K562 screen + RPE1 half A)/2. Strength-matched knockdowns in CD4 are expressed targets (at least 0.5 in both states) whose two guides both knocked down the target, the weaker reaching a Δlog1p of at most −0.25 and at least 70% of the stronger, in both states; effects were re-centered within the subset (615 targets in donor 1 and 470 in donor 2). The text and Fig. 5d give k = 20; with k = 10 the explained shares were 0.20 and 0.29. The two implementations gave the same conclusions (in CD4 donor 1, 0.446 and 0.469 under random halves and at most 0.004 under guide halves).

### Ceilings

We ran cell-eval version 0.8.2 (PyPI, 27 July 2026)<sup>7</sup> in an isolated environment on resting cells of donor 1, with at least 40 cells per guide: 150 unexpressed targets drawn at random and the 150 targets with the largest mean knockdown among those meeting the strength-matching rule in resting cells. Each perturbation was capped at 400 cells, 15,000 control cells were drawn, and cells were downsampled to equal numbers per guide (median 62 per side); perturbations with fewer than 10 cells per side were dropped. Both designs used the same cells. The random-half ceiling is the tool's `compute_ceiling` with seeds 0–2. The guide ceiling runs the tool's `compute` with one guide's cells as "real" and the other's as "prediction", in both directions, and applies the tool's Spearman–Brown function to the mean over perturbations; our aggregation reproduced the tool's own aggregate exactly. Results were summarized over perturbations present in all runs.

### Reliability labels

Following ref. 6, ρ is the Spearman–Brown value of the median between-half Pearson correlation, and perturbations with ρ ≥ 0.5 are labelled reliable. For random halves the median was taken over the three halvings r1, r2 and r1 ⊕ r2; guide halves have a single split. Undefined values were counted as not reliable. Resting cells of donor 1 were used.

### Model and simulation

A cell of perturbation p in context c from unit k has expression μ_pc(1 + g_k) + a_k + e, with unit state a_k ~ N(0, σ_a²I), unit gain g_k ~ N(0, σ_g²) and noise e ~ N(0, σ_e²I). For random halves, E[x_A · x_B] = |μ|²(1 + σ_g²/K_eff) + Gσ_a²/K_eff, with K_eff = 1/Σ_k w_k² for unit weights w_k; for independent halves it is |μ|². The simulation used G = 1,000 genes, 800 perturbations (20% null), true responses of rank 10, log-normal effect sizes s_p (log-mean log 0.3, s.d. 0.5) and context-specific parts s_p(0.5 z_p R W₂ + 0.3 d_p). Each perturbation and transduction has K = 6 units with Dirichlet(0.5) weights, σ_a = 0.2, σ_g = 0.6 and σ_e = 1, and each half draws 60 cells from the units (multinomial). Units are shared between contexts (additive part scaled by ρ = 0.5 in the second) or separate. In the fabrication block, the effect scale was 0.1, there was no true shared-to-specific structure and unit states lay along the shared response directions. In the attenuation block, the true explained share was about 0.70 and σ_a was 0.2, 0.6, 1.2, 1.5, 2.0 or 2.5, isotropic or along the shared response directions; σ_a of 1.5 and 2.0 give unit-to-signal energy ratios of 4.1–4.2 and 7.2–7.6, the range of the CD4 screen. The predictor was a two-fold version of the one above (k = 10, λ from 10⁻² to 10³), and results are means over three seeds.

### Statistics

Confidence intervals are 2.5th–97.5th bootstrap percentiles over 1,000 resamples of perturbations, with predictions held fixed; Supplementary Table S4 used 500 resamples. They describe variation over perturbations within one experiment; agreement across the available independent units (two donors in the CD4 single-cell data, four in its pseudobulk and two K562 screens) is reported separately. In CD4, the three random halvings share cells and are not independent resamples; they gave nearly identical estimates (an explained share of 0.468–0.471 in donor 1). Shares of positive null perturbations from separate extractions, halvings and minimum-cell filters agreed within 1.5 percentage points.

### Software

Analyses ran in Python 3.11 with NumPy 2.4, SciPy 1.17, h5py 3.16 and PyArrow 20; cell-eval 0.8.2 ran in its own environment (Python 3.11, AnnData 0.12); figures were drawn with Matplotlib 3.10.

### Data availability

All data are public. The CD4 Perturb-seq data are available at https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data, the data of Replogle et al. as figshare+ article 20029387 (https://doi.org/10.25452/figshare.plus.20029387), the data of Nadig et al. under GEO accession GSE264667, and the H1 data of the 2025 Virtual Cell Challenge through the Arc Virtual Cell Atlas (https://console.cloud.google.com/marketplace/product/bigquery-public-data/arc-institute).

### Code availability

The analysis code, the result files behind every figure and table, and a record of the reproduction of all reported numbers from the raw data are available at https://github.com/xuanruihuang/random-halves-not-replicates under the MIT licence; the version used here (v1.0.0) is archived at Zenodo<sup>18</sup>.

---

## Figure legends

**Figure 1 | Random halves share clones; independent units do not.** **a**, In a pooled screen, the cells carrying a guide descend from a few transduced cells; the descendants of one transduced cell form a clone (dashed rings, one colour per clone), and uneven expansion lets a few clones contribute most of the cells. Random halves of a perturbation's cells contain the same clones and agree whether or not the perturbation has an effect; halves formed by guide or donor share no clone and agree only through the perturbation's effect. **b**, Uses of random halves, each comparing half A with half B of the same cells: noise ceilings and reliability labels, decomposition into shared and context-specific parts, and evaluation of predictors of the context-specific part from the shared part. **c**, Datasets used here, with their independent units and null controls; the two K562 screens are the essential-gene and the genome-wide screen. Asterisk: in phenotype-selected libraries, the structure of unexpressed targets is guide-level rather than clonal. NT, non-targeting.

**Figure 2 | Null perturbations look reproducible under random halves.** **a**, Cross-half signal of unexpressed targets in resting CD4+ T cells of donor 1 (n = 1,239) under random (red) and guide (blue) halves; values clipped to [−15, 30]. **b**, Percentage of unexpressed targets with a positive cross-half signal in donors 1 (n = 1,239) and 2 (n = 872), under random halves (red) and guide or donor halves (blue; donors 1–2 against 3–4, n = 1,479). Error bars, bootstrap 95% intervals over targets (1,000 resamples); dashed line, 50%. **c**, Non-targeting guides in resting CD4+ T cells (donors 1 and 2, n = 739 and 578), in five cell-line CRISPRi screens<sup>1,2</sup> (K562 essential, 97; K562 GW, 507; RPE1, 113; Jurkat, 128; HepG2, 112) and in the H1 hESC data of the 2025 Virtual Cell Challenge<sup>4</sup> (31). Left, percentage with a positive cross-half signal under random halves (red; error bars, bootstrap 95% intervals over guides) and between two different guides at the same depth (blue; mean and range over 50 pairings). Right, mean cross-half signal of non-targeting guides under random halves, relative to that of targeting perturbations of the same screen; dashed line, 1. Numbers in Supplementary Table S1. GW, genome-wide.

**Figure 3 | Random halves share clones, which carry most of the random-half signal at any depth.** **a**, The TCR β rearrangement of each T cell fixes whether it uses the TRBC1 or the TRBC2 constant region, and its clone inherits the choice (dashed rings, one clone each). The cells of a perturbation that come from many clones have a TRBC1 share near the population average; cells that come from few clones deviate from it, and random halves of them share the deviation. **b**, z-score of each perturbation's mean TRBC1 fraction relative to random assignment of cells, observed (orange) and with perturbation labels permuted (grey); resting cells of donor 1, perturbations with at least 30 cells carrying TRBC counts (n = 6,610). **c**, Deviation of the TRBC1 fraction in half A against half B under random (left) and guide (right) halves; donor 1, resting; perturbations with at least 15 cells carrying TRBC counts in each half (n = 6,568 and 6,072). **d**, Random-half cross-half signal of unexpressed targets by quartile of clonal skew, the absolute z-score of panel b, TCR genes excluded; circles, donor 1 (n = 1,234); squares, donor 2 (n = 864); error bars, s.e.m. **e**, Cross-half signal (left) and between-half Pearson correlation (right) of unexpressed targets at a quarter and half of each perturbation's cells, under random (red) and guide (blue) splits; circles, donor 1 (n = 909 resting, 856 stimulated); squares, donor 2 (n = 520 and 523); solid lines, resting; dashed lines, stimulated. **f**, Composition of the random-half signal of expressed targets: shared only within one transduction (red), specific to the guide (pink) and reproduced by the other guide (blue); n = 787 (donor 1, resting), 847 (donor 1, stimulated) and 515 (donor 2, resting) expressed targets (Methods).

**Figure 4 | Random halves inflate benchmark ceilings and reliability labels.** **a**, cell-eval 0.8.2 ceilings (Spearman–Brown corrected) under random halves (red; mean of three seeds) and guide halves (blue; mean of both directions) for 300 perturbations of resting cells of donor 1, median 62 cells per side. Asterisks, metrics with chance level 0.5; dotted line, 0.67, the corrected chance level (Supplementary Fig. S3); open marker, no ceiling defined (split-half correlation at most 0), plotted at 0. Ranges in Supplementary Table S3. **b**, Percentage of perturbations labelled reliable (ρ ≥ 0.5) under the split-half rule of ref. 6, in resting cells of donor 1 (all, n = 6,615; unexpressed, 1,235; expressed, 1,097) and donor 2 (5,003; 869; 859); error bars, bootstrap 95% intervals over perturbations. DE, differential expression.

**Figure 5 | Shared units carry cross-context structure that independent units do not reproduce.** **a**, Left, one clone measured in two states and split into random halves: its deviation is present in all four cell sets and is related but not identical between the states, so it enters both the shared part computed from half A and the context-specific part computed from half B. Right, separately screened lines share no units, and each clone exists in one context only. **b**, TRBC1 deviation of each perturbation in resting against stimulated cells of donor 1, in units of the cell-level standard deviation; perturbations of the main analysis with at least 30 cells carrying TRBC counts in both states (n = 6,582; in donor 2, n = 4,963 and r = 0.94). **c**, Reproducible energy per target of the shared and context-specific parts under random (red) and guide (blue) halves; CD4 donor 1 (n = 6,615); numbers above bars, context-specific share. **d**, Shared-to-specific predictor (k = 20) under random (red) and guide (blue) halves in CD4 (donor 1, n = 6,615; donor 2, 5,003; strength-matched knockdowns, 615 and 470; error bars, bootstrap 95% intervals over targets, shown for the analyses of all targets), and for K562 and RPE1 screened separately, with the K562 part of the shared response from the same (red) or an independent screen (blue). **e**, Three-term model, means of three seeds. Left, no true structure: explained share under random halves (red) and independent units (blue). Right, true explained share 0.70 (dashed line): explained share across independent units against the unit-to-signal energy ratio; grey band, the two CD4 donors.

---

## Ethics

This study analyzed only publicly released, de-identified data.

## Author contributions

X.H. conceived the study, performed the analyses and wrote the manuscript.

## Competing interests

The author declares no competing interests.

---

## Supplementary information

- **Supplementary Table S1.** Non-targeting guides as null perturbations: random halves against different guides, and the TRBC clone marker.
- **Supplementary Table S2.** Same guide against different guides across independent donors: guide-specific effects and the shared-to-specific predictor.
- **Supplementary Table S3.** cell-eval 0.8.2 ceilings for all metrics (all, unexpressed and strong subsets; corrected and uncorrected; ranges over seeds and directions).
- **Supplementary Table S4.** Structure of unexpressed targets carried by the shared-to-specific predictor, after removing the template direction: CD4 donors 1 and 2 (random and guide halves), K562, RPE1, Jurkat and HepG2.
- **Supplementary Table S5.** Robustness of the CD4 results to target set, gene set, guide balance and number of components.
- **Supplementary Figure S1.** Three-term model: expectations and simulation scenarios.
- **Supplementary Figure S2.** Unexpressed targets in K562 screens: recurrence across two independent screens, energy per target, and the template artifact in the genome-wide screen.
- **Supplementary Figure S3.** The Spearman–Brown correction at chance level.

---

## References

1. Replogle, J. M. et al. Mapping information-rich genotype–phenotype landscapes with genome-scale Perturb-seq. *Cell* **185**, 2559–2575 (2022).
2. Nadig, A. et al. Transcriptome-wide analysis of differential expression in perturbation atlases. *Nat. Genet.* **57**, 1228–1237 (2025).
3. Zhu, R. et al. Genome-scale perturb-seq in primary human CD4+ T cells maps context-specific regulators of T cell programs and human immune traits. *Cell* https://doi.org/10.1016/j.cell.2026.08.002 (2026).
4. Roohani, Y. H. et al. Virtual Cell Challenge: toward a Turing test for the virtual cell. *Cell* **188**, 3370–3374 (2025).
5. Ahlmann-Eltze, C., Huber, W. & Anders, S. Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines. *Nat. Methods* **22**, 1657–1661 (2025).
6. Wang, X., Kuipers, J., Hugi, F., Platt, R. J. & Beerenwinkel, N. Reliable single-cell perturbations explain and improve model performance. Preprint at *bioRxiv* https://doi.org/10.64898/2026.08.11.744177 (2026).
7. Arc Institute. cell-eval, version 0.8.2. https://github.com/ArcInstitute/cell-eval (2026).
8. Spearman, C. Correlation calculated from faulty data. *Br. J. Psychol.* **3**, 271–295 (1910).
9. Brown, W. Some experimental results in the correlation of mental abilities. *Br. J. Psychol.* **3**, 296–322 (1910).
10. Molina, A. & Zhang, X. Perturbation response decomposition enables biologically aligned generalization to unseen perturbations and cellular contexts. Preprint at *bioRxiv* https://doi.org/10.64898/2026.07.24.740459 (2026).
11. Schmierer, B. et al. CRISPR/Cas9 screening using unique molecular identifiers. *Mol. Syst. Biol.* **13**, 945 (2017).
12. Sivakumar, S. et al. Benchmarking and optimizing Perturb-seq in differentiating human pluripotent stem cells. *Stem Cell Rep.* **20**, 102713 (2025).
13. Regalado, S. G. et al. Barcoded monoclonal embryoids are a potential solution to confounding bottlenecks in mosaic organoid screens. Preprint at *bioRxiv* https://doi.org/10.1101/2025.05.23.655669 (2025).
14. Squair, J. W. et al. Confronting false discoveries in single-cell differential expression. *Nat. Commun.* **12**, 5692 (2021).
15. Zimmerman, K. D., Espeland, M. A. & Langefeld, C. D. A practical solution to pseudoreplication bias in single-cell studies. *Nat. Commun.* **12**, 738 (2021).
16. Raju, P. C. Directional coherence and effect magnitude in single-cell CRISPR perturbation responses. Preprint at *arXiv* https://arxiv.org/abs/2604.16642 (2026).
17. Rosenbluh, J. et al. Complementary information derived from CRISPR Cas9 mediated gene deletion and suppression. *Nat. Commun.* **8**, 15403 (2017).
18. Huang, X. random-halves-not-replicates, version 1.0.0. *Zenodo* https://doi.org/10.5281/zenodo.23476350 (2026).

