# Single-cell DE methods literature review — justification for adding dreamlet

**Date:** 2026-04-28
**Author:** Castaldi lab
**Scope:** Methods justification for the LTRC snRNA-seq COPD case-vs-control DE
analysis (`results/differential_expression/COPD_status/`). Argues for adding
**dreamlet** as a sensitivity analysis alongside the primary limma/voom
pseudobulk pipeline, and for **not** using cell-level methods that ignore
donor structure.

---

## 1. Problem statement

The current primary DE script (`ltrc_scrna_copd_de_basal_abc_2026.Rmd`) uses
pseudobulk + `limma::voomWithQualityWeights` with a per-cell-type sample
filter `n_cells >= 10` (Zhang et al. 2026 *Nat Genet* convention).

Two concerns motivate the literature review:

1. **Differential cell counts per sample.** In our cohort, COPD cases have
   higher per-sample basal-cell counts than controls (median 13 vs 7) and
   higher ABC counts (median 23 vs 11). The `n_cells >= 10` filter passes
   44/92 cases vs 8/18 controls for Basal (48% vs 44%) and 65/92 vs 10/18
   for ABC (71% vs 56%). Differential filtering by group is a potential
   selection-bias concern — even though the average pass rates are similar,
   cases that pass have systematically more cells than controls that pass.
2. **Many genes (e.g. MUTYH) are zero in 30–55% of pseudobulks**, which the
   `filterByExpr` step removes from the genome-wide test (we force-retain
   for targeted candidates).

The desired property of any alternative method is: **use all donors with
≥1 cell of the cell type without inflating false positives, and
appropriately handle high-zero-fraction genes for UMI count data.**

## 2. The two methodological questions

This problem decomposes into two distinct technical questions:

- **Q1 — How should donor structure be modeled?** Each donor contributes
  many cells; cells from the same donor are not independent observations.
  Treating cells as independent is *pseudoreplication* and inflates false
  discoveries.
- **Q2 — Is structural zero-inflation needed for UMI snRNA-seq counts?**
  If yes, we need a zero-inflated method (ZINB-WaVE, MAST hurdle, iDESC).
  If no, a standard negative-binomial (NB) or limma/voom fit suffices.

Section 3 addresses Q1; Section 4 addresses Q2.

## 3. Cell-level methods that ignore donor structure are biased

The 2021–2025 benchmark literature is unanimous: methods that treat cells
as independent observations produce inflated type-I error in case-vs-control
designs with multiple donors per group.

### 3.1 Squair et al. 2021 *Nat Commun* — foundational evidence

Squair et al. 2021 (*Nat Commun* 12:5692,
[doi:10.1038/s41467-021-25960-2](https://doi.org/10.1038/s41467-021-25960-2))
tested 14 DE methods on real cross-replicate datasets where the null was
known. Their Abstract states:

> "Methods that ignore this inevitable variation are biased and prone to
> false discoveries. Indeed, the most widely used methods can discover
> hundreds of differentially expressed genes in the absence of biological
> differences."

In Results — *"DE analysis of single-cell data must account for biological
replicates"*:

> "Failing to account for biological replicates causes single-cell methods
> to systematically underestimate the variance of gene expression."

In Results — *"Pseudobulk methods outperform generic and specialized
single-cell DE methods"*:

> "All six of the top-performing methods shared a common analytical
> property. These methods aggregated cells within a biological replicate,
> to form so-called 'pseudobulks', before applying a statistical test."

The six top performers included pseudobulk variants of DESeq2, edgeR, and
limma-voom (their Fig. 2; the prose quoted above does not enumerate the
packages but the figure does).

### 3.2 Murphy & Skene 2022 *Nat Commun* — under MCC, pseudobulk wins

Murphy & Skene 2022 (*Nat Commun* 13:7851,
[doi:10.1038/s41467-022-35519-4](https://doi.org/10.1038/s41467-022-35519-4))
specifically rebuts an earlier Zimmerman et al. claim that mixed-effects
models outperform pseudobulk. Under the Matthews correlation coefficient
(a balanced metric robust to class imbalance):

> "we corrected these issues, reran the author's analysis and found that
> pseudobulk methods outperformed mixed models"

> "Our MCC analysis demonstrates that pseudobulk approaches achieve highest
> performance across individuals and cells variations"

Specifically on mixed-effect approaches:

> "Interestingly, we show that the two mixed model approaches ('Two-part
> hurdle: RE' and 'GLMM Tweedie') perform relatively poorly even compared
> to some pseudoreplication approaches"

### 3.3 Junttila et al. 2022 *Brief Bioinform* — pseudobulk > mixed models

Junttila et al. 2022 (*Brief Bioinform* 23(5):bbac286) tested 18 methods
and concluded:

> "Our results suggest that the pseudobulk methods performed generally
> best."

> "Overall, the pseudobulk methods outperformed the mixed models."

The mixed-model methods tested included NEBULA-LN, MAST_RE, and muscat_MM.
Notably:

> "MAST_RE achieved slightly better overall performance compared to
> NEBULA-LN, whereas muscat_MM was inferior"

> "precision and specificity were higher for the pseudobulk methods than
> for the mixed models"

So **muscat's mixed-model mode** (`muscat_MM`) underperforms even other
mixed models — consistent with reports elsewhere — and should not be used.

### 3.4 Implication for our analysis

Methods that operate on individual cells *without* a donor random effect
or pseudobulk aggregation — `Wilcoxon` on cells, vanilla `MAST`,
`ZINB-WaVE`-weighted edgeR, `Seurat::FindMarkers` defaults, `scVI`'s
built-in DE — are inappropriate here. They will inflate false positives.

## 4. Structural zero-inflation is unnecessary for UMI scRNA-seq

The 2019–2022 literature is clear that UMI-based scRNA-seq counts (which
includes our 10x snRNA-seq data) do **not** have the technical excess of
zeros that motivated zero-inflated methods on read-count data.

### 4.1 Townes et al. 2019 *Genome Biol* — direct empirical demonstration

Townes et al. 2019 (*Genome Biol* 20:295,
[doi:10.1186/s13059-019-1861-6](https://doi.org/10.1186/s13059-019-1861-6))
state in the Abstract:

> "Using negative controls, we show UMI counts follow multinomial sampling
> with no zero inflation."

In Background:

> "However, as shown below, the sampling distribution of UMI counts is not
> zero inflated and differs markedly from read counts, so application of
> read count models to UMI counts needs either theoretical or empirical
> justification."

In Results and Discussion:

> "The results suggest that while read counts appear zero-inflated and
> multimodal, UMI counts follow a discrete distribution with no zero
> inflation."

> "In particular, the empirical probability of a gene being zero across
> droplets was well calibrated to the theoretical prediction based on the
> multinomial model. This also demonstrates that UMI counts are not zero
> inflated…"

In Methods:

> "When applied to UMI counts, the multinomial, Dirichlet-multinomial, and
> Poisson (as approximation to multinomial) distributions fit best."

### 4.2 Jiang et al. 2022 *Genome Biol* — ZI modeling does not improve UMI DE

Jiang et al. 2022 (*Genome Biol* 23:31, "Statistics or biology: the
zero-inflation controversy about scRNA-seq data") confirm Townes' result
on multiple UMI datasets:

> "UMI counts are not zero-inflated when compared with the Poisson or NB
> distribution"

> "non-zero-inflated distributions (Poisson and NB) are chosen for almost
> all genes in the 10x Genomics and Drop-seq datasets"

> "For the two UMI-based protocols, Monocle3 and Seurat have comparable
> performance in terms of F₁ scores"

The implication for DE: tools that explicitly model structural ZI on UMI
data confer no measurable benefit.

### 4.3 Svensson 2020 *Nat Biotechnol* — paywalled (PDF needed for direct quote)

Svensson 2020 (*Nat Biotechnol* 38:147–150, "Droplet scRNA-seq is not
zero-inflated") reaches the same conclusion as Townes 2019. The Letter
itself is paywalled and PubMed has no abstract; the public precursor is a
2017 nxn.se blog post by the same author. **If you want a verbatim quote
from the published Letter, please provide the PDF.** The blog precursor
states (for orientation, not for citation):

> "These observed zeros are consistent with count statistics, and droplet
> scRNA-seq protocols are not producing higher numbers of 'dropouts' than
> expected because of technical artifacts."

For the methods document, Townes 2019 and Jiang 2022 are sufficient
citations — both peer-reviewed open-access primary sources making the
same claim.

### 4.4 Implication for our analysis

We should **not** use ZINB-WaVE, iDESC's zero-inflated component, or
MAST's hurdle for our snRNA-seq data. These add modeling complexity for
no benefit. iDESC's actual statistical contribution is its NB mixed
model with subject random effect — which **NEBULA** also provides, faster
and more actively maintained.

## 5. Two viable options that use all donors

Given Sections 3 and 4, the question reduces to: how do we use **all**
donors with ≥1 cell while properly handling the donor structure and
without ZI modeling? Two viable methods:

### 5.1 dreamlet (Hoffman et al., Bioconductor 2024)

dreamlet (Hoffman et al. bioRxiv 2023.03.17.533005, Bioconductor 2024;
[doi:10.1101/2023.03.17.533005](https://doi.org/10.1101/2023.03.17.533005))
builds pseudobulk for every (cell-type × donor) with ≥1 cell and applies
**precision weights** computed from cell counts so that low-cell-count
pseudobulks are appropriately down-weighted. The Methods section
("*Modeling measurement uncertainty in count models using a two-stage
weighting approach*") justifies the weights:

> "gene expression measurements aggregated across a larger number of
> cells more precisely represents expression in the full population of
> cells. Consequently, the precision of expression measurements of a gene
> across samples is directly related to the number of cells sequenced per
> specimen."

dreamlet supports donor-level random effects via the
`variancePartition`/`dream` framework:

> "Use of precision-weighted linear mixed models enables accounting for
> repeated measures study designs, including technical or biological
> replicates, high dimensional batch effects due to sample multiplexing"

> "Our motivation for using weighted linear regression to model
> transformed counts follows that of Law, et al, but becomes more
> pressing with repeated measures and complex study designs."

Note that "no hard cell-count filter is required" is **inferential** from
the precision-weighting design, not stated as a single sentence in the
preprint. The dreamlet vignette provides operational guidance.

### 5.2 NEBULA-HL (He et al. 2021 *Commun Biol*) — fast cell-level NB GLMM

NEBULA (He et al. 2021 *Commun Biol* 4:629,
[doi:10.1038/s42003-021-02146-6](https://doi.org/10.1038/s42003-021-02146-6))
is a true cell-level NB GLMM that scales to large datasets. It uses
NEBULA-LN (Laplace approximation) by default and switches to NEBULA-HL
(higher-order Laplace) for low-count genes:

> "We, therefore, developed an efficient higher-order LA method to
> correct for this skewness for analyzing low-expression genes"

> "We use NEBULA-LN as a first-line solution to estimate σ² and ϕ. If
> n̄ᵢϕ̂ or κ is smaller than a predefined threshold, we then resort to
> NEBULA-HL"

It is dramatically faster than `glmer.nb`:

> "These benchmarks were on average >100-fold and ~50-fold faster than
> glmer.nb with nAGQ = 0, respectively."

> "In our following analysis of ~34,000 excitatory neurons from 48
> subjects in the snRNA-seq data adopted from ref. 5, NEBULA accomplished
> an association analysis of ~16,000 genes for identifying marker genes
> in ~40 min, compared to ~67 hours using glmer.nb"

The benchmark uses `glmer.nb`; we are not citing a verbatim `glmmTMB`
comparison.

## 6. The Lee & Han 2024 equivalence result

Lee & Han 2024 (*Bioinformatics* 40(8):btae498) prove that
**properly-offset pseudobulk and a single-cell NB GLMM produce equivalent
estimates** in case-vs-control designs.

Discussion:

> "offset-pseudobulk has almost the same statistical properties, such as
> point estimates and standard errors, as GLMMs"

Results:

> "The logFC estimates of offset-pseudobulk and NB GLMM were nearly
> identical across all 5 cell types and 100 trials."

Materials and Methods:

> "the solution for regression (8) is analytically equivalent to the
> solution for [Equation (5)] exactly"

Abstract:

> "Offset-pseudobulk is substantially faster (>×10) and numerically more
> stable than GLMMs."

Results:

> "at N = 200 and 80 cells per subject, running DGE for 20 cell types
> (10K genes) takes only 15 minutes in pseudobulk but 3 hours in NB
> GLMM" (≈ 12× speedup)

This implies **dreamlet (precision-weighted pseudobulk + LMM via dream)
and NEBULA-HL (cell-level NB GLMM) should agree at the estimate level**
in case-vs-control designs of our type. Use dreamlet as the primary
sensitivity (faster, more stable, established voom-style precision
weights); NEBULA as a third sensitivity if disagreement appears.

## 7. Heumos et al. 2025 *Brief Bioinform* — atlas-scale recommendation

Heumos et al. 2025 (*Brief Bioinform* 26(4):bbaf397, "Single-cell
differential expression analysis between conditions within nested
settings") benchmarks specifically nested case-vs-control designs at
atlas scale. The Discussion/Conclusion:

> "For increasingly larger datasets, the user should consider DREAM."

Abstract:

> "For atlas-level analysis, permutation-based methods excel in
> performance but show poor runtime, suggesting to use DREAM as a
> compromise between quality and runtime."

The "DREAM" recommendation refers to the engine inside dreamlet
(`variancePartition::dream`).

The paper does **not** state in a single sentence that "methods ignoring
donor structure are dominated by structure-aware methods" — that is
implicit in the benchmark rankings. We should not cite a verbatim claim
to that effect.

## 8. Recommendation

**Primary DE analysis** (existing): pseudobulk + `voomWithQualityWeights`
with `n_cells >= 10` filter (Zhang 2026 paper convention), in
`ltrc_scrna_copd_de_basal_abc_2026.Rmd`.

**Sensitivity analysis (to add):** dreamlet — pseudobulk for every
(cell-type × donor) with ≥1 cell, precision weights from cell counts,
`dream`-style donor random effect optional. Same model:
`~ copd_analysis_group + sex + age + smoking_status`. Compare:

- Top hits intersection (high-confidence DE).
- Genes recovered only by dreamlet (samples that the `n_cells >= 10`
  filter excluded had real signal).
- FDR-rank concordance (Spearman of `t`-statistics).

**Optional third sensitivity:** NEBULA-HL on individual cells (no
aggregation). Only run if the dreamlet vs limma/voom comparison shows
substantive disagreement. Lee & Han 2024 predicts they will agree.

## 9. Methods to avoid (for our design)

- `Wilcoxon` on individual cells
- vanilla `MAST` without random effects
- `ZINB-WaVE`-weighted edgeR (designed for read-count data; ZI
  unnecessary for UMI per Townes 2019, Jiang 2022)
- `scVI` built-in DE (treats cells as independent)
- `Seurat::FindMarkers` defaults (treats cells as independent)
- `muscat`'s mixed-model mode (`muscat_MM`) — underperforms even other
  mixed models per Junttila 2022
- iDESC's zero-inflated formulation — for UMI data the ZI component is
  unnecessary; the donor-RE component is duplicated by NEBULA

## 10. References

1. Squair JW, Gautier M, Kathe C, et al. Confronting false discoveries
   in single-cell differential expression. *Nat Commun* 12, 5692
   (2021). [doi:10.1038/s41467-021-25960-2](https://doi.org/10.1038/s41467-021-25960-2).
2. Murphy AE, Skene NG. A balanced measure shows superior performance
   of pseudobulk methods in single-cell RNA-sequencing analysis.
   *Nat Commun* 13, 7851 (2022).
   [doi:10.1038/s41467-022-35519-4](https://doi.org/10.1038/s41467-022-35519-4).
3. Junttila S, Smolander J, Elo LL. Benchmarking methods for detecting
   differential states between conditions from multi-subject single-cell
   RNA-seq data. *Brief Bioinform* 23(5):bbac286 (2022).
   [doi:10.1093/bib/bbac286](https://doi.org/10.1093/bib/bbac286).
4. Townes FW, Hicks SC, Aryee MJ, Irizarry RA. Feature selection and
   dimension reduction for single-cell RNA-Seq based on a multinomial
   model. *Genome Biol* 20, 295 (2019).
   [doi:10.1186/s13059-019-1861-6](https://doi.org/10.1186/s13059-019-1861-6).
5. Svensson V. Droplet scRNA-seq is not zero-inflated. *Nat Biotechnol*
   38, 147–150 (2020).
   [doi:10.1038/s41587-019-0379-5](https://doi.org/10.1038/s41587-019-0379-5).
   *Paywalled — PDF needed for verbatim quote.*
6. Jiang R, Sun T, Song D, Li JJ. Statistics or biology: the
   zero-inflation controversy about scRNA-seq data. *Genome Biol* 23, 31
   (2022).
   [doi:10.1186/s13059-022-02601-5](https://doi.org/10.1186/s13059-022-02601-5).
7. Hoffman GE, Lee D, Bendl J, et al. Efficient differential expression
   analysis of large-scale single cell transcriptomics data using
   dreamlet. *bioRxiv* 2023.03.17.533005 (2023; published Bioconductor
   2024).
   [doi:10.1101/2023.03.17.533005](https://doi.org/10.1101/2023.03.17.533005).
8. Hoffman GE, Roussos P. Dream: powerful differential expression
   analysis for repeated measures designs. *Bioinformatics* 37(2):192–201
   (2021).
   [doi:10.1093/bioinformatics/btaa687](https://doi.org/10.1093/bioinformatics/btaa687).
9. He L, Davila-Velderrain J, Sumida TS, et al. NEBULA is a fast
   negative binomial mixed model for differential or co-expression
   analysis of large-scale multi-subject single-cell data.
   *Commun Biol* 4, 629 (2021).
   [doi:10.1038/s42003-021-02146-6](https://doi.org/10.1038/s42003-021-02146-6).
10. Lee H, Han B. Pseudobulk with proper offsets has the same statistical
    properties as generalized linear mixed models in single-cell
    case-control studies. *Bioinformatics* 40(8):btae498 (2024).
    [doi:10.1093/bioinformatics/btae498](https://doi.org/10.1093/bioinformatics/btae498).
11. Heumos L, et al. Single-cell differential expression analysis between
    conditions within nested settings. *Brief Bioinform* 26(4):bbaf397
    (2025).
    [doi:10.1093/bib/bbaf397](https://doi.org/10.1093/bib/bbaf397).
12. Crowell HL, Soneson C, Germain P-L, et al. muscat detects
    subpopulation-specific state transitions from multi-sample
    multi-condition single-cell transcriptomics data.
    *Nat Commun* 11, 6077 (2020).
    [doi:10.1038/s41467-020-19894-4](https://doi.org/10.1038/s41467-020-19894-4).
13. Liu Y, Zhao J, Adams TS, et al. iDESC: identifying differential
    expression in single-cell RNA sequencing data with multiple
    subjects. *BMC Bioinformatics* 24, 318 (2023).
    [doi:10.1186/s12859-023-05432-8](https://doi.org/10.1186/s12859-023-05432-8).
14. Heumos L, Schaar AC, Lance C, et al. Best practices for single-cell
    analysis across modalities. *Nat Rev Genet* 24, 550–572 (2023).
    [doi:10.1038/s41576-023-00586-w](https://doi.org/10.1038/s41576-023-00586-w).

## 11. Caveats in this review

The literature search (sub-agent on 2026-04-28) flagged the following
limitations that affect specific claims:

- **Squair 2021 Claim C** ("DESeq2, edgeR, limma-voom were among the top
  performers") — the explicit package list appears in Fig. 2 / tables;
  the prose quote above states "all six top performers were pseudobulk"
  but does not enumerate the three packages. Cite both prose and figure.
- **Hoffman et al. 2023 — "no hard cell-count filter required"** is
  inferential from the precision-weighting design, not stated verbatim
  in the preprint. The dreamlet vignette and `aggregateToPseudoBulk`
  documentation make this operational.
- **NEBULA speed** — the >100× / 50× speedup vs. `glmer.nb` is
  verbatim-quotable; a direct comparison vs. `glmmTMB` is in the
  supplement, not the abstract.
- **Heumos 2025 Claim B** ("methods ignoring donor structure are
  dominated by structure-aware methods") — implicit in the benchmark
  ranking, not stated verbatim. Avoid this phrasing or weaken to
  "consistent with."
- **Svensson 2020** — paywalled. PDF needed for verbatim quote. Townes
  2019 and Jiang 2022 are sufficient peer-reviewed substitutes for the
  same scientific claim.
