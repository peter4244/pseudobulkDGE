# pseudobulkDGE

Reusable pseudobulk differential gene expression for multi-subject
single-cell / single-nucleus RNA-seq.

`pseudobulkDGE` accepts either:

- a **Seurat object** (cell-level), or
- a **directory of per-sample 10x mtx triplets** (with a per-sample
  `cell_metadata.tsv` carrying cell-type labels)

aggregates to **per-(sample × cell-type) pseudobulk** using sparse
matrix multiplication, and runs both:

- **`limma::voomWithQualityWeights`** + `lmFit` + `eBayes(robust=TRUE)`
  with optional candidate-gene force-retention, and
- **dreamlet** (precision-weighted pseudobulk; Hoffman et al.,
  Bioconductor 2024)

with consistent output schemas (gene × celltype tibbles).

The package is engineered for the lab's report convention and is
distilled from the LTRC snRNA-seq COPD analysis pipeline
(`https://changit.bwh.harvard.edu/repjc/ltrc_2026_singlecell_dge`).

## Install

```r
# From Changit (within the network)
remotes::install_git(
  "https://changit.bwh.harvard.edu/repjc/pseudobulkDGE.git"
)

# Or from a local clone during development
remotes::install_local("/path/to/pseudobulkDGE")
```

Bioconductor dependencies (`dreamlet`, `variancePartition`,
`SingleCellExperiment`, `edgeR`, `limma`, `S4Vectors`,
`SummarizedExperiment`) install via `BiocManager`.

## Usage

### 1. Build pseudobulk

```r
library(pseudobulkDGE)

# From a Seurat object:
pb <- build_pseudobulk_from_seurat(
  seurat        = my_seurat,
  sample_id_col = "donor_id",
  cell_type_col = "cell_type"
)

# Or from a directory of per-sample 10x triplets:
pb <- build_pseudobulk_from_mtx_dir(
  root          = "/path/to/per_sample_dirs/",
  cell_type_col = "cell"     # column name in cell_metadata.tsv
)

# pb$matrices is a named list (cell_type -> genes x samples sparse matrix)
# pb$n_cells  is a long tibble of per-(cell_type x sample) counts
# pb$genes    is a vector of gene IDs (rownames)
# pb$samples  is a vector of sample IDs (column order)
```

### 2. Run limma/voom DE

```r
res_lv <- run_limma_voom_de(
  pseudobulk    = pb,
  meta_n_cells  = pb$n_cells,
  pheno         = sample_pheno,            # data.frame with sample_name col
  formula       = ~ group + sex + age,
  coef          = "groupcase",             # contrast direction
  cell_types    = c("Basal", "ABC"),
  min_cells     = 10                       # paper convention
)

res_lv$Basal$results        # per-gene topTable for Basal
res_lv$Basal$n_genes_kept   # genes after filterByExpr
res_lv$Basal$voom           # voom EList for diagnostics
res_lv$Basal$fit            # full lmFit + eBayes object
```

For prespecified candidate-gene tests, force-retain past
`filterByExpr` and report uncorrected P:

```r
res_targeted <- run_limma_voom_de(
  pseudobulk         = pb,
  meta_n_cells       = pb$n_cells,
  pheno              = sample_pheno,
  formula            = ~ group + sex + age,
  coef               = "groupcase",
  cell_types         = c("Basal", "ABC"),
  min_cells          = 10,
  force_retain_genes = c("MUTYH"),
  adjust_method      = "none"
)
```

### 3. Run dreamlet DE (no per-sample cell-count filter)

```r
res_dl <- run_dreamlet_de(
  pseudobulk   = pb,
  meta_n_cells = pb$n_cells,
  pheno        = sample_pheno,
  formula      = ~ group + sex + age,
  coef         = "groupcase",
  cell_types   = c("Basal", "ABC"),
  min_cells    = 1                          # all donors with >=1 cell
)

res_dl$results$Basal     # tibble of gene-level results
res_dl$n_in_design       # samples kept per cell type
```

For targeted candidates with dreamlet, relax the gene filter so the
candidate is retained:

```r
res_dl_targeted <- run_dreamlet_de(
  pseudobulk    = pb, meta_n_cells = pb$n_cells, pheno = sample_pheno,
  formula       = ~ group + sex + age, coef = "groupcase",
  cell_types    = c("Basal", "ABC"),
  min_cells     = 1,
  min_count     = 1, min_samples = 4, min_prop = 0.05,
  adjust_method = "none"
)
```

## Conventions

- **Pheno table** must have a `sample_name` column matching the
  pseudobulk column names (sample IDs).
- **Contrast direction** is given by `coef`. With a 2-level factor
  `factor(group, levels = c("control","case"))` the case-vs-control
  coefficient is `"groupcase"` and **positive logFC = up in cases**.
- For UMI scRNA-seq, structural zero-inflation modeling is **not**
  used (NB / negative binomial fits well per Townes 2019, Svensson
  2020, Cao 2022).

## Methods justification

The dreamlet vs. limma/voom design and the choice of avoiding ZINB
methods are grounded in the 2021–2025 single-cell DE benchmark
literature (Squair 2021 *Nat Commun*, Murphy & Skene 2022, Junttila
2022, Lee & Han 2024, Heumos 2025), with primary-source quotes
catalogued in the LTRC repo at
`results/differential_expression/COPD_status/de_methods_literature_review_2026.md`.

## Status

Pre-1.0. API may change. Use a fixed git ref or commit hash for
reproducibility.

## License

MIT.
