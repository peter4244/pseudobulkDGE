#' Per-cell-type pseudobulk DE with limma + voomWithQualityWeights
#'
#' For each cell type in `cell_types`, builds an `edgeR::DGEList` from the
#' (genes x samples) pseudobulk, restricts to samples with at least
#' `min_cells` cells of that cell type, applies TMM normalization, drops
#' low-count genes via `edgeR::filterByExpr` (with optional candidate
#' force-retention), runs `limma::voomWithQualityWeights` (or plain
#' `voom`), fits the model, and returns the empirical-Bayes top-table.
#'
#' Standard contrast extraction; the contrast direction is given by the
#' user-supplied `coef` (e.g. `"groupcase"` for case vs control where the
#' factor reference is `control`).
#'
#' @param pseudobulk Either the named list of (genes x samples) sparse
#'   matrices returned by [build_pseudobulk_from_seurat()] /
#'   [build_pseudobulk_from_mtx_dir()] under `$matrices`, or that whole
#'   list (the function will look for `$matrices`).
#' @param meta_n_cells A long-format tibble or data.frame with columns
#'   `celltype`, `sample_name`, and `n_cells`, identifying the
#'   per-(sample x cell-type) cell count. Typically the `$n_cells` slot
#'   from the IO functions.
#' @param pheno A data.frame keyed by `sample_name` containing the
#'   covariates referenced by `formula`.
#' @param formula Right-hand-side model formula, e.g.
#'   `~ group + sex + age`. Must reference columns of `pheno`.
#' @param coef Coefficient name to extract for the topTable (e.g.
#'   `"groupcase"` if `group` is a 2-level factor with `control` as the
#'   reference). Must match a column of the design matrix produced by
#'   `model.matrix(formula, pheno)`.
#' @param cell_types Character vector of cell types to test. Must be a
#'   subset of `names(pseudobulk$matrices)` (or `names(pseudobulk)` if
#'   already a list of matrices).
#' @param min_cells Minimum per-(sample x cell-type) cell count to
#'   retain a sample for that cell type. Default 10 (Zhang 2026 / Squair
#'   2021 convention).
#' @param force_retain_genes Optional character vector of gene symbols
#'   to force-retain past `filterByExpr`. Useful for prespecified
#'   candidate-gene tests of low-count genes.
#' @param adjust_method Multiple-testing correction method passed to
#'   `topTable` (default `"BH"`). Use `"none"` for prespecified
#'   single-gene candidate tests.
#' @param robust_eb Whether to use `eBayes(robust = TRUE)` (recommended).
#' @param with_quality_weights Whether to use `voomWithQualityWeights`
#'   (recommended for variable per-sample cell counts) or plain `voom`.
#' @param confint Whether to include 95% CIs in the topTable.
#'
#' @return A named list, one element per cell type, each a list with:
#'   \itemize{
#'     \item `samp`: tibble of samples retained for this cell type
#'       (with `n_cells_ct` column added).
#'     \item `n_genes_input`: number of genes before `filterByExpr`.
#'     \item `n_genes_kept`: number of genes after.
#'     \item `dgelist`: filtered `DGEList`.
#'     \item `voom`: voom output (`EList`).
#'     \item `design`: design matrix.
#'     \item `fit`: `lmFit` + `eBayes` fit.
#'     \item `coef_name`: the coefficient that was extracted.
#'     \item `results`: tibble topTable for `coef`.
#'   }
#'
#' @export
run_limma_voom_de <- function(pseudobulk,
                              meta_n_cells,
                              pheno,
                              formula,
                              coef,
                              cell_types,
                              min_cells           = 10L,
                              force_retain_genes  = NULL,
                              adjust_method       = "BH",
                              robust_eb           = TRUE,
                              with_quality_weights = TRUE,
                              confint             = TRUE) {
  pb <- if (is.list(pseudobulk) && !is.null(pseudobulk$matrices))
    pseudobulk$matrices else pseudobulk
  stopifnot(is.list(pb), all(cell_types %in% names(pb)))

  if (!"sample_name" %in% colnames(pheno)) {
    stop("`pheno` must have a `sample_name` column")
  }
  meta_n_cells <- as.data.frame(meta_n_cells)
  for (col in c("celltype", "sample_name", "n_cells")) {
    if (!col %in% colnames(meta_n_cells)) {
      stop(sprintf("`meta_n_cells` must have column '%s'", col))
    }
  }

  results <- vector("list", length(cell_types))
  names(results) <- cell_types

  for (ct in cell_types) {
    n_ct <- meta_n_cells[meta_n_cells$celltype == ct,
                         c("sample_name", "n_cells"),
                         drop = FALSE]
    colnames(n_ct)[2] <- "n_cells_ct"
    samp <- merge(pheno, n_ct, by = "sample_name", all.x = FALSE)
    samp <- samp[samp$n_cells_ct >= min_cells, , drop = FALSE]
    samp <- samp[order(samp$sample_name), , drop = FALSE]

    n_case <- NA_integer_
    if (length(all.vars(formula)) > 0) {
      first_var <- all.vars(formula)[1]
      if (first_var %in% colnames(samp) && is.factor(samp[[first_var]])) {
        # Best-effort guard: ensure both factor levels present at >=2 samples.
        tab <- table(samp[[first_var]])
        if (length(tab) < 2L || any(tab < 2L)) {
          stop(sprintf(
            "Cell type '%s': insufficient samples in factor '%s' after n_cells >= %d filter (counts: %s).",
            ct, first_var, min_cells,
            paste(sprintf("%s=%d", names(tab), as.integer(tab)),
                  collapse = ", ")))
        }
      }
    }

    counts_ct <- pb[[ct]][, samp$sample_name, drop = FALSE]
    d <- edgeR::DGEList(counts = as.matrix(counts_ct), samples = samp)
    d$samples$lib.size <- colSums(d$counts)
    d <- edgeR::calcNormFactors(d, method = "TMM")

    design <- stats::model.matrix(formula, data = samp)
    if (!coef %in% colnames(design)) {
      stop(sprintf("coef '%s' not in design columns: %s",
                   coef, paste(colnames(design), collapse = ", ")))
    }

    keep_default <- edgeR::filterByExpr(d, design = design)
    keep <- if (length(force_retain_genes)) {
      keep_default | (rownames(d) %in% force_retain_genes)
    } else {
      keep_default
    }
    if (!any(keep)) {
      stop(sprintf(
        "filterByExpr removed all genes for cell type '%s'", ct))
    }
    n_genes_input <- nrow(d)
    d <- d[keep, , keep.lib.sizes = FALSE]

    v <- if (with_quality_weights) {
      limma::voomWithQualityWeights(d, design = design, plot = FALSE)
    } else {
      limma::voom(d, design = design, plot = FALSE)
    }
    fit <- limma::lmFit(v, design)
    fit <- limma::eBayes(fit, robust = robust_eb)

    tt <- limma::topTable(fit, coef = coef, number = Inf,
                          sort.by = "P", confint = confint,
                          adjust.method = adjust_method)
    tt <- tibble::rownames_to_column(tt, "gene")
    tt <- tibble::as_tibble(tt)

    results[[ct]] <- list(
      celltype       = ct,
      samp           = tibble::as_tibble(samp),
      n_genes_input  = n_genes_input,
      n_genes_kept   = nrow(d),
      dgelist        = d,
      voom           = v,
      design         = design,
      fit            = fit,
      coef_name      = coef,
      results        = tt
    )
  }

  results
}
