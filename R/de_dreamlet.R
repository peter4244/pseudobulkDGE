#' Per-cell-type pseudobulk DE with dreamlet
#'
#' Builds a pseudobulk `SingleCellExperiment` mimicking
#' `muscat::aggregateToPseudoBulk()`'s structure (one assay per cell
#' type, `int_colData(sce)$n_cells` carrying per-(sample x cell-type)
#' cell counts) directly from already-aggregated pseudobulk matrices,
#' then runs `dreamlet::processAssays()` followed by
#' `dreamlet::dreamlet()` to fit a precision-weighted linear model with
#' voom-style mean-variance modeling.
#'
#' Use `min_cells = 1` (the default) to test every (cell-type x donor)
#' pseudobulk with at least one cell — addressing the selection-bias
#' concern of a hard `n_cells >= 10` filter. Use a relaxed gene filter
#' (`min.count = 1`, `min.samples = 4`, `min.prop = 0.05`) when you
#' want to force-retain a low-count candidate gene.
#'
#' @param pseudobulk Either the named list of (genes x samples) sparse
#'   matrices returned by [build_pseudobulk_from_seurat()] /
#'   [build_pseudobulk_from_mtx_dir()] under `$matrices`, or that whole
#'   list (the function will look for `$matrices`).
#' @param meta_n_cells Long-format tibble/data.frame with columns
#'   `celltype`, `sample_name`, `n_cells`. Typically the `$n_cells` slot
#'   from the IO functions.
#' @param pheno A data.frame keyed by `sample_name` containing the
#'   covariates referenced by `formula`.
#' @param formula Right-hand-side model formula, e.g.
#'   `~ group + sex + age`. Must reference columns of `pheno`.
#' @param coef Coefficient name to extract for the topTable.
#' @param cell_types Character vector of cell types to test. Must be a
#'   subset of `names(pseudobulk$matrices)` (or `names(pseudobulk)`).
#' @param min_cells Minimum cells per (sample x cell-type) for the
#'   sample to enter the dreamlet test. Default 1.
#' @param min_count,min_samples,min_prop Parameters for
#'   `processAssays`'s gene filter. Defaults match dreamlet's own
#'   defaults (5 / 4 / 0.4); relax (e.g. 1 / 4 / 0.05) for targeted
#'   candidate-gene tests.
#' @param normalize_method Normalization method passed to
#'   `processAssays`. Default `"TMM"`.
#' @param adjust_method Multiple-testing correction (default `"BH"`).
#'   Use `"none"` for prespecified single-gene candidate tests.
#' @param confint Whether to include 95% CIs in the topTable.
#'
#' @return A named list with:
#'   \itemize{
#'     \item `fit`: the dreamlet fit object (one
#'       `variancePartition::MArrayLM2` per cell type).
#'     \item `coef_name`: the coefficient that was extracted.
#'     \item `results`: a named list, one per cell type, of tibbles
#'       (gene-level topTables with `assay`, `gene`, `logFC`, `CI.L`,
#'       `CI.R`, `AveExpr`, `t`, `P.Value`, `adj.P.Val`, ...).
#'     \item `n_in_design`: tibble with one row per cell type giving the
#'       number of samples and cases/controls (or rather, the design
#'       matrix dim) used by dreamlet for that cell type.
#'   }
#'
#' @export
run_dreamlet_de <- function(pseudobulk,
                            meta_n_cells,
                            pheno,
                            formula,
                            coef,
                            cell_types,
                            min_cells       = 1L,
                            min_count       = 5L,
                            min_samples     = 4L,
                            min_prop        = 0.4,
                            normalize_method = "TMM",
                            adjust_method   = "BH",
                            confint         = TRUE) {
  if (!requireNamespace("dreamlet", quietly = TRUE) ||
      !requireNamespace("variancePartition", quietly = TRUE)) {
    stop("dreamlet and variancePartition required; install from Bioconductor.")
  }

  pb <- if (is.list(pseudobulk) && !is.null(pseudobulk$matrices))
    pseudobulk$matrices else pseudobulk
  stopifnot(is.list(pb), all(cell_types %in% names(pb)))
  if (!"sample_name" %in% colnames(pheno)) {
    stop("`pheno` must have a `sample_name` column")
  }

  sample_names <- pheno$sample_name

  # Subset matrices to the cell types and samples in `pheno`.
  pb_list <- lapply(setNames(cell_types, cell_types), function(ct) {
    M <- pb[[ct]][, sample_names, drop = FALSE]
    as.matrix(M)
  })

  # Build per-sample n_cells list (named integer vector over cell types).
  meta_n_cells <- as.data.frame(meta_n_cells)
  n_cells_list <- lapply(sample_names, function(s) {
    out <- setNames(integer(length(cell_types)), cell_types)
    rows <- meta_n_cells$sample_name == s &
            meta_n_cells$celltype %in% cell_types
    if (any(rows)) {
      sub <- meta_n_cells[rows, , drop = FALSE]
      out[sub$celltype] <- as.integer(sub$n_cells)
    }
    out
  })
  names(n_cells_list) <- sample_names

  # Build the pseudobulk SCE in muscat-compatible shape.
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays  = pb_list,
    colData = S4Vectors::DataFrame(
      pheno,
      cluster_id = factor(sample_names),
      sample_id  = factor(sample_names),
      row.names  = sample_names
    )
  )
  SingleCellExperiment::int_colData(sce)$n_cells <- n_cells_list
  S4Vectors::metadata(sce)$agg_pars <- list(
    assay = "counts", by = c("cluster_id", "sample_id"),
    fun   = "sum",    scale = FALSE
  )
  S4Vectors::metadata(sce)$aggr_means <- tibble::tibble(
    cluster_id = factor(rep(cell_types, each = length(sample_names))),
    sample_id  = factor(rep(sample_names, times = length(cell_types)))
  )
  S4Vectors::metadata(sce)$experiment_info <- data.frame(
    sample_id = sample_names,
    n_cells   = vapply(sample_names,
                       function(s) sum(n_cells_list[[s]]),
                       integer(1)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  processed <- dreamlet::processAssays(
    sce,
    formula          = formula,
    min.cells        = min_cells,
    min.count        = min_count,
    min.samples      = min_samples,
    min.prop         = min_prop,
    normalize.method = normalize_method,
    quiet            = TRUE
  )
  fit <- dreamlet::dreamlet(processed, formula = formula,
                            min.cells = min_cells, quiet = TRUE)

  all_top <- dreamlet::topTable(
    fit, coef = coef, number = Inf, sort.by = "P",
    confint = confint, adjust.method = adjust_method
  )
  all_top <- as.data.frame(all_top)

  results <- lapply(setNames(cell_types, cell_types), function(ct) {
    df <- all_top[all_top$assay == ct, , drop = FALSE]
    df <- df[order(df$P.Value), , drop = FALSE]
    colnames(df)[colnames(df) == "ID"]    <- "gene"
    colnames(df)[colnames(df) == "assay"] <- "celltype"
    tibble::as_tibble(df)
  })

  n_in_design <- do.call(rbind, lapply(cell_types, function(ct) {
    d <- fit[[ct]]$design
    case_col <- coef
    n_total <- nrow(d)
    n_case  <- if (case_col %in% colnames(d)) sum(d[, case_col] == 1) else NA_integer_
    n_ctrl  <- if (case_col %in% colnames(d)) sum(d[, case_col] == 0) else NA_integer_
    tibble::tibble(celltype = ct,
                   n_samples_in_design = n_total,
                   n_at_case_level     = n_case,
                   n_at_ref_level      = n_ctrl)
  }))

  list(
    fit         = fit,
    coef_name   = coef,
    results     = results,
    n_in_design = n_in_design
  )
}
