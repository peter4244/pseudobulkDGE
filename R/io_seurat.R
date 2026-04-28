#' Build pseudobulk count matrices from a Seurat object
#'
#' Aggregates a Seurat object (cell-level) into per-(sample x cell-type)
#' pseudobulk count matrices, one per cell type. The result mirrors what
#' `muscat::aggregateToPseudoBulk()` produces for cell-level
#' SingleCellExperiment input — a list of sparse (genes x samples)
#' matrices keyed by cell type, plus a long-format n_cells table.
#'
#' @param seurat A Seurat object. `seurat[[assay]]@counts` is used as the
#'   source counts matrix.
#' @param sample_id_col Name of the column in `seurat@meta.data` giving
#'   the per-cell sample identifier (the unit of pseudobulk aggregation
#'   on the column dimension). Each unique value becomes one column in
#'   every output matrix.
#' @param cell_type_col Name of the column in `seurat@meta.data` giving
#'   the per-cell cell-type label (the unit of pseudobulk aggregation on
#'   the cell-type dimension). Each unique value becomes one matrix.
#' @param assay Which Seurat assay to use. Defaults to `"RNA"`.
#' @param sample_order Optional character vector specifying the desired
#'   column order in the output matrices. Defaults to
#'   `sort(unique(meta[[sample_id_col]]))`.
#'
#' @return A named list:
#' \itemize{
#'   \item `matrices`: named list of (genes x samples) sparse matrices,
#'     one per cell type.
#'   \item `n_cells`: a long tibble (`celltype`, `sample_name`, `n_cells`)
#'     giving the per-(sample x cell-type) cell count.
#'   \item `genes`: character vector of gene rownames (identical for every
#'     matrix).
#'   \item `samples`: character vector of sample column names (identical
#'     for every matrix).
#' }
#'
#' @export
build_pseudobulk_from_seurat <- function(seurat,
                                         sample_id_col,
                                         cell_type_col,
                                         assay = "RNA",
                                         sample_order = NULL) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package required; install via install.packages('Seurat').")
  }
  meta <- seurat@meta.data
  for (col in c(sample_id_col, cell_type_col)) {
    if (!col %in% colnames(meta)) {
      stop(sprintf("Column '%s' not found in seurat@meta.data", col))
    }
  }
  counts <- Seurat::GetAssayData(seurat, assay = assay, layer = "counts")
  if (ncol(counts) != nrow(meta)) {
    stop("counts has ", ncol(counts), " cells but meta.data has ",
         nrow(meta), " rows")
  }

  samples <- if (is.null(sample_order)) {
    sort(unique(as.character(meta[[sample_id_col]])))
  } else {
    sample_order
  }

  # Aggregate per-sample, then stack per-celltype.
  per_sample_pb     <- vector("list", length(samples))
  names(per_sample_pb) <- samples
  per_sample_n_cell <- vector("list", length(samples))
  names(per_sample_n_cell) <- samples

  for (s in samples) {
    cells_s <- which(as.character(meta[[sample_id_col]]) == s)
    if (!length(cells_s)) {
      next
    }
    counts_s <- counts[, cells_s, drop = FALSE]
    types_s  <- as.character(meta[[cell_type_col]][cells_s])
    per_sample_pb[[s]] <- aggregate_pseudobulk(counts_s, types_s)
    per_sample_n_cell[[s]] <- table(types_s)
  }

  stacked <- stack_per_celltype(
    per_sample_pb[!vapply(per_sample_pb, is.null, logical(1))],
    all_samples   = samples,
    gene_universe = rownames(counts)
  )

  ct_levels <- names(stacked$matrices)
  n_cells <- do.call(rbind, lapply(ct_levels, function(ct) {
    nc <- vapply(samples, function(s) {
      tab <- per_sample_n_cell[[s]]
      if (is.null(tab) || !(ct %in% names(tab))) 0L else as.integer(tab[ct])
    }, integer(1))
    tibble::tibble(celltype = ct, sample_name = samples, n_cells = nc)
  }))

  list(matrices = stacked$matrices,
       n_cells  = n_cells,
       genes    = rownames(counts),
       samples  = samples)
}
