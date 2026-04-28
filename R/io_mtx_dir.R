#' Build pseudobulk count matrices from a directory of per-sample 10x mtx
#'
#' Aggregates a directory of per-sample 10x triplet matrices into
#' per-(sample x cell-type) pseudobulk count matrices. Expected layout:
#'
#' ```
#' <root>/
#'   <sample1>/
#'     mtx_triplet/
#'       matrix.mtx.gz
#'       barcodes.tsv.gz
#'       features.tsv.gz
#'     cell_metadata.tsv     # has cell_barcode + cell-type column(s)
#'   <sample2>/
#'     ...
#' ```
#'
#' One row per cell barcode in `cell_metadata.tsv`; the cell-type column
#' is named by `cell_type_col`. The cell barcode column must be named
#' `cell_barcode` (matches the bundled 10x barcode order).
#'
#' @param root Path to the directory of per-sample subdirectories.
#' @param cell_type_col Name of the column in each sample's
#'   `cell_metadata.tsv` giving the cell-type label.
#' @param samples Optional character vector of sample directory names
#'   (relative to `root`) to include. Defaults to all subdirectories of
#'   `root`.
#' @param mtx_subdir Subdirectory under each sample containing the
#'   `matrix.mtx.gz` triplet. Defaults to `"mtx_triplet"`.
#' @param cell_metadata_filename Name of the per-sample cell-metadata
#'   TSV file. Defaults to `"cell_metadata.tsv"`.
#' @param features_col Which column of `features.tsv.gz` to use as the
#'   gene identifier (rownames of the output matrices). Defaults to 1
#'   (typically the gene symbol in 10x output).
#'
#' @return Same shape as `build_pseudobulk_from_seurat()`: a named list
#'   with `matrices`, `n_cells`, `genes`, `samples`.
#'
#' @export
build_pseudobulk_from_mtx_dir <- function(root,
                                          cell_type_col,
                                          samples                = NULL,
                                          mtx_subdir             = "mtx_triplet",
                                          cell_metadata_filename = "cell_metadata.tsv",
                                          features_col           = 1L) {
  if (!dir.exists(root)) stop("root does not exist: ", root)
  if (is.null(samples)) {
    samples <- list.dirs(root, recursive = FALSE, full.names = FALSE)
  }
  samples <- sort(samples)

  per_sample_pb     <- vector("list", length(samples))
  names(per_sample_pb) <- samples
  per_sample_n_cell <- vector("list", length(samples))
  names(per_sample_n_cell) <- samples
  gene_universe <- NULL

  for (s in samples) {
    d <- file.path(root, s, mtx_subdir)
    if (!dir.exists(d)) {
      warning(sprintf("Skipping '%s': %s not found", s, mtx_subdir))
      next
    }
    M  <- Matrix::readMM(gzfile(file.path(d, "matrix.mtx.gz")))
    bc <- readLines(gzfile(file.path(d, "barcodes.tsv.gz")))
    ft <- utils::read.delim(gzfile(file.path(d, "features.tsv.gz")),
                            header = FALSE, stringsAsFactors = FALSE)
    rownames(M) <- ft[[features_col]]
    colnames(M) <- bc

    if (is.null(gene_universe)) {
      gene_universe <- rownames(M)
    } else if (!identical(rownames(M), gene_universe)) {
      stop(sprintf("gene order mismatch in sample '%s' (different reference?)",
                   s))
    }

    cm_path <- file.path(root, s, cell_metadata_filename)
    if (!file.exists(cm_path)) {
      warning(sprintf("Skipping '%s': %s not found", s, cell_metadata_filename))
      next
    }
    cm <- utils::read.delim(cm_path, stringsAsFactors = FALSE)
    if (!"cell_barcode" %in% colnames(cm)) {
      stop(sprintf("'%s' has no cell_barcode column for sample '%s'",
                   cell_metadata_filename, s))
    }
    if (!cell_type_col %in% colnames(cm)) {
      stop(sprintf("Cell-type column '%s' missing in '%s' for sample '%s'",
                   cell_type_col, cell_metadata_filename, s))
    }
    if (!all(bc %in% cm$cell_barcode)) {
      stop(sprintf("Some matrix barcodes are not in %s for sample '%s'",
                   cell_metadata_filename, s))
    }
    cm <- cm[match(bc, cm$cell_barcode), , drop = FALSE]
    M  <- methods::as(M, "CsparseMatrix")

    types_s <- as.character(cm[[cell_type_col]])
    per_sample_pb[[s]]     <- aggregate_pseudobulk(M, types_s)
    per_sample_n_cell[[s]] <- table(types_s)
  }

  keep <- !vapply(per_sample_pb, is.null, logical(1))
  stacked <- stack_per_celltype(
    per_sample_pb[keep],
    all_samples   = samples,
    gene_universe = gene_universe
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
       genes    = gene_universe,
       samples  = samples)
}
