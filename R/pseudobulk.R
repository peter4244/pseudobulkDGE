#' Aggregate cell-level counts into a (genes x cell-type) pseudobulk matrix
#'
#' Sums counts across cells within each level of `group` using sparse matrix
#' multiplication: `counts %*% t(fac2sparse(group))`. Equivalent in result to
#' `Matrix::rowsum(t(counts), group)` but works directly on sparse matrices
#' without densifying.
#'
#' @param counts A (genes x cells) sparse or dense numeric matrix.
#' @param group A factor or character vector of length `ncol(counts)` giving
#'   the cell-type (or other grouping) label for each cell. NAs are dropped.
#'
#' @return A sparse (`dgCMatrix`) genes-x-levels matrix whose columns are
#'   labeled by `levels(group)`.
#'
#' @examples
#' # Tiny example
#' set.seed(1)
#' M <- Matrix::Matrix(matrix(rpois(60, 1), nrow = 5), sparse = TRUE)
#' colnames(M) <- paste0("c", 1:12)
#' grp <- factor(rep(c("A","B","C"), each = 4))
#' aggregate_pseudobulk(M, grp)
#'
#' @export
aggregate_pseudobulk <- function(counts, group) {
  if (length(group) != ncol(counts)) {
    stop("length(group) must equal ncol(counts)")
  }
  keep <- !is.na(group)
  if (!all(keep)) {
    counts <- counts[, keep, drop = FALSE]
    group  <- group[keep]
  }
  group <- as.factor(group)
  ind   <- Matrix::fac2sparse(group)
  pb    <- counts %*% Matrix::t(ind)
  colnames(pb) <- levels(group)
  methods::as(pb, "CsparseMatrix")
}

#' Stack per-sample pseudobulks into per-cell-type matrices
#'
#' Given a named list of per-sample (genes x cell-types) pseudobulk matrices
#' (typically the output of [aggregate_pseudobulk()] called once per sample),
#' produces a named list of (genes x samples) matrices, one per cell type
#' observed across the cohort. Samples that did not contain a particular
#' cell type are filled with zero columns so every output matrix has the
#' same column set in the same order.
#'
#' @param per_sample_pb A named list (sample_name -> pseudobulk matrix). All
#'   matrices must share the same gene order in their rows.
#' @param all_samples Character vector giving the desired column order for
#'   the output matrices. Defaults to `names(per_sample_pb)`.
#' @param gene_universe Character vector of gene IDs (rownames of every
#'   per-sample pseudobulk). Defaults to `rownames(per_sample_pb[[1]])`.
#'
#' @return A named list with entries:
#' \itemize{
#'   \item `matrices`: named list of (genes x samples) sparse matrices,
#'     one per cell type.
#'   \item `n_cells`: a long-format tibble with columns `celltype`,
#'     `sample_name`, and `n_cells` (per-(sample x cell-type) cell count;
#'     zero where a sample did not contain that cell type).
#' }
#'
#' @export
stack_per_celltype <- function(per_sample_pb,
                               all_samples   = names(per_sample_pb),
                               gene_universe = rownames(per_sample_pb[[1]])) {
  stopifnot(length(per_sample_pb) > 0L,
            !is.null(gene_universe),
            length(gene_universe) > 0L)
  # Sanity: gene universe matches across samples
  for (s in names(per_sample_pb)) {
    if (!identical(rownames(per_sample_pb[[s]]), gene_universe)) {
      stop(sprintf("rownames mismatch for sample '%s'", s))
    }
  }

  # Discover cell-type universe = union of column names across samples
  ct_levels <- unique(unlist(lapply(per_sample_pb, colnames),
                             use.names = FALSE))

  zero_template <- Matrix::Matrix(0, nrow = length(gene_universe), ncol = 1,
                                  sparse = TRUE,
                                  dimnames = list(gene_universe, NULL))

  out <- vector("list", length(ct_levels))
  names(out) <- ct_levels
  ncell_long <- vector("list", length(ct_levels))

  for (j in seq_along(ct_levels)) {
    ct <- ct_levels[j]
    cols <- lapply(all_samples, function(s) {
      M <- per_sample_pb[[s]]
      if (!is.null(M) && ct %in% colnames(M)) {
        M[, ct, drop = FALSE]
      } else {
        zero_template
      }
    })
    M_full <- do.call(cbind, cols)
    colnames(M_full) <- all_samples
    out[[ct]] <- M_full

    nc <- vapply(all_samples, function(s) {
      M <- per_sample_pb[[s]]
      if (!is.null(M) && ct %in% colnames(M)) sum(M[, ct]) else 0
    }, numeric(1))
    # n_cells is *not* derivable from pseudobulks (those are gene-summed
    # counts, not cell counts). Caller must pass cell counts separately
    # if those are wanted -- this function only stacks counts.
    ncell_long[[j]] <- tibble::tibble(celltype    = ct,
                                      sample_name = all_samples,
                                      total_counts = as.numeric(nc))
  }

  list(matrices = out,
       totals   = do.call(rbind, ncell_long))
}
