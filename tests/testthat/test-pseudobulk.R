test_that("aggregate_pseudobulk sums match raw counts", {
  set.seed(1)
  M <- Matrix::Matrix(matrix(as.numeric(rpois(60, 1)), nrow = 5),
                      sparse = TRUE)
  colnames(M) <- paste0("c", 1:12)
  rownames(M) <- paste0("g", 1:5)
  grp <- factor(rep(c("A", "B", "C"), each = 4))

  pb <- aggregate_pseudobulk(M, grp)
  expect_s4_class(pb, "CsparseMatrix")
  expect_equal(dim(pb), c(5L, 3L))
  expect_equal(colnames(pb), c("A", "B", "C"))

  # Group sums should match independent rowsum-style aggregation
  M_dense <- as.matrix(M)
  for (lvl in c("A", "B", "C")) {
    expect_equal(as.numeric(pb[, lvl]),
                 as.numeric(rowSums(M_dense[, grp == lvl, drop = FALSE])),
                 info = paste("level", lvl))
  }
})

test_that("aggregate_pseudobulk drops NA group entries", {
  M <- Matrix::Matrix(matrix(1, 3, 4), sparse = TRUE)
  colnames(M) <- paste0("c", 1:4)
  grp <- c("A", NA, "B", "A")
  pb <- aggregate_pseudobulk(M, grp)
  # cells 2 dropped; A = c1+c4, B = c3
  expect_equal(as.numeric(pb[, "A"]), c(2, 2, 2))
  expect_equal(as.numeric(pb[, "B"]), c(1, 1, 1))
})

test_that("aggregate_pseudobulk errors on length mismatch", {
  M <- Matrix::Matrix(matrix(1, 3, 4), sparse = TRUE)
  expect_error(aggregate_pseudobulk(M, c("A", "B")), "ncol")
})

test_that("stack_per_celltype fills missing samples with zero columns", {
  set.seed(2)
  M1 <- Matrix::Matrix(matrix(as.numeric(rpois(20, 1)), nrow = 5),
                       sparse = TRUE)
  rownames(M1) <- paste0("g", 1:5); colnames(M1) <- c("A", "B", "C", "D")
  M2 <- Matrix::Matrix(matrix(as.numeric(rpois(15, 1)), nrow = 5),
                       sparse = TRUE)
  rownames(M2) <- paste0("g", 1:5); colnames(M2) <- c("A", "B", "C")  # no D
  pb_list <- list(s1 = M1, s2 = M2)
  out <- stack_per_celltype(pb_list)

  expect_setequal(names(out$matrices), c("A", "B", "C", "D"))
  # D matrix: s1 column should have M1[, "D"]; s2 column should be all zeros
  expect_equal(as.numeric(out$matrices[["D"]][, "s1"]),
               as.numeric(M1[, "D"]))
  expect_equal(as.numeric(out$matrices[["D"]][, "s2"]),
               rep(0, 5))
})
