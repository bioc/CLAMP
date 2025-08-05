library(testthat)
library(Matrix)

test_that("gmtListToSparseMat builds a sparse matrix with correct dims and values", {
  nested <- list(
    lib1 = list(A = letters[1:3], B = letters[2:4]),
    lib2 = list(C = letters[5:6])
  )
  mat <- gmtListToSparseMat(nested)

  expect_s4_class(mat, "dgCMatrix")
  expect_equal(ncol(mat), 3)           # three sets: A, B, C
  expect_equal(sort(rownames(mat)), sort(unique(unlist(unlist(nested)))))

  # membership checks
  expect_equal(mat["a", "A"], 1)
  expect_equal(mat["b", "B"], 1)
  expect_equal(mat["e", "C"], 1)
})
