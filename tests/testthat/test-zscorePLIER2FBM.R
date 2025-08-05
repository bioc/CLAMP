library(testthat)
library(bigstatsr)

test_that("zscorePLIER2FBM applies Z‐score in‐place on an FBM", {
  mat <- matrix(c(2,4,6,8), nrow = 2, byrow = TRUE)
  fbm <- FBM(nrow(mat), ncol(mat), init = mat)
  stats <- list(
    row_means     = rowMeans(mat),
    row_variances = apply(mat, 1, var)
  )

  expect_message(
    zscorePLIER2FBM(fbm, stats, chunk_size = 1),
    "Applying Z-score transformation"
  )
  out <- fbm[,]
  expected <- sweep(sweep(mat, 1, stats$row_means, "-"),
                    1, sqrt(stats$row_variances), "/")
  expect_equal(out, expected)
})
