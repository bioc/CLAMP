library(testthat)

test_that("zscorePLIER2 centers and scales each row", {
  Y <- matrix(c(2,4,6,8), nrow = 2, byrow = TRUE,
              dimnames = list(c("g1","g2"), c("s1","s2")))
  stats <- data.frame(
    mean     = rowMeans(Y),
    variance = apply(Y,1,var),
    row.names = rownames(Y),
    stringsAsFactors = FALSE
  )
  Z <- zscorePLIER2(Y, stats)

  # expected: each row becomes (original-mean)/sd ⇒ [-1, 1]
  expect_equal(Z, matrix(c(-1,1,-1,1), nrow = 2, byrow = TRUE))
})
