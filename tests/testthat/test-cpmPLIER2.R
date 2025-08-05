library(testthat)

test_that("cpmPLIER2 computes counts per million correctly", {
  cnt <- matrix(c(1,2,3,4), nrow = 2)
  # col sums = (1+3)=4, (2+4)=6
  got <- cpmPLIER2(cnt)
  expect_equal(got[,1], cnt[,1]/4 * 1e6)
  expect_equal(got[,2], cnt[,2]/6 * 1e6)
})
