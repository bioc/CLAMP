test_that("cpmPLIER2 computes counts per million correctly", {
    cnt <- matrix(c(1, 2, 3, 4), nrow = 2)
    cs <- colSums(cnt)
    got <- cpmPLIER2(cnt)
    expect_equal(got[, 1], cnt[, 1] / cs[1] * 1e6)
    expect_equal(got[, 2], cnt[, 2] / cs[2] * 1e6)
})
