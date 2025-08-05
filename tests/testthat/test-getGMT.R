library(testthat)

test_that("getGMT reads a GMT file into a named list", {
  tmp <- tempfile(fileext = ".gmt")
  writeLines(c(
    "Set1\tfoo\tA\tB\tA",
    "Set2\tbar\tC\tD"
  ), tmp)

  gmt <- getGMT(tmp, name = "foo", cache_dir = tempdir(), redownload = TRUE)
  expect_true(is.list(gmt))
  expect_equal(sort(names(gmt)), c("Set1","Set2"))
  expect_equal(gmt$Set1, sort(c("A","B")))
  expect_equal(gmt$Set2, sort(c("C","D")))
})
