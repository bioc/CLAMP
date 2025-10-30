test_that("CLAMPfull and returns B, Z, U", {
  set.seed(1)

  data("dataWholeBlood", package = "CLAMP")
  data("xCell",          package = "CLAMP")

  matchedPaths <- getMatchedPathwayMatList(
    xCell,
    new.genes = rownames(dataWholeBlood),
    min.genes = 3
  )
  
  base <- CLAMPbase(
    Y = dataWholeBlood,
    trace = FALSE,
    adaptive.p = 0.05
  )
  
  expect_type(base, "list")
  expect_true(all(c("B", "Z") %in% names(base)))
})