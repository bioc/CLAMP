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
  
  full <- CLAMPfull(
    Y = dataWholeBlood,
    priorMat = matchedPaths,
    clamp.base.result = base,
    trace = FALSE,
    max.iter = 1,
    max.U.updates = 0
  )

  expect_type(full, "list")
  expect_true(all(c("B", "Z", "U") %in% names(full)))
})