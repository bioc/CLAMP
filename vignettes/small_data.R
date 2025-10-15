## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, cache=F)
library(PLIER2)

## -----------------------------------------------------------------------------
data("dataWholeBlood")
data("majorCellTypes")

# Scale each gene to mean 0 and variance 1
dataWholeBlood <- tscale(dataWholeBlood)

## -----------------------------------------------------------------------------
# Download pathway and cell marker libraries from Enrichr
gmtList <- list(
  CellMarkers = getGMT(
    "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=CellMarker_2024",
    "CellMarker_2024"
  ),
  KEGG = getGMT(
    "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human",
    "KEGG_2021_Human"
  )
)

# Combine into a single sparse matrix
pathMatCell <- gmtListToSparseMat(gmtList)

# Load additional xCell reference matrix
data("xCell")

# Match pathways to the gene space of whole blood
matchedPathsWB <- getMatchedPathwayMatList(
  pathMatCell, xCell,
  new.genes = rownames(dataWholeBlood),
  min.genes = 2
)

## -----------------------------------------------------------------------------
# Compute truncated SVD for initialization
set.seed(1);dataWholeBlood.svd=rsvd(dataWholeBlood)

## -----------------------------------------------------------------------------
# Fit baseline PLIER
suppressMessages({
  wb.plier.base <- PLIERbase(
    dataWholeBlood,
 #   k = 25,
    svdres = dataWholeBlood.svd,
    trace = FALSE,
    adaptive.p = 0.05
  )
})

## -----------------------------------------------------------------------------
# Fit variance-prior extended PLIER model
suppressMessages({
  wb.plier.full.vp <- PLIERfullVP(
    dataWholeBlood,
    priorMat = matchedPathsWB,
    plier.base.result = wb.plier.base,
    trace = TRUE, use_cpp = T
  )
})

## ----fig.width=6, fig.height=5------------------------------------------------
output <- compareBs(
  wb.plier.base,
  wb.plier.full.vp,
  celltypeTargets,
  "PLIERbase",
  "PLIERfullVP",
  m = "s"
)

output$plot

## -----------------------------------------------------------------------------
knitr::kable(output$df, format = "markdown")

## -----------------------------------------------------------------------------
dataWholeBloodFBM=bigstatsr::as_FBM(dataWholeBlood)
suppressMessages({
  wb.plier.full.vp.fbm <- PLIERfullVP(
    dataWholeBloodFBM,
    priorMat = matchedPathsWB,
    plier.base.result = wb.plier.base,
    trace = TRUE, use_cpp = T
  )
})


## ----fig.width=6, fig.height=5------------------------------------------------
output <- compareBs(
  wb.plier.base,
  wb.plier.full.vp.fbm,
  celltypeTargets,
  "PLIERbase",
  "PLIERfullVP",
  m = "s"
)

output$plot

