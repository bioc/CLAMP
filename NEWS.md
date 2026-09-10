# CLAMP 0.99.8

## Bug fixes
- Fixed whole-blood examples in the introductory vignette and README to
  use the supplied normalized expression values directly, avoiding redundant
  CPM normalization and log2 transformation that caused the pkgdown build
  to fail during `CLAMPfull()`.

# CLAMP 0.99.7

## Improvements
- Added `log2_transform = TRUE` to `preprocessCLAMP()`,
  `preprocessCLAMPFBM()`, and `cleanFBM()`. Set it to `FALSE` to skip
  log2 transformation while retaining missing-value handling and filtering.

## Bug fixes
- Changed `preprocessCLAMP()` to use population variance, matching
  `preprocessCLAMPFBM()` for variance filtering and downstream scaling.

# CLAMP 0.99.6

## Bug fixes
- Corrected the pathway dot-plot heading and Visualization examples in the
  introductory vignette so all six plots are rendered.

# CLAMP 0.99.5

## Improvements
- Added automated BiocCheck and cross-platform R CMD check workflows.
- Added automatic pkgdown deployment after updates to `devel`.

# CLAMP 0.99.4

## Bug fixes
- Made `preprocessCLAMP()` conditionally log2-transform and replace missing
  values before filtering, matching `preprocessCLAMPFBM()`.
- Made `CLAMPplotU()` handle missing FDR values when no pathways pass the
  requested thresholds.

# CLAMP 0.99.0

## New features
- Initial Bioconductor submission.
- Implements prior-informed latent variable decomposition for gene expression.
- Supports file-backed matrices (`bigstatsr::FBM`) for large datasets.
- Added CPM, z-score, and filtering functions for preprocessing.
- Projection of new datasets into pre-trained models.
- Vignettes with example workflows and detailed documentation.

## Improvements
- Optimized in-place operations for large matrices.
- Added cross-validation utilities and adaptive sparsity support.
- Roxygen2 documentation for all exported functions.

## Bug fixes
- Fixed edge case in CPM normalization for zero-count columns.
- Corrected NA handling in preprocessing.
