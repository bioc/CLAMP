#' @keywords internal
#' @noRd
#'
#' @importFrom stats cor sd var t.test smooth quantile median p.adjust coef
#' @importFrom utils download.file data flush.console globalVariables
#' @importFrom Matrix Matrix sparseMatrix colSums colMeans t
#' @importFrom ggplot2 ggplot aes annotate geom_point geom_abline labs theme_minimal
#' @importFrom ggrepel geom_text_repel
#' @importFrom rlang .data
#' @importFrom bigstatsr big_cprodMat big_prodMat big_apply rows_along FBM
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom rsvd rsvd
#' @importFrom irlba irlba
#'
NULL

# silence NSE notes from ggplot/data.frame column names
utils::globalVariables(c("Cor1", "Cor2", "Label"))
