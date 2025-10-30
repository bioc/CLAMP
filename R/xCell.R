#' xCell cell-signature matrix
#'
#' A numeric matrix of cell-type signatures derived from xCell,
#' used to represent the transcriptional profiles of distinct immune
#' and stromal cell populations.
#'
#' @format A numeric matrix with G genes (rows) and C cell types (columns).
#'   Row names are gene symbols; column names are xCell cell-type labels.
#' @usage data(xCell)
#' @keywords datasets
#' @docType data
#' @name xCell
#' @aliases xCell
#' @examples
#' data(xCell)
#' @return A matrix of xCell-derived reference expression signatures.
"xCell"
