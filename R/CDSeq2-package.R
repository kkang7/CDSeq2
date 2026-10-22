#' CDSeq2: Individual-level cell-type expression deconvolution
#'
#' \code{CDSeq2-R-package} Individual-level complete deconvolution method.
#' @author Kai Kang \email{kangk@@uncw.edu}
#' @section Expression deconvolution:
#' CDSeq2 estimates donor-specific cell-type expression and RNA fractions
#' from bulk RNA sequencing and supports optional gene-type prioritization.
#' @import Rcpp foreach doParallel parallel
#' @importFrom Rcpp evalCpp
#' @importFrom utils head
#' @useDynLib CDSeq2
#' @keywords internal
"_PACKAGE"
NULL
