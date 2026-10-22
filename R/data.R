#' Example CDSeq Data
#'
#' A list containing example data used for demonstrating cdseq functions.
#' The object includes matrices and vectors that represent a simulated bulk RNA-seq experiment.
#'
#' @format ## `cdseq_example_data`
#' \describe{
#'   \item{bulk_data}{A numeric matrix of gene expression (genes × samples).}
#'   \item{cstGEP_true}{True cell type-specific gene expression profile for each individual.}
#'   \item{reference_scRNA_count}{Single cell reference.}
#'   \item{reference_scRNA_annotation}{Cell type annotations for reference_scRNA_count.}
#' }
#'
#' @usage data(cdseq_example_data)
#' @docType data
#' @keywords datasets
#' @name cdseq_example_data
#' @examples
#' data(cdseq_example_data)
#' dim(cdseq_example_data$bulk_data)
NULL
