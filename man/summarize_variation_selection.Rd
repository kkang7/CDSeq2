\name{summarize_variation_selection}
\alias{summarize_variation_selection}
\title{Summarize gene-type prioritization coverage and model-domain diagnostics}
\usage{summarize_variation_selection(selected)}
\arguments{\item{selected}{Output from \code{select_variation_profiles()}.}}
\description{Reports selection coverage, missing scores, unsupported cell types and out-of-range features without changing scores or selection.}
\value{One row per cell type with total, eligible and retained counts, retained fractions and model-domain diagnostics. The score is not a probability of association correctness.}
\seealso{\code{\link{select_variation_profiles}}, \code{\link{score_variation}}}
