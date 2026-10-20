#' Summarize a gene-type prioritization decision
#'
#' Reports coverage and model-domain diagnostics without estimating the
#' probability that an association is correct. No selection or scores are changed.
#' @param selected Output from [select_variation_profiles()].
#' @return One row per cell type, with total, eligible, retained, missing-score,
#'   unsupported and out-of-range counts and retained fractions. Model identity
#'   and context, if recorded, are attached as attributes.
#' @export
summarize_variation_selection <- function(selected) {
  required <- c('gene','cell_type','variation_score','retained')
  if (!is.data.frame(selected) || !all(required %in% names(selected)) ||
      !is.logical(selected$retained) || anyNA(selected$retained) ||
      !is.numeric(selected$variation_score) || anyNA(selected$gene) ||
      anyNA(selected$cell_type) || anyDuplicated(paste(selected$gene,selected$cell_type,sep='\r')))
    stop('Supply unmodified output from select_variation_profiles().',call.=FALSE)
  if (any(selected$retained & !is.finite(selected$variation_score)))
    stop('A retained profile has no finite variation score.',call.=FALSE)
  ans <- lapply(unique(selected$cell_type),function(ct) {
    d <- selected[selected$cell_type==ct,,drop=FALSE]
    eligible <- sum(is.finite(d$variation_score)); kept <- sum(d$retained)
    outside <- if ('outside_training_range' %in% names(d)) d$outside_training_range else rep(NA,nrow(d))
    unsupported <- if ('status' %in% names(d)) sum(d$status=='unsupported_cell_type',na.rm=TRUE) else NA_integer_
    data.frame(cell_type=ct,total=nrow(d),eligible=eligible,retained=kept,
      fraction_total=kept/nrow(d),fraction_eligible=if(eligible)kept/eligible else NA_real_,
      missing_score=nrow(d)-eligible,unsupported=unsupported,
      outside_training_range=sum(outside%in%TRUE),
      retained_outside_training_range=sum(d$retained & outside%in%TRUE),
      unknown_range=sum(is.na(outside)),stringsAsFactors=FALSE)
  })
  out <- if(length(ans))do.call(rbind,ans)else data.frame()
  attr(out,'model_id') <- attr(selected,'model_id')
  attr(out,'context') <- attr(selected,'context')
  out
}
