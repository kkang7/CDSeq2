fit_tutorial <- function(cohort,data,seed) {
  fit <- CDSeq2::cdseq2(bulk_data=cohort$bulk,
    reference_scRNA_count=data$reference,
    reference_scRNA_annotation=data$annotation,
    reference_scRNA_donor=data$reference_donors,
    reference_aggregation='pooled_counts',reference_precision=10000,
    mcmc_iterations=100,save_last_n_mc=40,cpu_number=1,n_chains=1L,
    seed=seed,return_mcmc_arrays=FALSE,return_assignments=FALSE,verbose=FALSE)
  genes <- rownames(cohort$bulk);samples <- colnames(cohort$bulk)
  types <- colnames(data$reference_mean)
  # Paper-compatible posterior mean allocated reads, not estGEP_RB parameters.
  counts <- array(0,c(length(genes),length(types),length(samples)),list(genes,types,samples))
  fractions <- matrix(0,length(samples),length(types),dimnames=list(samples,types))
  for(i in seq_along(samples)) {
    counts[,,i] <- fit$cdseq_result[[i]]$estGEP_count
    fractions[i,] <- as.numeric(fit$cdseq_result[[i]]$estProp)
  }
  stopifnot(max(abs(apply(counts,c(1,3),sum)-cohort$bulk))<1e-6,
    max(abs(rowSums(fractions)-1))<1e-8)
  list(counts=counts,fractions=fractions,fit=fit)
}
