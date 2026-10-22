# Assisted-by: OpenAI Codex. Synthetic checks of inference invariants.
library(CDSeq2)
set.seed(7)
ref <- matrix(rpois(120, 5), 20, 6,
    dimnames = list(paste0('g', seq_len(20)), paste0('c', seq_len(6))))
bulk <- ref[, 1:2] + ref[, 4:5]
colnames(bulk) <- c('s1', 's2')
run_fit <- function() cdseq2(bulk, reference_scRNA_count = ref,
    reference_scRNA_annotation = rep(c('A', 'B'), each = 3),
    reference_aggregation = 'pooled_counts', reference_precision = 100,
    mcmc_iterations = 20, save_last_n_mc = 10, cpu_number = 1,
    seed = 91, verbose = FALSE)
rng_before <- .Random.seed
fit <- run_fit()
stopifnot(identical(rng_before, .Random.seed))
fit_again <- run_fit()
for (i in seq_len(ncol(bulk))) {
    a <- fit$cdseq_result[[i]]
    b <- fit_again$cdseq_result[[i]]
    stopifnot(isTRUE(all.equal(a$estGEP_count, b$estGEP_count)),
        max(abs(rowSums(a$estGEP_count) - bulk[, i])) < 1e-8,
        abs(sum(a$estProp) - 1) < 1e-8,
        all(is.finite(a$estGEP_RB)), all(a$estGEP_count >= 0))
}
# Split storage must still replay the frozen model's original prediction controls.
check_deg_scorer(load_deg_scorer())
