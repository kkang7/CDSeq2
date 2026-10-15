association_tutorial <- function(counts,selection,predictor) {
  samples <- dimnames(counts)[[3]]
  stopifnot(identical(names(predictor),samples))
  expression <- log2(sweep(counts,c(2,3),apply(counts,c(2,3),sum),'/')*1e6+.5)
  rows <- list()
  for(k in seq_len(dim(counts)[2])) for(g in seq_len(dim(counts)[1])) {
    y <- expression[g,k,]
    fit <- summary(lm(y~predictor))$coefficients
    rows[[length(rows)+1L]] <- data.frame(gene=dimnames(counts)[[1]][g],
      cell_type=dimnames(counts)[[2]][k],effect=fit[2,1],p=fit[2,4])
  }
  results <- do.call(rbind,rows)
  # Define the full testing family before optional profile selection.
  results$q <- p.adjust(replace(results$p,!is.finite(results$p),1),method='BH')
  key <- function(x)paste(x$gene,x$cell_type,sep='|')
  ix <- match(key(results),key(selection));stopifnot(!anyNA(ix))
  results$retained <- selection$retained[ix]
  results$called_without_prioritization <- is.finite(results$p)&results$q<=.05
  results$called_with_prioritization <- results$called_without_prioritization & results$retained
  results
}
