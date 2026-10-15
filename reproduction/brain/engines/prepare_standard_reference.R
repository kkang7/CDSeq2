# Pooled-reference preparation extracted unchanged from the audited input checks.
# No BayesPrism inference, output selection or ranking analysis is defined here.
.bpc_check_source <- function(path,spec) {
  recorded <- spec$input_md5
  j <- which(basename(names(recorded))==basename(path))
  .eqs_assert(length(j)==1L && identical(.eqs_md5(path),unname(recorded[j])),
    'Input no longer matches the saved eQTL experiment: ',path)
}
.bpc_read_inputs <- function(run_dir,sim_dir,meta,d) {
  .eqs_assert(identical(meta$spec$version,'three_type_eqtl_1.0.0'), 'Expected the completed three-type v1 experiment.')
  paths <- c(input=file.path(run_dir,'inputs','direct_inputs.rds'),reference=file.path(run_dir,'inputs','reference.rds'))
  for(p in paths) .bpc_check_source(p,meta$spec)
  inp <- .eqs_load(paths['input']); ref <- .eqs_load(paths['reference'])
  a <- inp$validation$truth; axes <- dimnames(a)
  .eqs_assert(length(dim(a))==3L&&length(axes)==3L, 'Missing validation component axes.')
  for(v in axes) .eqs_assert(length(v)>0L&&!anyNA(v)&&!anyDuplicated(v)&&all(nzchar(v)), 'Invalid component identifiers.')
  genes <- axes[[1L]]; types <- axes[[2L]]; samples <- axes[[3L]]
  for(nm in c('training','validation')) {
    ids <- inp[[nm]]$donors
    .eqs_assert(length(ids)==dim(inp[[nm]]$truth)[3L]&&!anyNA(ids)&&!anyDuplicated(ids)&&all(nzchar(ids)), 'Invalid biological donor IDs.')
  }
  refs <- inp$reference_donors
  .eqs_assert(length(refs)>0L&&!anyNA(refs)&&!anyDuplicated(refs), 'Missing reference donor IDs.')
  .eqs_assert(!any(inp$training$donors%in%inp$validation$donors)&&
    !any(c(inp$training$donors,inp$validation$donors)%in%refs), 'Biological donor overlap.')
  rc <- ref$counts; ann <- as.character(ref$annotation)
  .eqs_assert(is.matrix(rc)&&is.numeric(rc)&&setequal(rownames(rc),genes)&&length(ann)==ncol(rc)&&
    !anyNA(ann)&&setequal(ann,types)&&all(is.finite(rc))&&all(rc>=0), 'Invalid original single-nucleus reference.')
  # One state per original type: summed raw counts are exactly the input that
  # BayesPrism's count.matrix route collapses to before reference normalization.
  pooled <- vapply(types,function(ct)rowSums(rc[genes,ann==ct,drop=FALSE]),numeric(length(genes)))
  dimnames(pooled) <- list(genes,types)
  mu <- sweep(pooled,2,colSums(pooled),'/')
  .eqs_assert(all(rowSums(pooled)>0)&&all(colSums(pooled)>0)&&
    setequal(rownames(inp$reference_mean),genes)&&setequal(colnames(inp$reference_mean),types)&&
    max(abs(mu-inp$reference_mean[genes,types,drop=FALSE]))<1e-12, 'Pooled reference differs from original support reference.')
  .eqs_assert(setequal(unique(d$cell_type),types), 'Saved tests and reference cell types differ.')
  list(reference=t(pooled),mu=mu,genes=genes,types=types,samples=samples,
    donor_checks=data.frame(check=c('reference','training','validation','reference_training_overlap',
      'reference_validation_overlap','training_validation_overlap'),
      n=c(length(refs),length(inp$training$donors),length(inp$validation$donors),0L,0L,0L)),
    reference_cells=data.frame(cell_type=types,n_cells=as.integer(table(factor(ann,levels=types)))),
    checksums=tools::md5sum(paths))
}
