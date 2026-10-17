# Compare numerical inputs, not RDS compression, paths, or serialization metadata.
rosmap_signatures <- function(b, sce, v, req) {
  h <- function(x) digest::digest(x, algo='sha256', serialize=TRUE, serializeVersion=2)
  arr <- function(x) list(dim=dim(x), axes=dimnames(x), values=as.numeric(x))
  rc <- methods::as(SummarizedExperiment::assay(sce,'counts'),'dgCMatrix')
  md <- as.data.frame(SummarizedExperiment::colData(sce))
  c(training=h(list(genes=as.character(b$gene_panel), types=as.character(b$major_celltypes),
      bulk=arr(b$pseudobulk_counts), truth=lapply(b$truth_expression_counts,arr),
      cells=arr(b$truth_cell_counts), samples=as.character(b$target_metadata$sample),
      donors=as.character(b$target_metadata$donor), refs=as.character(b$reference_subjects))),
    reference=h(list(dim=dim(rc),axes=dimnames(rc),i=rc@i,p=rc@p,x=as.numeric(rc@x),
      subject=as.character(md$subject),type=as.character(md$cellType))),
    scaffold=h(list(counts=arr(v$counts),samples=as.character(v$metadata$sample),
      donors=as.character(v$metadata$donor),cells=arr(v$n_cells),
      pooled=arr(v$pooled_reference),refs=as.character(v$reference_donors))),
    request=h(list(mean=arr(req$reference_mean),refs=as.character(req$reference_donors),fit=req$fit)))
}
rosmap_verify_objects <- function(b,sce,v,req,expected_file) {
  expected <- read.csv(expected_file,stringsAsFactors=FALSE)
  observed <- rosmap_signatures(b,sce,v,req)
  bad <- names(observed)[observed != expected$sha256[match(names(observed),expected$input)]]
  if(length(bad) || anyNA(expected$sha256[match(names(observed),expected$input)]))
    stop('Reconstruction differs from manuscript inputs: ',paste(bad,collapse=', '),
      '. Do not use it for exact paper replay. Inspect clustering versions and source files.')
  invisible(observed)
}
