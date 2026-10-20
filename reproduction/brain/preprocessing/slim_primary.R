# Primary preparation only reads raw counts and metadata; discard analysis caches.
read_rosmap_clustered <- function(path) {
  object<-readRDS(path)
  before<-digest::digest(list(counts=methods::slot(object[['RNA']],'counts'),metadata=object[[]]),algo='sha256')
  assay<-object[['RNA']]
  counts<-methods::slot(assay,'counts')
  methods::slot(assay,'data')<-Matrix::sparseMatrix(i=integer(),j=integer(),x=numeric(),dims=dim(counts),dimnames=dimnames(counts))
  methods::slot(assay,'scale.data')<-matrix(numeric(),0,0)
  object[['RNA']]<-assay
  object@reductions<-list();object@graphs<-list();object@neighbors<-list();object@commands<-list()
  after<-digest::digest(list(counts=methods::slot(object[['RNA']],'counts'),metadata=object[[]]),algo='sha256')
  stopifnot(identical(before,after))
  object
}
body(load_rosmap_components) <- parse(text=gsub(
  'readRDS(cfg$input_files[[component]])', 'read_rosmap_clustered(cfg$input_files[[component]])',
  paste(deparse(body(load_rosmap_components)),collapse='\n'),fixed=TRUE))[[1]]
