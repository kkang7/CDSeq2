# Upgrade legacy Seurat slots without changing the observed RNA counts.
read_rosmap_source <- function(path) {
  object<-readRDS(path)
  if(!inherits(object,'Seurat'))stop('Expected a Seurat object: ',path)
  raw <- methods::slot(methods::slot(object,'assays')[['RNA']],'counts')
  before <- digest::digest(list(raw=raw,metadata=methods::slot(object,'meta.data')),algo='sha256')
  object<-SeuratObject::UpdateSeuratObject(object)
  after_raw<-methods::slot(methods::slot(object,'assays')[['RNA']],'counts')
  after<-digest::digest(list(raw=after_raw,metadata=methods::slot(object,'meta.data')),algo='sha256')
  if(!identical(before,after))stop('Seurat object update changed raw counts or metadata: ',path)
  message('Verified unchanged counts and cell metadata after object update: ',basename(path))
  object
}
