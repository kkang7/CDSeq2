# Install dependencies into a separate library; never activate this project globally.
# Rscript reproduction/environment/restore.R /absolute/path/to/benchmark-library
a<-commandArgs(TRUE);stopifnot(length(a)==1L)
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);here<-dirname(normalizePath(script))
lib<-normalizePath(a[1],mustWork=FALSE);dir.create(lib,recursive=TRUE,showWarnings=FALSE)
if(!requireNamespace('renv',quietly=TRUE))install.packages('renv',lib=lib,repos='https://cloud.r-project.org')
.libPaths(c(lib,.libPaths()))
renv::restore(project=here,lockfile=file.path(here,'renv.lock'),library=lib,prompt=FALSE)
message('Dependencies restored. Install the archived CDSeq2 source into this library for exact benchmark replay.')
