# From repository root: Rscript tools/install_dependencies.R [benchmark]
options(repos=c(CRAN='https://cloud.r-project.org'))
cran <- c('Rcpp','RcppArmadillo','foreach','doParallel','Matrix','mgcv','remotes','BiocManager')
if('benchmark'%in%commandArgs(TRUE))cran<-c(cran,'snowfall','NMF','gplots')
missing <- cran[!vapply(cran,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))install.packages(missing)
bioc <- 'limma'
if('benchmark'%in%commandArgs(TRUE))bioc<-c(bioc,'scran','BiocParallel','Rsamtools')
missing <- bioc[!vapply(bioc,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))BiocManager::install(missing,ask=FALSE,update=FALSE)
if('benchmark'%in%commandArgs(TRUE)&&!requireNamespace('BayesPrism',quietly=TRUE))
  remotes::install_github('Danko-Lab/BayesPrism/BayesPrism@19052052a6f30833b27e2148294459b0ba2f923e',upgrade='never')
message('Dependencies installed. Run R CMD INSTALL . to install this CDSeq2 source.')
message('For exact paper replay, compare versions with the benchmark dependency_versions.csv.')
