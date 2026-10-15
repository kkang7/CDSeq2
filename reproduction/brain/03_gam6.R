script <- sub("^--file=", "", grep("^--file=",commandArgs(),value=TRUE)[1])
source(file.path(dirname(normalizePath(script)),"bootstrap.R"))
initialize(); engine("gam6.R")
# Already in a fresh Rscript process; use the engine directly, retaining .libPaths().
.g6_run(direct,gam6)
