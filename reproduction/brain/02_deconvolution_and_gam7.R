script <- sub("^--file=", "", grep("^--file=",commandArgs(),value=TRUE)[1])
source(file.path(dirname(normalizePath(script)),"bootstrap.R"))
direct_action("run")
direct_action("analyze")
