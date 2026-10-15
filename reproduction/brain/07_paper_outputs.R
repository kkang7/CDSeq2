script <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
source(file.path(dirname(normalizePath(script)),'bootstrap.R'))
initialize()
code <- normalizePath(file.path(bundle,'../paper/code'))
out <- file.path(output_root,'paper_outputs')
status <- system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',
  shQuote(file.path(code,'regenerate_figures.R')),shQuote(direct),shQuote(out),
  shQuote(file.path(output_root,'standard_bayesprism'))))
if(status!=0L)stop('Paper output regeneration failed.')
message('Brain benchmark outputs: ',out,'; independent-LD results are a separate reproduction part.')
