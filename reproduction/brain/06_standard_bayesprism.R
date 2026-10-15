# Current paper comparison: native preprocessing, NO BayesPrism output ranking.
script <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
source(file.path(dirname(normalizePath(script)),'bootstrap.R'))
initialize();engine('eqtl_helpers.R');engine('prepare_standard_reference.R')
stopifnot(file.exists(file.path(eqtl,'COMPLETED')))
d <- read.csv(file.path(eqtl,'all_test_details.csv'),stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(d[c('scenario','engine','gene','cell_type')]),
  all(c('CDSeq2','observed')%in%d$engine),!any(d$engine=='BayesPrism'))
meta <- .eqs_load(file.path(eqtl,'experiment.rds'))
prep <- .bpc_read_inputs(direct,eqtl,meta,d)
dir.create(bp,recursive=TRUE,showWarnings=FALSE)
.eqs_save(prep,file.path(bp,'prepared_reference.rds'))
.eqs_csv(d,file.path(bp,'comparison_test_details.csv'))
# Recover the original CDSeq2 top-40% membership on the common finite
# three-ranking universe. No BayesPrism ranking is computed or used.
cd <- subset(d,engine=='CDSeq2')
mem <- do.call(rbind,lapply(split(cd,interaction(cd$scenario,cd$cell_type,drop=TRUE)),function(z){
  eligible <- which(is.finite(z$scenario_GAM6_score)&is.finite(z$scenario_contribution_score)&
    is.finite(z$scenario_expected_count_score))
  eligible <- eligible[order(-z$scenario_GAM6_score[eligible],z$gene[eligible],method='radix')]
  z$selected <- FALSE;z$selected[head(eligible,floor(.4*length(eligible)+1e-9))] <- TRUE
  z$selection_policy <- 'own_support_or_GAM6';z$ranking <- 'GAM6';z$retention <- .4
  z[c('scenario','engine','gene','cell_type','selected','selection_policy','ranking','retention')]
}))
dir.create(file.path(bp,'analysis'),showWarnings=FALSE)
write.csv(mem,file.path(bp,'analysis/selection_membership.csv'),row.names=FALSE)
cache <- file.path(output_root,'standard_bayesprism')
code <- normalizePath(file.path(bundle,'../paper/code'))
for(i in 1:9) {
  status <- system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',
    shQuote(file.path(code,'run_standard_bayesprism.R')),shQuote(direct),shQuote(cache),i,workers))
  if(status!=0L)stop('Native BayesPrism failed in scenario ',i)
}
message('Completed nine standard BayesPrism scenarios in ',cache)
