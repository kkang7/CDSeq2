make_tutorial_benchmark <- function(cohort,inferred,data,label) {
  CDSeq2::prepare_variation_benchmark(bulk=cohort$bulk,
    reference_mean=data$reference_mean,rna_fractions=inferred$fractions,
    estimated_counts=inferred$counts,truth_counts=cohort$truth,
    donor_ids=cohort$donor_ids,reference_donor_ids=unique(data$reference_donors),
    cohort_id=label,oracle='direct_observed_component_counts')
}
prioritize_tutorial <- function(data,training,validation) {
  train <- make_tutorial_benchmark(data$training,training,data,'synthetic_training')
  test <- make_tutorial_benchmark(data$validation,validation,data,'synthetic_validation')
  scorer <- CDSeq2::train_variation_scorer(train,context=list(
    application='Synthetic documentation tutorial',tissue='Synthetic brain-like mixture',
    deconvolution_method='CDSeq2',fit_settings=list(mcmc_iterations=100,
      save_last_n_mc=40,reference_precision=10000,estimator='estGEP_count + estProp')),
    folds=3L,min_profiles=50L,seed=20261015,model='GAM6')
  evaluation <- CDSeq2::validate_variation_scorer(scorer,test)
  features <- CDSeq2::prepare_variation_features(data$validation$bulk,
    data$reference_mean,validation$fractions,scorer=scorer)
  scores <- CDSeq2::score_variation(features,scorer,range_action='flag')
  selected <- CDSeq2::select_variation_profiles(scores,fraction=.2)
  stopifnot(all(table(selected$cell_type[selected$retained])==48),
    !any(data$training$donor_ids %in% data$validation$donor_ids),
    all(evaluation$donor_checks$training_target_overlap==0))
  list(scorer=scorer,evaluation=evaluation,selected=selected)
}
