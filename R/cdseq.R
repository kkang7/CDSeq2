#' Complete deconvolution and cell type-specific differential expression using bulk RNAseq.
#'
#' \code{cdseq} deconvolves individual-level bulk RNASeq and enables cell type-specific differential expression analysis.
#'
#' @param bulk_data RNA-Seq read counts matrix. Columns represent samples and rows represent genes.
#'
#' @param mcmc_iterations number of iterations for the Gibbs sampler; default value is 700.
#'
#' @param cpu_number number of cpu cores that can be used for parallel computing;
#' Default is NULL and cdseq will detect the available number of cores on the device and use number of all cores - 1 for parallel computing.
#'
#' @param reference_scRNA_count single cell reference count matrix.
#'
#' @param reference_scRNA_annotation single cell reference cell type annotations. A vector with length equals to ncol(reference_scRNA).
#'
#' @param reference_proportion cell type proportion reference. It can be any positive scalar or vector of length equals to the number of cell types.
#'
#' @param reference_precision a scalar parameter which intuitively controls the level of confidence in single reference.
#'
#' @param smooth_ref Logical. Whether perform smoothing on single cell reference.
#'
#' @param smoothing_alpha numerical value for smoothing procedure.
#'
#' @param rng random number generator, 1 for Mersenne Twister, 2 for Xoshiro256plus.
#'
#' @param save_last_n_mc number of MCMC samples to save, default is 100.
#'
#' @param adjust_prior logical value: TRUE means adjust prior and FALSE means not adjust prior. Default is FALSE.
#'
#' @param verbose if TRUE, then print progress message to the console. Default is FALSE.
#'
#' @importFrom stats cor
#'
#' @importFrom foreach %dopar%
#'
#' @export
#'
#' @return cdseq returns estimates of both cell-type-specific gene expression profiles and cell-type proportions for each individual.


cdseq <- function(bulk_data,
                  mcmc_iterations = 700,
                  #tuning_iterations = 1,
                  #dilution_factor = 1,
                  cpu_number = 1,
                  reference_scRNA_count = NULL,
                  reference_scRNA_annotation = NULL,
                  reference_proportion = 1,
                  #monitor_prop = TRUE,
                  #low_prop_cutoff = 1e-4,
                  #absent_prop_cutoff = 1e-5,
                  reference_precision = 1e10,
                  #precision_increments = 1,
                  smooth_ref = TRUE,
                  smoothing_alpha = 1,
                  #nsamples = c("tissue"),
                  rng = 1,
                  save_last_n_mc = 100,
                  adjust_prior = FALSE,
                  verbose = FALSE) {

  cat("cdseq - reference-based complete deconvolution and cell type-specific differential expression analysis.\n")
  start_time <- Sys.time()

  ##################################################################
  ##                         check inputs                         ##
  ##################################################################

  #check number of arguments
  #if(nargs()>11){stop("cdseq function requires at most 11 arguments.")}

  #check bulk_data
  if(is.null(bulk_data)){stop("bulk_data is NULL. Please provide unnormalized raw counts data.")}
  if(!is.matrix(bulk_data)){stop("The input bulk_data has to be a matrix.")}
  if(anyDuplicated(colnames(bulk_data)) != 0) {stop("Please make sure bulk_data has unique column names.")}

  #Check mcmc_iterations
  if(is.null(mcmc_iterations)){stop(" mcmc_iterations, the number of MCMC iterations, is missing.")}
  if(mcmc_iterations<1){stop(" mcmc_iterations, the number of MCMC iterations, has to be greater than 1.")}

  #Check dilution_factor
  # if(length(dilution_factor)>1){stop("dilution_factor has to be an positive integer, NOT a vector.")}
  # if(dilution_factor<0){stop("dilution_factor has to be positive.")}

  #Check reference_scRNA
  if(is.null(reference_scRNA_count)){stop("Please provide single cell raw count matrix.")}
  if(is.null(reference_scRNA_annotation)){stop("Please provide single cell annotations.")}
  if(is.null(rownames(reference_scRNA_count))){stop("Please provide gene names in reference single cell data.")}
  if(nrow(reference_scRNA_count)!=nrow(bulk_data)){stop("Dimenson of reference_scRNA and bulk_data do not match.")}
  if(!identical(rownames(reference_scRNA_count), rownames(bulk_data))){stop("Genes of reference_scRNA and bulk_data do not match.")}
  if(save_last_n_mc>mcmc_iterations){stop("save_last_n_mc has to be less than mcmc_iterations.")}
  if(is.null(reference_proportion)){stop("Please provide reference_proportion.")}

  ##################################################################
  ##               save sample names and gene names               ##
  ##################################################################

  # smoothing reference GEP to avoid zero estimates
  if(smooth_ref){
    reference_scRNA_count <- reference_scRNA_count + smoothing_alpha
    reference_scRNA = t(rowsum(t(reference_scRNA_count), group = reference_scRNA_annotation))
    reference_scRNA = t(reference_scRNA)/colSums(reference_scRNA)
  }else{
    reference_scRNA = t(rowsum(t(reference_scRNA_count), group = reference_scRNA_annotation))
    reference_scRNA = t(reference_scRNA)/colSums(reference_scRNA)
  }

  beta = reference_scRNA * reference_precision

  #alpha_initial = 5 # hard coded alpha
  alpha_initial = reference_proportion #* reference_precision

  cell_type_number = length(unique(reference_scRNA_annotation))
  gene_symbols <- row.names(bulk_data)
  sample_names <- colnames(bulk_data)
  cell_types <- NULL

  if(is.null(gene_symbols)){stop("Please provide gene names")}
  if(is.null(sample_names)){stop("Please provide column names for bulk_data.")}

  #Check if gene name is consistent in bulk_data and reference_scRNA_count
  if(!is.null(gene_symbols) && !is.null(row.names(reference_scRNA_count))){
    if( !identical(gene_symbols,colnames(reference_scRNA)) ){
      stop("gene names in bulk_data are different from the gene names in single cell reference.")
    }
  }
  cell_types = row.names(reference_scRNA)
  reference_scRNA = as.matrix(reference_scRNA)


  #################################################################
  ##                        data dilution                        ##
  #################################################################

  # only for speeding up the MCMC when necessary
  #bulk_data <- ceiling(bulk_data/dilution_factor)


  ##################################################################
  ##                individual level deconvolution                ##
  ##################################################################

  nsamples <- colnames(bulk_data)
  samples  <- unique(nsamples)

  # Optionally check that cols are uniquely named:
  if (length(nsamples) != length(samples)) {
    stop("Please make sure bulk_data has unique column names.")
  }

  # Split out each sample’s columns into its own matrix
  individual_sample <- lapply(samples, function(s) {
    bulk_data[, nsamples == s, drop = FALSE]
  })
  names(individual_sample) <- samples

  ###########################
  # parallel computing setup
  ###########################
  if(is.null(cpu_number)){
    doParallel::registerDoParallel(1)
    message("cpu_number is not provided. cdseq uses 1 cpu.\n")
  }else{
    if(!is.numeric(cpu_number) || !(is.atomic(cpu_number) && length(cpu_number)==1) || is.matrix(cpu_number)){stop("cpu_number has to be an integer that is greater than or equal to 1.")}
    #if(cpu_number >  parallel::detectCores()-1 | cpu_number < 1){stop("cpu_number is invalid. To detect how many cores are available, use detectCores().")}
    if(cpu_number%%1!=0){
      cpu_number <- ceiling(cpu_number)
      warning("cpu_number is not an integer. cdseq uses ceiling(cpu_integer).")
    }
    doParallel::registerDoParallel(cpu_number)
  }
  if(verbose){verbose_int=1}else{verbose_int=0}
  celltype_assignment <- NULL
  cellTypeAssignSplit <- NULL

  printout <- 1
  i=1
  # lgpst = rep(0,tuning_iterations)
  # lglike = rep(0,tuning_iterations)
  # lgpst2 = rep(0,tuning_iterations)
  # lglike2 = rep(0,tuning_iterations)
  lgpst = rep(0,1)
  lglike = rep(0,1)
  lgpst2 = rep(0,1)
  lglike2 = rep(0,1)


  #Gibbs sampler and store its running time
  start.time = Sys.time()

  # adjust prior
  if(adjust_prior){
    alpha = matrix(alpha_initial,nrow=cell_type_number,ncol = ncol(bulk_data))
    adjust_prior_result <- gibbsSampler(alpha,beta,bulk_data,cell_type_number,mcmc_iterations, printout, 0, 0 ,verbose_int, rng, 5)
    estGEP_vec<-adjust_prior_result$csGEP_vec
    estSSp_vec<-adjust_prior_result$SSP_vec

    #vector to matrix
    samplesize<-ncol(individual_sample[[i]])
    estGEP_mat<-t(matrix(estGEP_vec,nrow = cell_type_number,ncol = nrow(individual_sample[[i]]))) # gene by cell_type
    estGEP_read2<-t(t(estGEP_mat)/colSums(estGEP_mat))#gene by cell_type
    beta = t(estGEP_read2 * reference_precision )
  }

  #gibbsRunningTime<-system.time(
  allresult <- foreach::foreach(i=1:ncol(bulk_data), .inorder = FALSE,.verbose=FALSE) %dopar% {

    alpha = matrix(alpha_initial,nrow=cell_type_number,ncol = ncol(individual_sample[[i]]))
    # for (itr in 1:tuning_iterations) {
    #   cat("==== Tuning iteration : ", itr," ====\n")
    #   if(itr==1){
    #     result <- gibbsSampler(alpha,beta,individual_sample[[i]],cell_type_number,mcmc_iterations, printout, 0, 0 ,verbose_int, rng, save_last_n_mc)
    #   }
    #   if(itr>1){
    #     # update beta and alpha
    #     reference_precision = reference_precision #* precision_increments
    #     beta = t(estGEP_read2 * reference_precision )
    #     #cat("reference_precision=",reference_precision,"; dim(beta)=",dim(beta),"\n")
    #     alpha = estProp2*reference_precision
    #     result <- gibbsSampler(alpha,beta,individual_sample[[i]],cell_type_number,mcmc_iterations, printout, 0,0,verbose_int, rng,save_last_n_mc)#1, result$cellTypeBag
    #   }
    #
    #   #outputs are two vectors. estGEP_vec is gene x cell type; estSSp_vec is sample x cell type
    #   estGEP_vec<-result$csGEP_vec
    #   estSSp_vec<-result$SSP_vec
    #
    #   #vector to matrix
    #   samplesize<-ncol(individual_sample[[i]])
    #   estGEP_mat<-t(matrix(estGEP_vec,nrow = cell_type_number,ncol = nrow(individual_sample[[i]]))) # gene by cell_type
    #   estSSp_mat<-matrix(estSSp_vec,nrow = cell_type_number,ncol = samplesize) #cell_type by sample_size
    #
    #   # predicted proportions and GEPs
    #   estProp<-t(t(estSSp_mat+alpha)/colSums(estSSp_mat+alpha))  #cell_type by sample_size
    #   estGEP_read<-t(t(estGEP_mat+ t(beta))/colSums(estGEP_mat+ t(beta)))#gene by cell_type
    #
    #   lgpst[itr] = logpost(estProp, estGEP_read, individual_sample[[i]], alpha ,beta)
    #   lglike[itr] = loglikelihood(estProp, estGEP_read, individual_sample[[i]], alpha ,beta)
    #
    #   #  estimates using mcmc state
    #   estProp2<-t(t(estSSp_mat)/colSums(estSSp_mat))  #cell_type by sample_size
    #   estGEP_read2<-t(t(estGEP_mat)/colSums(estGEP_mat))#gene by cell_type
    #
    #
    #   lgpst2[itr] =logpost(estProp2, estGEP_read2, individual_sample[[i]], alpha ,beta)
    #   lglike2[itr] = loglikelihood(estProp2, estGEP_read2, individual_sample[[i]], alpha ,beta)
    #
    #   # store the 3d array
    #   estGEP_array = array(result$csGEP_all,c(cell_type_number, nrow(individual_sample[[i]]), save_last_n_mc),
    #                        dimnames = list(cell_types,gene_symbols,1:save_last_n_mc))
    #   estSSP_array = array(result$SSP_all,c(cell_type_number, samplesize, save_last_n_mc),
    #                        dimnames = list(cell_types,colnames(individual_sample[[i]]),1:save_last_n_mc))
    #
    # }

    result <- gibbsSampler(alpha,beta,individual_sample[[i]],cell_type_number,mcmc_iterations, printout, 0, 0 ,verbose_int, rng, save_last_n_mc)
    #outputs are two vectors. estGEP_vec is gene x cell type; estSSp_vec is sample x cell type
    estGEP_vec<-result$csGEP_vec
    estSSp_vec<-result$SSP_vec

    #vector to matrix
    samplesize<-ncol(individual_sample[[i]])
    estGEP_mat<-t(matrix(estGEP_vec,nrow = cell_type_number,ncol = nrow(individual_sample[[i]]))) # gene by cell_type
    estSSp_mat<-matrix(estSSp_vec,nrow = cell_type_number,ncol = samplesize) #cell_type by sample_size

    # predicted proportions and GEPs
    estProp<-t(t(estSSp_mat+alpha)/colSums(estSSp_mat+alpha))  #cell_type by sample_size
    estGEP_read<-t(t(estGEP_mat+ t(beta))/colSums(estGEP_mat+ t(beta)))#gene by cell_type

    lgpst = logpost(estProp, estGEP_read, individual_sample[[i]], alpha ,beta)
    lglike = loglikelihood(estProp, estGEP_read, individual_sample[[i]], alpha ,beta)

    #  estimates using mcmc state
    estProp2<-t(t(estSSp_mat)/colSums(estSSp_mat))  #cell_type by sample_size
    estGEP_read2<-t(t(estGEP_mat)/colSums(estGEP_mat))#gene by cell_type


    lgpst2 =logpost(estProp2, estGEP_read2, individual_sample[[i]], alpha ,beta)
    lglike2 = loglikelihood(estProp2, estGEP_read2, individual_sample[[i]], alpha ,beta)

    # store the 3d array
    estGEP_array = array(result$csGEP_all,c(cell_type_number, nrow(individual_sample[[i]]), save_last_n_mc),
                         dimnames = list(cell_types,gene_symbols,1:save_last_n_mc))
    estSSP_array = array(result$SSP_all,c(cell_type_number, samplesize, save_last_n_mc),
                         dimnames = list(cell_types,colnames(individual_sample[[i]]),1:save_last_n_mc))


    return(list(estProp = estProp, estGEP = estGEP_read, lgpst = lgpst,lglike = lglike ,
                estProp2 = estProp2, estGEP2 = estGEP_read2, lgpst2 = lgpst2,lglike2 = lglike2 ,
                estSSp_mat = estSSp_mat,estGEP_mat = estGEP_mat,
                estGEP_array = estGEP_array, estSSP_array = estSSP_array,
                cellTypeAssignSplit = result$cellTypeAssignSplit, cellTypeBag = result$cellTypeBag,
                sample_names = colnames(individual_sample[[i]])))

  }#foreach end
  #)#gibbs time
  end.time = Sys.time()
  gibbsRunningTime = end.time - start.time

  #stop parallel
  doParallel::stopImplicitCluster()

  # add sample names to the cdseq_result list
  allsample = 1:length(allresult)
  names(allsample) = lapply(allresult,"[[","sample_names") |> unlist()
  cdseq_result = lapply(allsample,FUN = function(x){
    estProp <- allresult[[x]]$estProp
    estGEP <- allresult[[x]]$estGEP
    lgpst <- allresult[[x]]$lgpst
    lglike <- allresult[[x]]$lglike
    celltype_assignment<-allresult[[x]]$celltype_assignment
    cellTypeAssignSplit <- allresult[[x]]$cellTypeAssignSplit
    cellTypeBag = allresult[[x]]$cellTypeBag
    #processIDs[1] <-allresult[[1]]$processID
    estProp2 <- allresult[[x]]$estProp2
    estGEP2 <- allresult[[x]]$estGEP2
    lgpst2 <- allresult[[x]]$lgpst2
    lglike2 <- allresult[[x]]$lglike2
    estGEP_mat = allresult[[x]]$estGEP_mat
    estSSp_mat = allresult[[x]]$estSSp_mat
    estSSP_array = allresult[[x]]$estSSP_array
    estGEP_array = allresult[[x]]$estGEP_array

    #Final output
    tmp_output <-list(estProp=estProp,
                      estGEP=estGEP,
                      logpost = lgpst,
                      loglikelihood = lglike,
                      estProp2=estProp2,
                      estGEP2=estGEP2,
                      logpost2 = lgpst2,
                      loglikelihood2 = lglike2,
                      gibbsRunningTime = gibbsRunningTime,
                      estGEP_mat = estGEP_mat,
                      estSSp_mat = estSSp_mat,
                      estSSP_array = estSSP_array,
                      estGEP_array = estGEP_array,
                      cell_type_assignment = celltype_assignment,
                      cellTypeAssignSplit = cellTypeAssignSplit,
                      cellTypeBag = cellTypeBag)

    colnames(tmp_output$estGEP)<-cell_types
    rownames(tmp_output$estGEP)<-gene_symbols

    rownames(tmp_output$estProp)<-cell_types
    colnames(tmp_output$estProp)<-allresult[[x]]$sample_names

    colnames(tmp_output$estGEP2)<-cell_types
    rownames(tmp_output$estGEP2)<-gene_symbols

    rownames(tmp_output$estProp2)<-cell_types
    colnames(tmp_output$estProp2)<-allresult[[x]]$sample_names

    dimnames(tmp_output$cellTypeAssignSplit)[[1]] <- gene_symbols
    dimnames(tmp_output$cellTypeAssignSplit)[[2]] <- allresult[[x]]$sample_names
    dimnames(tmp_output$cellTypeAssignSplit)[[3]] <- cell_types

    return(tmp_output)
  })

  # keep all the parameters
  parameters <- list( beta = beta,
                      #alpha = alpha,
                      mcmc_iterations = mcmc_iterations,
                      cpu_number = cpu_number,
                      reference_scRNA = reference_scRNA)


  end_time <- Sys.time()
  cdseq_running_time <- as.numeric( end_time - start_time, units = "hours")
  message(sprintf("cdseq completed successfully using %.4f hours\n",cdseq_running_time))
  return(list(cdseq_result = cdseq_result, parameters = parameters))

}#end function

