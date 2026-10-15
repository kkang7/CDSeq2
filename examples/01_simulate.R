# Entirely synthetic teaching data; no participant data are used.
simulate_tutorial <- function() {
  set.seed(20261015)
  genes <- sprintf('Gene_%03d', seq_len(240))
  types <- c('Excitatory neurons', 'Inhibitory neurons', 'Astrocytes')
  mu <- matrix(rgamma(720, shape=1.5), 240, 3, dimnames=list(genes, types))
  for(k in 1:3) mu[((k-1)*60+1):(k*60),k] <- mu[((k-1)*60+1):(k*60),k]*8
  mu <- sweep(mu,2,colSums(mu),'/')
  annotation <- rep(types,each=12)
  reference <- vapply(seq_along(annotation),function(j)
    as.numeric(rmultinom(1,4000,mu[,annotation[j]])),numeric(240))
  dimnames(reference) <- list(genes,sprintf('Reference_cell_%02d',1:36))
  reference_donors <- rep(sprintf('Reference_donor_%02d',1:6),6)
  reference_mean <- vapply(types,function(ct)rowSums(reference[,annotation==ct,drop=FALSE]),numeric(240))
  reference_mean <- sweep(reference_mean,2,colSums(reference_mean),'/')
  dimnames(reference_mean) <- list(genes,types)
  make_cohort <- function(label,n=20) {
    samples <- sprintf('%s_%02d',label,seq_len(n))
    phenotype <- as.numeric(scale(rnorm(n)))
    genotype <- rbinom(n,2,.35)
    counts <- array(0,c(240,3,n),list(genes,types,samples))
    theta <- matrix(rgamma(n*3,shape=4),n,3,dimnames=list(samples,types))
    theta <- theta/rowSums(theta)
    for(i in seq_len(n)) for(k in 1:3) {
      z <- log(mu[,k])+rnorm(240,0,runif(240,.1,.65))
      if(k==1) {z[1:8] <- z[1:8]+.4*phenotype[i];z[9:16] <- z[9:16]+.45*(genotype[i]-.7)}
      counts[,k,i] <- as.numeric(rmultinom(1,round(30000*theta[i,k]),exp(z)))
    }
    list(bulk=apply(counts,c(1,3),sum),truth=counts,
      donor_ids=setNames(samples,samples),
      phenotype=setNames(phenotype,samples),genotype=setNames(genotype,samples))
  }
  list(reference=reference,annotation=annotation,reference_donors=reference_donors,
    reference_mean=reference_mean,training=make_cohort('Training'),validation=make_cohort('Validation'))
}
