// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>

#include <Rcpp.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>
//#include <numeric> // will use std::accumulate function, for sum up the elements in mixtureSamples

#include "cokus.h"
#include "Xoshiro256plus.h"

using namespace std;
using namespace Rcpp;

//' This is the Gibbs sampler.
//' \code{GibbsSampler} returns estimated GEPs and cell type proportions.
//' @param ALPHA hyperparameter for cell type proportion.
//' @param BETA hyperparameter for cell-type-specific GEPs.
//' @param mixtureSamples bulk RNA-seq data in form of read counts.
//' @param n_types number of cell types.
//' @param NN number of MCMC iteration.
//' @param OUTPUT MCMC progress output control.
//' @param initialization if 1, initialize cell type assignment with provided value, if 0, initialize cell type assignment with random values.
//' @param CT initial values for cell type assignment.
//' @param verbose if greater than or euqal to 1, then print working progress in console, otherwise do not print in console.
//' @param RNG choose one of the random number generators, 1 for MT, 2 for Xoshiro256plus
//' @param NMC number of MCMC samples tos save
//' @return random integers  uniformly distributed in 0..(2^32 - 1).
// [[Rcpp::export]]
List gibbsSampler(NumericMatrix ALPHA, NumericMatrix BETA, NumericMatrix mixtureSamples, int n_types, int NN, int OUTPUT, int initialization, std::vector<unsigned int> CT, int verbose, int RNG, int NMC)
 {
   //clock_t start,finish;
   //start=clock();
   /*
    *  initialize variables
    */
   int genes = mixtureSamples.rows(), samples = mixtureSamples.cols();
   int n = 0;
   //n = std::accumulate(mixtureSamples.begin(),mixtureSamples.end(),0);
   vector<unsigned int> sid, gid;// unsigned int for saving space
   for(int i = 0; i < genes; i ++)
   {
     for(int j = 0; j < samples; j ++)
     {
       n += mixtureSamples(i, j);
       for(int k = 0; k < mixtureSamples(i, j); k ++)
       {
         sid.push_back(j);
         gid.push_back(i);
       }
     }
   }

   //vector<unsigned int>  cellTypeTot(n_types),csGEP_vec(genes*n_types), SSP_vec(n_types*samples);
   vector<double> probs(n_types);
   vector<unsigned int> order(n), cellTypeBag(n);// unsigned int for saving space

   vector<double> WBETA(n_types);

   arma::uvec  cellTypeTot(n_types), csGEP_vec(genes*n_types), SSP_vec(n_types*samples);
   cellTypeTot.fill(0);csGEP_vec.fill(0);SSP_vec.fill(0);

   // this matrix stores reads-cell-type-assignment for each sample
   // use RcppArmadillo cube data structure as 3D matrix container
   arma::cube cellTypeAssignSplit(genes,samples,n_types);
   cellTypeAssignSplit.fill(0);// initialize with zeros
   // keep all mcmc estimates
   arma::umat csGEP_all(genes*n_types, NMC);
   csGEP_all.fill(0);
   arma::umat SSP_all(n_types*samples, NMC);
   SSP_all.fill(0);

   int wi,di,i,ii,j,cellType, rp, iter, wioffset, dioffset;
   double totprob, r, max;
   //ofstream logfile (CDSeq_tmp_log, std::ios_base::app);
   //auto time_now = std::chrono::system_clock::now();
   //std::time_t cdseq_time = std::chrono::system_clock::to_time_t(time_now);

   // seeding
   int SEED = 3;
   seedMT( 1 + SEED * 2 ); // seeding only works on uneven numbers
   xoshiro256seed(1);

   for (i = 0; i < n; i ++)
   {
     //cout << i << "/" << n << endl;
     wi = gid[i];
     di = sid[i];
     // pick a random cell type 0..n_types-1
     if(initialization==0){
        if(RNG==1) {cellType = (int) ( (double) randomMT() * (double) n_types / (double) (4294967296.0 + 1.0) );} // To prevent overflow
        if(RNG==2) {cellType = (int) ( (double) xoshiro256rng() * (double) n_types) / (double) (1.8446744e+19);} // To prevent overflow
        cellTypeBag[i] = cellType; // assign this gene read to this cell type
     }else{
        // initialize cell type assignment
        cellType = (int) CT[i];
        cellTypeBag[i] = cellType;
     }
     csGEP_vec[wi*n_types + cellType]++;
     SSP_vec[di*n_types + cellType]++;
     cellTypeTot[cellType]++; // increment cellTypeTot matrix
     cellTypeAssignSplit(wi,di,cellType)++;
     order[i]=i;
   }
   //for (i = 0; i < n; i ++) order[i]=i; // fill with increasing series

   // Fisher-Yates Shuffle Algorithm
   for (i = 0; i < (n-1); i ++) {
     // pick a random integer between i and nw
     if(RNG==1){rp = i + (int) ((double) (n-i) * (double) randomMT() / (double) (4294967296.0 + 1.0));}
     if(RNG==2){rp = i + (int) ((double) (n-i) * (double) xoshiro256rng() / (double) (1.8446744e+19) );}

     // switch contents on position i and position rp
     swap(order[rp], order[i]);
   }

   //WBETA = (double) (accumulate(BETA.begin(), BETA.end(), 0));

   for(i=0; i<n_types; i++){
     NumericVector x = BETA( i, _ ); //copy the i th row
     // I had a bug for large precision and it turned out here is where it goes wrong
     // the std lib cpp fct accumulate is instantiated with the type of the third argument which was used 0 and it is an type of int
     // this caused overflow issues for large BETA values and accumulate returns negative values for large BETA values
     // current solution is to use 0.0 instead of 0 to make it a double type
     WBETA[i] = (double) (accumulate(x.begin(), x.end(), 0.0));
     //Rcout<<"i="<<i<<"; WBETA[i]="<< WBETA[i]<<endl;

   }

   for (iter = 0; iter < NN; iter ++)
   {
     Rcpp::checkUserInterrupt();
     if (OUTPUT >= 1 && verbose>=1) {
       if (((iter+1) % 10) == 0) Rprintf( "%-4d of %-4d MCMC runs (finished %2d%%, press esc at any time to terminate the MCMC iteration)\r" , iter+1 , NN , 100*(iter+1)/NN);
       if (iter == NN-1) Rprintf( "%-4d of %-4d MCMC runs (finished %2d%%)\nGibbs sampler completed succefully\n" , iter + 1 , NN, 100*(iter+1)/NN );
     }
     //if(OUTPUT==0){
       //if (iter % 10 == 0){cout<<"Running in parallel:"<< 100*iter/NN<<"% of the job is finished\r";}
       //if (iter % 10 == 0){Rcpp::Rcout<<"Cell type number is"<<setw(4)<<n_types<<", data block ["<<setw(4)<<anatomical_idx<<"] finished "<<setw(4)<<100*iter/NN<<"% of MCMC iterations (worker process ID "<<setw(7)<<processID<<")\r";}
       //if(iter == NN){Rcpp::Rcout<<"Cell type number is"<<setw(4)<<n_types<<", data block ["<<setw(4)<<anatomical_idx<<"] finished "<<setw(4)<<100*iter/NN<<"% of MCMC iterations (worker process ID "<<setw(7)<<processID<<")"<<std::endl;}
       //if (logfile.is_open()){
        // if (iter == 0){logfile<<std::ctime(&cdseq_time);}
        // if (iter % 10 == 0){logfile<<"Cell type number is"<<setw(4)<<n_types<<", data block ["<<setw(4)<<anatomical_idx<<"] finished "<<setw(4)<<100*iter/NN<<"% of MCMC iterations (worker process ID "<<setw(7)<<processID<<")\n";}
        // if (iter == NN){logfile << "\n\n";}
         //if(iter == NN){logfile<<"Cell type number is"<<setw(4)<<n_types<<", data block ["<<setw(4)<<anatomical_idx<<"] finished "<<setw(4)<<100*iter/NN<<"% of MCMC iterations (worker process ID "<<setw(7)<<processID<<")\n";}
       //}
     //}

     for (ii = 0; ii < n; ii++) {
       i = order[ii]; // current gene read to assess
       wi  = gid[i]; // current gene index
       di  = sid[i]; // current sample index
       cellType = cellTypeBag[i]; // current cell type assignment to gene read
       cellTypeTot[cellType]--;  // substract this from counts
       cellTypeAssignSplit(wi,di,cellType)--;

       wioffset = wi*n_types;
       dioffset = di*n_types;

       csGEP_vec[wioffset+cellType]--;
       SSP_vec[dioffset+cellType]--;

       totprob = (double) 0;
       for (j = 0; j < n_types; j++) {
         probs[j] = ((double) csGEP_vec[wioffset + j] + (double) BETA(j,wi))/( (double) cellTypeTot[j]+ (double) WBETA[j])*( (double) SSP_vec[dioffset+ j] + (double) ALPHA(j,di));
         totprob += probs[j];
       }
       // sample a cell type from the distribution
       if(RNG==1){r = (double) totprob * (double) randomMT() / (double) 4294967296.0;}
       if(RNG==2){r = (double) totprob * (double) xoshiro256rng() / (double) 1.8446744e+19;}
       max = probs[0];
       cellType = 0;
       while (r > max) {
         cellType ++;
         max += probs[cellType];
       }

       cellTypeBag[i] = cellType; // assign current gene read i to cell type j
       csGEP_vec[wioffset + cellType ]++; // and update counts
       SSP_vec[dioffset + cellType ]++;
       cellTypeTot[cellType] ++;
       cellTypeAssignSplit(wi,di,cellType)++;//keep read-cell-type-assignment for each cell type in each sample
     }

     // keep NMC mcmc runs
     if(iter > NN - NMC - 1){
        csGEP_all.col(iter - NN + NMC) = csGEP_vec;
        SSP_all.col(iter - NN + NMC) = SSP_vec;
     }
   }

   //logfile.close();
   List result;
   result["csGEP_vec"] = csGEP_vec;
   result["SSP_vec"] = SSP_vec;
   result["cellTypeTot"] = cellTypeTot;
   result["cellTypeAssignSplit"] = cellTypeAssignSplit;
   result["cellTypeBag"] = cellTypeBag;
   result["SSP_all"] = SSP_all;
   result["csGEP_all"] = csGEP_all;
   return result;
 }






