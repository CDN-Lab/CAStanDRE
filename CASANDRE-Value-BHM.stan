// BAYESIAN HIERARCHICAL EXTENSION TO CASANDRE MODEL FOR FITTING CONFIDENCE PARAMETERS IN PREFERENCE-BASED TASKS (Plate et al., 2025)
//
// Bayesian hierarchical model extending the Confidence as a Noisy Decision Reliability Estimate (CASANDRE) likelihood model (Boundy-Singer et al., 2022)
// Modified for preference-based tasks with 2 alternative forced choice decisions and trial-by-trial binary confidence reports
//
// INPUTS
//
// ns             - number of subjects in the analysis
// nt             - largest number of trials completed by any one subject
// resps          - an array of matrices, one matrix per subject, containing the trial-by-trial decisions and confidence reports, with one response per row
//                  the five possible responses are indicated such that:
//                  [1 0 0 0] indicates high confidence in an option 1 choice
//                  [0 1 0 0] indicates low confidence in an option 1 choice
//                  [0 0 1 0] indicates low confidence in an option 2 choice
//                  [0 0 0 1] indicates high confidence in an option 2 choice
//                  [0 0 0 0] indicates a missed trial or that no trial took place
// vals           - a matrix of the subjective value differences between the two options (option 2 - option 1) for each subject (column) and trial (row)
// sampn          - an integer declaring the number of samples of cumulative density used to calculate the cutoffs for each choice-confidence pair
//                  the higher the number of samples, the more accurate
//                  for tractable run times, suggest 30
// sampx          - equally spaced cumulative densities between 0 and 1, the number of samples of which is declared by sampn
//
// OUTPUTS
//
// lp__           - log density up to a constant
// accept_stat__  - Metropolis acceptance probability
// stepsize__     - numerical integration integrator step size
// treedepth__    - depth of binary tree used by No-U-Turn Sampler algorithm
// n_leapfrog__   - number of leapfrog steps utilized in numerical integration
// divergent__    - number of divergent transitions that occured during sampling
// energy__       - the value of the Hamiltonian function at the state sampled for the given iteration
//
// ssm            - group decision sensitivity location hyperparameter
// sss            - group decision sensitivity shape hyperparameter
// mum            - group meta-uncertainty location hyperparameter
// mus            - group meta-uncertainty shape hyperparameter
// nccm           - group option 1 confidence criterion location hyperparameter
// nccs           - group option 1 confidence criterion shape hyperparameter
// pccm           - group option 2 confidence criterion location hyperparameter
// pccs           - group option 2 confidence crtierion shape hyperparameter
//
// guess          - subject lapse rate parameter
// sens           - subject decision sensitivity parameter
// meta           - subject meta-uncertainty parameter
// nconf          - subject option 1 confidence criterion parameter
// pconf          - subject option 2 confidence criterion parameter

functions {

 //Custom quantile function for lognormal distribution
 matrix log_normal_qf(vector x, vector m, vector s) {

  //takes a vector sampling the x-axis in terms of cumulative density, along with vectors of mu and sigma parameters
  return exp(rep_matrix(s' * sqrt(2), size(x)) .* rep_matrix(inv_erfc(-2 * x + 2), size(m)) + rep_matrix(m', size(x)));

 }

 //Likelihood model for CASANDRE
 real casandre_log(array[] matrix resps, vector guess, matrix sds, matrix sm, vector nconf, vector pconf) {

  //declare and initialize log likelihood variable
  real llhC;
  llhC = 0;

  //loop over participants to calculate log likelihood
  for (n in 1:size(resps)) {

   // declare the index-sizing variable
   int m[rows(resps[n])];

   // build the logical vector of non-zero trials
   for (tt in 1:size(m)) {
      m[tt] = sum(resps[n][tt]) == 1;
   }

   //declare the likelihood variable
   matrix[rows(sm),4] lhC;

   // check if subject has any non-zero data before running
   if (sum(resps[n]) != 0) {

    // declare incrementing and index variables
    int t;
    int ind[sum(m)];

    // initialize incrementing variable
    t = 1;

    //declare and calculate mu parameters for normal cdfs
    matrix[rows(sm),rows(sds)] avgs;
    avgs = rep_matrix(sm[,n],rows(sds)).*rep_matrix(sds[,n]',rows(sm));

    //loop over trials  
    for (tr in 1:rows(sm)){

     // check if trial has any non-zero data before running
     if (sum(resps[n][tr]) != 0) {

      //declare sample vector
      matrix[3,rows(sds)] raws;

      //sample the cumulative density of normals along the transformed x-axis for each confidence bin
      for (rws in 1:rows(sds)) {
       raws[1,rws] = normal_cdf(-nconf[n],avgs[tr,rws],sds[rws,n]);
      }
      for (rws in 1:rows(sds)) {
       raws[2,rws] = normal_cdf(0,avgs[tr,rws],sds[rws,n]);
      }
      for (rws in 1:rows(sds)) {
       raws[3,rws] = normal_cdf(pconf[n],avgs[tr,rws],sds[rws,n]);
      }

      //declare the cumulative likelihood variable
      vector[3] ratiodist;

      //take the mean of the sampled densities
      ratiodist[1] = mean(raws[1,:]);
      ratiodist[2] = mean(raws[2,:]);
      ratiodist[3] = mean(raws[3,:]);

      //implement a lapse rate parameter to absorb unmodeled noise, and calculate likelihoods of each choice possibility in the trial
      lhC[tr,1] = ((guess[n]/4) + ((1-guess[n])*ratiodist[1]));
      lhC[tr,2] = ((guess[n]/4) + ((1-guess[n])*(ratiodist[2]-ratiodist[1])));
      lhC[tr,3] = ((guess[n]/4) + ((1-guess[n])*(ratiodist[3]-ratiodist[2])));
      lhC[tr,4] = ((guess[n]/4) + ((1-guess[n])*(1-ratiodist[3])));
 
      // add this trial to the index
      ind[t] = tr;

      // increment the index
      t = t + 1;

     }

    }

   //multiply the choice matrix by the log of the likelihood matrix
   llhC += sum(columns_dot_product(resps[n][ind], log(lhC[ind])));

   }
  
  }

  //return the log likelihood for all data given the sampled parameters
  return llhC;
  
 }

 //Partial sum function
 real partial_sum_casandre_log(array[] matrix slice_n_resps, int start, int end, vector guess, matrix sds, matrix sm, vector nconf, vector pconf) {

 //dynamically slice data into the log likelihood function
 return casandre_log(slice_n_resps, guess[start:end], sds[:,start:end], sm[:,start:end], nconf[start:end], pconf[start:end]);

 }
  
}

data {

 // declare the number of subjects
 int ns;

 // declare the number of trials in the largest unpadded data set
 int nt;

 // declare the response matrices for value-based domain
 array[ns] matrix[nt,4] resps;

 // declare the experimental values matrices for value-based domain
 matrix[nt,ns] vals;

 // declare the number of x-axis samples and the vector sampling the x-axis in terms of cumulative density
 int sampn;
 vector[sampn] sampx;
 
}

parameters {

 //Parameters
 vector<lower=0,upper=1>[ns] guess; //lapse rate
 vector<lower=0>[ns] sens; //stimulus sensitivities
 vector<lower=0>[ns] meta; //meta-uncertainty
 vector<lower=0>[ns] nconf; //negative confidence criteria
 vector<lower=0>[ns] pconf; //positive confidence criteria

 //Hyperparameters
 real ssm; //mu hyperparameter for stimulus sensitivity lognormal prior
 real<lower=0> sss; //sigma hyperparameter for stimulus sensitivity lognormal prior
 real mum; //mu hyperparameter for meta-uncertainty lognormal prior
 real<lower=0> mus; //sigma hyperparameter for meta-uncertainty lognormal prior
 real nccm; //mu hyperparameter for negative stimulus criterion lognormal prior
 real<lower=0> nccs; //sigma hyperparameter for negative stimulus criterion lognormal prior
 real pccm; //mu hyperparameter for positive stimulus criterion lognormal prior
 real<lower=0> pccs; //sigma hyperparameter for positive stimulus criterion lognormal prior

}

model {

 // Hyperpriors
 ssm ~ normal(0,1);
 sss ~ lognormal(0,1);
 mum ~ normal(0,1);
 mus ~ lognormal(0,1);
 nccm ~ normal(0,1);
 nccs ~ lognormal(0,1);
 pccm ~ normal(0,1);
 pccs ~ lognormal(0,1);

 // Priors
 guess ~ beta(1,197.0/3.0);
 sens ~ lognormal(ssm,sss);
 meta ~ lognormal(mum,mus);
 nconf ~ lognormal(nccm,nccs);
 pconf ~ lognormal(pccm,pccs);

 //declare the transformed x-axis matrices
 matrix[sampn,ns] xtrans;
 matrix[sampn,ns] sds;

 //calculate the transformed x-axis matrix
 xtrans = log_normal_qf(sampx,log(1/sqrt(meta.^2+1)),sqrt(log(meta.^2+1)));
 sds = 1./xtrans;

 // Likelihood model (with local variable calculations)

 //declare matrices for calculating data transformations
 matrix[nt,ns] om;
 matrix[nt,ns] sm;

 //orientations in terms of standard deviations
 om =  vals.*rep_matrix(sens',nt);

 //transform the data (kept for cross-code consistency)
 sm = om;
 
 //hand the data and parameters to the reduce sum wrap function for dynamic within-chain parallelization
 target += reduce_sum(partial_sum_casandre_log, resps, 1, guess, sds, sm, nconf, pconf);

}