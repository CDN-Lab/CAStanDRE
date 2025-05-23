// BAYESIAN HIERARCHICAL EXTENSION TO NONLINEAR UTILITY HYPERBOLIC DISCOUNTING MODEL (Plate et al., 2025)
//
// Bayesian hierarchical model extending the nonlinear utility hyperbolic discounting model (Silvia Lopez-Guzman et al., 2018) with Luce decision rule (Luce, 1959; Holt & Laury, 2002)
// Built for delay discounting tasks with 2 alternative forced choice decisions
//
// INPUTS
//
// ns             - number of subjects in the analysis
// nt             - largest number of trials completed by any one subject
// resps          - an array of matrices, one matrix per subject, containing the trial-by-trial decisions, with one response per row
//                  the three possible responses are indicated such that:
//                  [1 0] indicates the sooner choice
//                  [0 1] indicates the later choice
//                  [0 0] indicates a missed trial or that no trial took place
// vnow           - a matrix of the value of the sooner option for each subject (column) and trial (row)
// vlater         - a matrix of the value of the later option for each subject (column) and trial (row)
// del            - a matrix of the delay period for each subject (column) and trial (row)
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
// km             - group discount rate location hyperparameter
// ks             - group discount rate shape hyperparameter
// mm             - group stochasticity location hyperparameter
// ms             - group stochasticity shape hyperparameter
//
// kapp_raw       - pre-regularization, untransformed, subject discount rate parameter
// mu_raw         - pre-regularization, untransformed, subject stochasticity parameter
//
// kapp           - subject delay discount parameter (parameter of interest)
//                  group hyperparameter regularized, anti-logged, kapp_raw
// mu             - subject stochasticity parameter (parameter of interest)
//                  group hyperparameter regularized, anti-logged, mu_raw

functions {

 //Custom function for calculating log likelihhood from choices and log odds of each choice
 real utility_log(array[] matrix resps, matrix log_pr, int ns, int nt) {

  //initialize log likelihood
  real llhC;
  llhC = 0.0;

  //cycle through participants
  for (n in 1:ns) {

   //initialize choice probability and likelihood
   vector[nt] log_p;
   matrix[nt,2] llh;

   //log trap
   for (t in 1:nt) {
    if (log_pr[t,n] == log(1.0)) {
     log_p[t] = log1m(machine_precision());
    } else if (log_pr[t,n] == log(0.0)) {
     log_p[t] = log(machine_precision());
    } else {
     log_p[t] = log_pr[t,n];
    }
   }

   //build likelihood matrix
   llh = append_col(log1m_exp(log_p),log_p);

   //calculate and sum log likelihood
   llhC += sum(columns_dot_product(resps[n],llh));

  }

  //return total log likelihood
  return llhC;

 }

}

data {

 //number of participants
 int ns;

 //number of trials
 int nt;

 //response matrices
 array[ns] matrix[nt,2] resps;

 //value later, value now, and delay matrices (one subject per column)
 matrix[nt,ns] vlater;
 matrix[nt,ns] vnow;
 matrix[nt,ns] del;

 //Utility-CRDM-BHM model fitted alpha values
 vector[ns] alph;

}

parameters {

 // Hyperparameters
 real km; //kappa "mu" hyperparameter
 real<lower=0> ks; //kappa "sigma" hyperparameter
 real mm; //mus "mu" hyperparameter
 real<lower=0> ms; //mu "sigma" hyperparameter

 // Parameters
 vector[ns] kapp_raw; //standard normal basis for hyperparameterized lognormal kappa
 vector[ns] mu_raw; //standard normal basis for hyperparameterized lognormal mu

}

transformed parameters {

 vector<lower=0>[ns] kapp = exp(km + ks*kapp_raw); //discounting rate
 vector<lower=0>[ns] mu = exp(mm + ms*mu_raw); //decision stochasticity

}

model {

 // Hyperpriors
 km ~ normal(-3.0,1.0);
 ks ~ lognormal(0.0,1.0);
 mm ~ normal(0.0,1.0);
 ms ~ lognormal(0.0,1.0);

 // Raw Priors
 kapp_raw ~ normal(0.0,1.0);
 mu_raw ~ normal(0.0,1.0);

 // Initialize local variables
 matrix[nt,ns] log_svlater;
 matrix[nt,ns] log_svnow;
 matrix[nt,ns] log_pr;

 // Calculate log subjective value difference matrix (number of trials x number of subjects)
 log_svlater = rep_matrix(alph',nt).*log(vlater) - log1p(rep_matrix(kapp',nt).*del);
 log_svnow = rep_matrix(alph',nt).*log(vnow);

 // Calculate log decision rule
 for (n in 1:ns) {
  for (t in 1:nt) {
   log_pr[t,n] = log_svlater[t,n]/mu[n]-log_sum_exp(log_svlater[t,n]/mu[n],log_svnow[t,n]/mu[n]);
  }
 }

 //Likelihood model 
 target += utility_log(resps, log_pr, ns, nt);

}