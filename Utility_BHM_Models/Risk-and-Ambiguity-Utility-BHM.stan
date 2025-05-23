// BAYESIAN HIERARCHICAL EXTENSION TO RISK AND AMBIGUITY UTILITY MODEL (Plate et al., 2025)
//
// Bayesian hierarchical model extending the risk and ambiguity utility model (Gilboa & Schmeidler, 1989) with Luce decision rule (Luce, 1959; Holt & Laury, 2002)
// Built for risky and ambiguous decision-making tasks with 2 alternative forced choice decisions with gain trials
//
// INPUTS
//
// ns             - number of subjects in the analysis
// nt             - largest number of trials completed by any one subject
// resps          - an array of matrices, one matrix per subject, containing the trial-by-trial decisions, with one response per row
//                  the three possible responses are indicated such that:
//                  [1 0] indicates the safe choice
//                  [0 1] indicates the lottery choice
//                  [0 0] indicates a missed trial or that no trial took place
// psure          - a matrix of the probability of the safe option for each subject (column) and trial (row)
// plott          - a matrix of the probability of the lottery option for each subject (column) and trial (row)
// asure          - a matrix of the ambiguity level of the safe option for each subject (column) and trial (row)
// alott          - a matrix of the ambiguity level of the lottery option for each subject (column) and trial (row)
// vsure          - a matrix of the value offered upon selection of the safe option for each subject (column) and trial (row)
// vlott          - a matrix of the value offered upon selection of the lottery option, if won, for each subject (column) and trial (row)
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
// am             - group risk aversion location hyperparameter
// as             - group risk aversion shape hyperparameter
// bm             - group ambiguity aversion location hyperparameter
// bs             - group ambiguity aversion scale hyperparameter
// mm             - group stochasticity location hyperparameter
// ms             - group stochasticity shape hyperparameter
//
// alph_raw       - pre-regularization, untransformed, subject risk aversion parameter
// mu_raw         - pre-regularization, untransformed, subject stochasticity parameter
//
// alph           - subject risk aversion parameter (parameter of interest)
//                  group hyperparameter regularized, anti-logged, alph_raw
// bet            - subject ambiguity aversion parameter
// mu             - subject stochasticity parameter (parameter of interest)
//                  group hyperparameter regularized, anti-logged, mu_raw


functions {

 //Custom function for calculating log likelihood from choices and log odds of each choice
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

 //risk, ambiguity, and value matrices (one subject per column)
 matrix[nt,ns] psure;
 matrix[nt,ns] plott;
 matrix[nt,ns] asure;
 matrix[nt,ns] alott;
 matrix[nt,ns] vsure;
 matrix[nt,ns] vlott;

}

parameters {

 // Hyperparameters
 real am; //alpha "mu" hyperparameter
 real<lower=0> as; //alpha "sigma" hyperparameter
 real bm; //beta "mu" hyperparameter
 real<lower=0> bs; //beta "sigma" hyperparameter
 real mm; //mu "mu" hyperparameter
 real<lower=0> ms; //mu "sigma" hyperparameter

 // Raw Parameters
 vector[ns] alph_raw; //standard normal basis for hyperparameterized lognormal alpha
 vector[ns] mu_raw; //standard normal basis for hyperparameterized lognormal mu

 // Parameters
 vector<lower=-1.35, upper=1.35>[ns] bet; //ambiguity tolerance

}

transformed parameters {

 vector<lower=0>[ns] alph = exp(am + as*alph_raw); //risk tolerance
 vector<lower=0>[ns] mu = exp(mm + ms*mu_raw); //decision stochasticity

}

model {

 // Hyperpriors
 am ~ normal(log(0.65),1.0);
 as ~ lognormal(0.0,1.0);
 bm ~ normal(0.0,1.0);
 bs ~ lognormal(0.0,1.0);
 mm ~ normal(0.0,1.0);
 ms ~ lognormal(0.0,1.0);

 // Raw Priors
 alph_raw ~ normal(0.0,1.0);
 mu_raw ~ normal(0.0,1.0);

 // Priors
 bet ~ normal(bm,bs);

 // Initialize local variables
 matrix[nt,ns] log_svlott;
 matrix[nt,ns] log_svsure;
 matrix[nt,ns] log_pr;

 // Calculate log subjective value difference matrix (number of trials x number of subjects)
 log_svlott = log(plott - rep_matrix(bet',nt).*(alott/2.0)) + rep_matrix(alph',nt).*log(vlott);
 log_svsure = log(psure - rep_matrix(bet',nt).*(asure/2.0)) + rep_matrix(alph',nt).*log(vsure);

 // Calculate log decision rule
 for (n in 1:ns) {
  for (t in 1:nt) {
   log_pr[t,n] = log_svlott[t,n]/mu[n]-log_sum_exp(log_svlott[t,n]/mu[n],log_svsure[t,n]/mu[n]);
  }
 }

 //Likelihood model 
 target += utility_log(resps, log_pr, ns, nt);

}