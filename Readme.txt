Computational characterization of metacognitive ability in subjective decision-making
Corey R. Plate, Dhruv Govil, Charles Y Zheng, Zoe M. Boundy-Singer, Corey M. Ziemba, Silvia Lopez-Guzman

Description:

This coding package features Stan (cmdstan 2.30.1) code for fitting Bayesian hierarchical extensions to two models of subjective utility, the Risk and Ambiguity Utility Model (Gilboa & Schmeidler, 1989) and the Nonlinear Utility Hyperbolic Discounting Model (Silvia Lopez-Guzman et al., 2018) with Luce decision rule (Luce, 1959; Holt & Laury, 2002), and two versions of a two-stage process model of decision-making and metacognition, the Confidence AS A Noisy Decision Reliability Estimate (CASANDRE) Model (Boundy-Singer et al., 2022), in perceptual- or value-based decision-making.

Contents:

Stan models
 Risk-and-Ambiguity-Utility-BHM.stan
 Nonlinear-Utility-Hyperbolic-Discounting-BHM.stan
 CASANDRE-Perceptual-BHM.stan
 CASANDRE-Value-BHM.stan

Sample data
 Risk-and-Ambiguity-Utility-Sample.mat
 Risk-and-Ambiguity-Utility-Sample.json
 Nonlinear-Utility-Hyperbolic-Discounting-Sample.mat
 Nonlinear-Utility-Hyperbolic-Discounting-Sample.json
 CASANDRE-Perceptual-Sample.mat
 CASANDRE-Perceptual-Sample.json
 CASANDRE-Value-Sample.mat
 CASANDRE-Value-Sample.json

Utility function
 jsonwrite.m

Documentation
 Readme.txt

Usage:

* Risk-and-Ambiguity-Utility-BHM.stan
Fits a model to 2AFC task data where subjects make trial-by-trial decisions between safe options and risky or ambiguous options in the gains domain. The model inputs are the task parameters and subject decisions saved as a structural variable in a .json file. The model outputs are one .csv file for each sampling chain where each column is a diagnostic or a posterior over subject-level risk tolerance, ambiguity tolerance, and choice stochasticity parameters, and the group hyperparameters.

* Nonlinear-Utility-Hyperbolic-Discounting-BHM.stan
Fits a model to 2AFC task data where subjects make trial-by-trial decisions between immediate rewards and rewards with a delay in the gains domain. The model inputs are the task parameters and choices saved as a structural variable in a .json file. The model outputs are one .csv file for each sampling chain where each column is a diagnostic or a posterior over subject-level discounting rate, and choice stochasticity parameters, and the group hyperparameters.

* CASANDRE-Perceptual-BHM.stan / CASANDRE-Value-BHM.stan
Fits a model to 2AFC and confidence task data, where subjects make trial-by-trial decisions and binary confidence reports. The model inputs are the objective or subjective magnitude of the stimulus (e.g. orientation, subjective value difference) and the subject decisions and confidence ratings saved as a structural variable in a .json file. The model outputs are one .csv for each sampling chain where each column is a diagnostic or a posterior over subject-level lapse rate, decision sensitivity, decision criterion, meta-uncertainty, and a confidence criterion for each decision option, and the group hyperparameters.

The value-based version of the CASANDRE BHM is used when the magnitude of the stimulus is obtained by fitting a model that gives the internal decision variable. The only difference in model output is that the model does not fit a decision criterion, which is fixed at zero.

Both of these models have been parallelized for multithreading to promote tractable run times.
