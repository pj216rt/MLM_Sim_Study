The purpose of this project is to derive a novel method to identify and estimate latent variables in a two level longitudinal model.  

This work is comprised of four main .R files.  The first serves as a repository for functions that will be used throughout the project.  
The second generates simulated data.  The third runs STAN sampling on these generated datasets using the cmdstan R package.  
The fourth analyzes the sampled results.

Our belief is that subjects' individual growth patterns arise from a multivariate Gaussian distribution (level 1) governed by 
population-specific parameters modeled via regression on level 2 explanatory covariates.  We are focused on a fully Bayesian 
approach to accurately identify variables.  Primarily, we are motivated by the work of van Erp, Oravecz and Muth, and guided by the 
STAN development team.

A large part of our work involves evaluating different priors in our models.  To that end, we consider multiple 
different priors including but not limited to ridge, LASSO, horseshoe, elastic net, and discrete normal mixtures.
 