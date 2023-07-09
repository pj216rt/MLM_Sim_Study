//new elastic net stan file
data {
  //Two sets of data training and testing
  //Training
  int<lower=0> N_obs_train; // number of observations
  int<lower=0> N_pts_train; // number of participants
  
  int<lower=0> K;   // number of level 2 predictors
  int<lower=0> L;   // number of level 1 predictors
  
  int pid_train[N_obs_train]; // participant id vector.  Vector will identify each participant in dataset
  
  matrix[N_obs_train,L] x_train; // level 1 predictors
  matrix[N_pts_train, K] x2_train;   // level 2 predictor matrix
  
  vector[N_obs_train] y_train;      // outcome vector
  
  //testing
  int<lower=0> N_obs_test; // number of observations
  matrix[N_obs_test,K*L] test_data; 
}

parameters {
  vector [L] beta_p [N_pts_train];
  vector<lower=0>[L] tau;      // prior scale
  matrix[K,L] gamma_raw; //regression parameters (level 2)
  cholesky_factor_corr[L] Lcorr; // cholesky factor (L_u matrix for R)
  real<lower=0> sigma2; // population variance
  
  //penalty parameter lambda
  matrix<lower=1>[K, L] tau2;
	real<lower=0> lambda1;
	real<lower=0> lambda2;
}

//We need to work on this
transformed parameters {
  matrix[K,L] gamma;
  gamma = sqrt(((sigma2*(tau2-1))/(lambda2*tau2))) * gamma_raw;
  
  matrix [N_pts_train, L] beta;
  beta = x2_train*gamma;
  
  //create vector from array by looping over the rows
  vector [L] beta_array [N_pts_train];
  for (i in 1:N_pts_train){
    beta_array[i] = to_vector(beta[i]);
  }
  
  corr_matrix[L] R; // correlation matrix
  cov_matrix[L] Sigma; // VCV matrix
  R = multiply_lower_tri_self_transpose(Lcorr); // R = Lcorr * Lcorr'
  Sigma = quad_form_diag(R, tau); // quad_form_diag: diag_matrix(sig) * R * diag_matrix(sig)

  vector[N_obs_train] mu;
  for(i in 1:N_obs_train) {
    mu[i] = (x_train[i] * (beta_p[pid_train[i]])); // * error is on this line
  }
  
  real<lower=0> sigma; //population SD
  sigma = sqrt(sigma2);
}

model {
  to_vector(gamma_raw) ~ normal(0, 1); //implies beta ~ normal(0, sqrt(1/(lambda2/sigma2 * (tau[j]/(tau[j]-1)))))
  to_vector(tau2) ~ gamma(0.5, (8*lambda2*sigma2)/(lambda1^2));
	lambda1 ~ cauchy(0, 1);
	lambda2 ~ cauchy(0, 1);
  
  Lcorr ~ lkj_corr_cholesky(2.0); // prior for cholesky factor of a correlation matrix
  tau ~ inv_gamma(1,7);
  sigma ~ inv_gamma(1,7);
  beta_p ~ multi_normal(beta_array, Sigma);
  y_train ~ normal(mu, sigma);
}

generated quantities {
  vector[N_obs_test] y_new;
  for (n in 1:N_obs_test){
    y_new[n] = normal_rng(test_data[n]*to_vector(gamma), sigma);
  }
}
