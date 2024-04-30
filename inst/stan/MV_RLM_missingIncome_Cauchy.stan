data {
  int<lower=0> n;
  int<lower=0> p;
  int<lower=0> n_obs; // number of observed income
  int<lower=0> J;
  int<lower=0> K;  // number of predictor variables for income
  int<lower=1> d;
  matrix[n, p] y;
  vector[n_obs] z_obs; // observed income
  matrix[n, J] x;
  matrix[n, K] u;  // predictors
  matrix[n,d] C;
}
transformed data {
  int n_mis = n - n_obs;
  matrix[n_obs,K] u_obs = u[1:n_obs];
  matrix[n_mis,K] u_mis = u[(n_obs+1):];
}
parameters {
  vector[p] delta;
  simplex[J] weights[p];
  vector[d] alpha[p];
  real<lower=0> sigma[p];
  real<lower=0, upper=1> prob;
  real mu_delta;
  real<lower=0> sigma_delta;
  // parameters for missing values model:
    vector[K] lambda;
  vector[n_mis] z_mis;
  real<lower=0> sigma_z;
}
model {
  vector[n] z = append_row(z_obs, z_mis);
  matrix[n,J] xnew = x;
  xnew[,1] = z + xnew[,1];
  // the priors for the main model
  prob ~ uniform(0,1);
  mu_delta ~ cauchy(0, 2.5);
  sigma_delta ~ lognormal(1,1);
  for(k in 1:p )
    weights[k] ~ dirichlet( rep_vector( 1.0/J, J ) );
  sigma ~ lognormal(1,1);
  for(k in 1:p)
    alpha[k] ~ cauchy(0, 2.5);
  // priors for the missing value model:
    sigma_z ~ lognormal(1,1);
  lambda ~ normal(0, 1);
  for(i in 1:n_mis)
    z_mis[i] ~ normal(u_mis[i] * lambda, sigma_z);
  // likelihood for missing value model:
    for(i in 1:n_obs)
      z_obs[i] ~ normal(u_obs[i] * lambda, sigma_z);
  // Likelihood main model:
    for(k in 1:p)
      target += normal_lpdf(y[,k]| delta[k] * xnew * weights[k] + C * alpha[k], sigma[k] ) +
      log_sum_exp( log(prob) +normal_lpdf(delta[k]|0, 0.001), log(1-prob) + normal_lpdf(delta[k]|mu_delta, sigma_delta) );
}





