data {
  int<lower=1> n;
  int<lower=1> p;
  int<lower=1> J;
  int<lower=1> d;
  matrix[n,p] y;
  matrix[n,J] x;
  matrix[n,d] C;
}


parameters {
  vector[p] delta;
  simplex[J] weights[p];
  vector[d] alpha[p];
  real<lower=0> sigma[p];
  real<lower=0, upper=1> prob;
  real mu_delta;
  real<lower=0> sigma_delta;
}

model {
  // the priors 
  prob ~ uniform(0,1);
  mu_delta ~ normal(0, 10);
  sigma_delta ~ lognormal(1,1);
  for(k in 1:p )
    weights[k] ~ dirichlet( rep_vector( 1.0/J, J ) );
  sigma ~ lognormal(1,1);
  for(k in 1:p)
    alpha[k] ~ normal(0, 1);
  //The likelihood
  for(k in 1:p)
    target += normal_lpdf(y[,k]| delta[k] * x * weights[k] + C * alpha[k], sigma[k] ) +
    log_sum_exp( log(prob) +normal_lpdf(delta[k]|0, 0.001), log(1-prob) + normal_lpdf(delta[k]|mu_delta, sigma_delta) );
}
