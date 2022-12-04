data {
  int<lower=1> n;
  int<lower=1> L;
  int<lower=1> d;
  vector[n] y;
  matrix[n,L] eta;
  real<lower=0> beta_par;
  real<lower=0> alpha_par;
  matrix[n,d] C;
}


parameters {
  real delta;
  vector[d] alpha;
  simplex[L] beta;
  real<lower=0> sigma;
}

model {
  // the priors 
  delta ~ normal(0, 10);
  beta ~ dirichlet( rep_vector( beta_par / L, L ) );
  sigma ~ lognormal(0,1);
  for ( j in 1:d )
    alpha[j] ~ normal(0, alpha_par);
  //The likelihood
  for (i in 1:n)
    y[i] ~ normal( delta * eta[i] * beta + C[i] * alpha, sigma^2 );
}



