data {
  int<lower=1> n;
  int<lower=1> L;
  int<lower=1> d;
  int<lower=0,upper=1> y[n];
  matrix[n,L] eta;
  real<lower=0> beta_par;
  real<lower=0> alpha_par;
  matrix[n,d] C;
}


parameters {
  real delta;
  vector[d] alpha;
  simplex[L] beta;
}

model {
  // the priors 
  delta ~ normal(0, 10);
  beta ~ dirichlet( rep_vector( beta_par / L, L ) );
  for ( j in 1:d )
    alpha[j] ~ normal(0, alpha_par);
  //The likelihood
    y ~ bernoulli_logit( delta * eta * beta + C * alpha );
}



