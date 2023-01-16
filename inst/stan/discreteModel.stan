data {
  int<lower=1> n;
  int<lower=1> J;
  vector[n] y;
  matrix[n,J] x;
}


parameters {
  real delta;
  simplex[J] weights;
  real<lower=0> sigma;
}

model {
  // the priors 
  delta ~ normal(0, 10);
  weights ~ dirichlet( rep_vector( 1, J ) );
  sigma ~ lognormal(0,1);
  //The likelihood
  for (i in 1:n)
    y[i] ~ normal( delta * x[i] * weights, sigma^2 );
}




