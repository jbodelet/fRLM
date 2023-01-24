data {
  int<lower=1> n;
  int<lower=1> J;
  int<lower=1> d;
  vector[n] y;
  matrix[n,J] x;
  matrix[n,d] C;
}


parameters {
  real delta;
  simplex[J] weights;
  vector[d] alpha;
  real<lower=0> sigma;
  real<lower=0> nu;
}

model {
  // the priors 
  delta ~ normal(0, 10);
  weights ~ dirichlet( rep_vector( 1.0/J, J ) );
  sigma ~ lognormal(1,1);
  nu ~ gamma(2, 0.1);
  for ( j in 1:d )
    alpha[j] ~ normal(0, 1);
  //The likelihood
  for (i in 1:n)
    y[i] ~ student_t( nu, delta * x[i] * weights + C[i] * alpha, sigma^2 );
}




