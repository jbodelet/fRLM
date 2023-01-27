data {
  int<lower=0> n;
  int<lower=0> n_obs;
  int<lower=0> J;
  int<lower=0> K;
  int<lower=1> d;
  vector[n] y;
  vector[n_obs] z_obs; // observed income
  matrix[n, J] x;
  matrix[n, K] u;
  matrix[n,d] C;
}
transformed data {
  int n_mis = n - n_obs;
  matrix[n_obs,K] u_obs = u[1:n_obs];
  matrix[n_mis,K] u_mis = u[(n_obs+1):];
}
parameters {
  real delta;
  simplex[J] weights;
  vector[K] lambda;
  vector[n_mis] z_mis;
  real<lower=0> sigma;
  real<lower=0> sigma_z;
  vector[d] alpha;
  real<lower=0> nu;
}
model {
  vector[n] z = append_row(z_obs, z_mis);
  matrix[n,J] xnew = x;
  xnew[,1] = z + xnew[,1];
  // prior
  delta ~ normal(0, 10);
  weights ~ dirichlet( rep_vector( 1.0/J, J ) );
  sigma ~ lognormal(1,1);
  sigma_z ~ lognormal(1,1);
  for ( k in 1:K )
    lambda[k] ~ normal(0, 1);
  for ( j in 1:d )
    alpha[j] ~ normal(0, 1);
  for(i in 1:n_mis)
    z_mis[i] ~ normal(u_mis[i] * lambda, sigma_z^2);
  nu ~ gamma(2, 0.1);
  // likelihood for observed Y
  for (i in 1:n_obs)
    z_obs[i] ~ normal(u_obs[i] * lambda, sigma_z^2);
  for (i in 1:n)
    y[i] ~ student_t( nu, delta * xnew[i] * weights + C[i] * alpha, sigma^2 );
}






