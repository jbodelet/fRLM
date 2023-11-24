data {
  int<lower=1> n;
  int<lower=1> L;
  int<lower=1> p;
  int<lower=1> d;
  vector[n] y;
  matrix[n,L] eta[p];
  matrix[n,d] C;
  real<lower=0> alpha_par;
  real<lower=0> beta_par;
}

parameters {
  vector[p] delta;
  vector[d] alpha;
  simplex[L] beta[p];
  real<lower=0> sigma;
}

model {
  // the priors
  delta ~ normal(0, 10);
  sigma ~ lognormal(0,1);
  for (j in 1:d)
    alpha[j] ~ normal(0, alpha_par);
  for (k in 1:p)
    beta[k] ~ dirichlet(rep_vector(beta_par / L, L));
  // The likelihood
  vector[n] mu;
  mu = rep_vector(0, n);
  for (k in 1:p) {
    mu += delta[k] * eta[k] * beta[k];
  }
  mu += C * alpha;
  for (i in 1:n)
    y[i] ~ normal(mu[i], sigma);
}
