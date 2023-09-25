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
  int<lower=0> T; // Number of time points for the random effects
  int<lower=0> S; // Number of levels for the grouping variables for the random effects
  int<lower=1, upper=S> Sigma[n, T]; // Matrix of random effect assignments for each period of time
}

parameters {
  vector[p] delta;
  vector[d] alpha;
  simplex[L] beta[p];
  real<lower=0> sigma;
  vector[S-1] xi_unconstrained; // since we have an intercept in C, we constrain the random effects to sum to 0, so we have S-1 degrees of freedom
  real<lower=0, upper=1> rho; // Discounting factor for the random effects (if T>1)
  real<lower=0> sigma_xi; // variance of random effects
}


transformed parameters {
  vector[S] xi;

  for (s in 1:(S-1))
    xi[s] = xi_unconstrained[s];

  xi[S] = -sum(xi_unconstrained); // Ensure the last element makes xi sum to zero
}

model {
  // the priors
  delta ~ normal(0, 10);
  sigma ~ lognormal(0,1);
  sigma_xi ~ lognormal(0,1);
  for (j in 1:d)
    alpha[j] ~ normal(0, alpha_par);
  for (k in 1:p)
    beta[k] ~ dirichlet(rep_vector(beta_par / L, L));

  xi_unconstrained ~ normal(0, sigma_xi);  // Prior for random effects

  // The random effects accumulated over time
  vector[n] random_effects;
  for (i in 1:n) {
    random_effects[i] = 0;
    for (t in 1:T) {
      random_effects[i] += rho^(T-t) * xi[Sigma[i,t]];
    }
  }

  // The likelihood
  vector[n] mu;
  mu = rep_vector(0, n);
  for (k in 1:p) {
    mu += delta[k] * eta[k] * beta[k];
  }
  mu += C * alpha + random_effects;
  for (i in 1:n)
    y[i] ~ normal(mu[i], sigma);
}
