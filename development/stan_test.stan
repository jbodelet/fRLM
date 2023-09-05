data {
  int<lower=0> N; // number of observations
  vector[N] x; // predictor
  vector[N] y; // outcome
}

parameters {
  real alpha; // intercept
  real beta; // slope
  real<lower=0> sigma; // standard deviation of the noise
}

model {
  y ~ normal(alpha + beta*x, sigma);
}
