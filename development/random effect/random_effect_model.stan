data {
  int<lower=0> N;              // number of observations
  int<lower=0> T;              // number of time points
  int<lower=0> S;              // number of schools/groups
  int<lower=1> p;              // number of covariates
  int<lower=1,upper=S> Sigma[N,T];   // Matrix of school assignments over time
  matrix[N,p] C;               // Covariate matrix
  vector[N] y;                 // Outcome
}

parameters {
  vector[p] beta;              // Coefficients for covariates
  real alpha;                  // Intercept
  vector[S] xi;                // School effects
  real<lower=0> sigma;         // Error standard deviation
  real<lower=0, upper=1> rho;  // Discounting factor
}

model {
  vector[N] school_effect;

  for (i in 1:N) {
    school_effect[i] = 0;
    for (t in 1:T) {
      school_effect[i] += rho^(T-t) * xi[Sigma[i,t]];
    }
  }

  alpha ~ normal(0, 5);           // Prior for alpha
  beta ~ normal(0, 5);            // Priors for beta
  xi ~ normal(0, 5);              // Priors for school effects
  sigma ~ cauchy(0, 5);           // Prior for sigma

  y ~ normal(C * beta + alpha + school_effect, sigma);  // Likelihood
}
