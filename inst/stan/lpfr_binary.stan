data {
  int<lower=1> n;
  int<lower=1> L;
  int<lower=1> K;
  int<lower=1> d;
  int<lower=1> Nmax;
  int<lower=0,upper=1> y[n];
  matrix[n,d] C;
  matrix[Nmax,K] phi_mat[n];
  matrix[K, L] J;
  // DATA for gaussian processes:
  int<lower=1> Nvec[n];
  real tobs[n, Nmax];
  vector[Nmax] xobs[n];
}
parameters {
  real delta;
  vector[d] alpha;
  simplex[L] beta;
  real<lower=0> sigma;
  real<lower=0> sigma_x;
  // GP parameters:
  vector<lower=0, upper=1>[K] xi[n];
}

model {
  for(i in 1:n){
    xi[i] ~ uniform(0,1);
  }
  sigma_x ~ lognormal(0,1);
  // fRLM priors:
  delta ~ normal(0, 10);
  beta ~ dirichlet( rep_vector( 1.0 / L, L ) );
  alpha ~ std_normal();
  // GP liklihood:
  for(i in 1:n ){
      xobs[i,:Nvec[i] ] ~ normal(phi_mat[i][:Nvec[i] ] * xi[i], sigma_x ); // try beta distribution
  }
  // the likelihood
  for (i in 1:n) {
    real lin_pred = delta * xi[i]' * J * beta + C[i] * alpha;
    // print("Linear predictor for ", i, ": ", lin_pred);
    y[i] ~ bernoulli_logit(lin_pred);
  }

}
