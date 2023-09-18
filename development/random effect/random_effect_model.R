library(rstan)

# Simulation parameters
N <- 500  # Number of students
T <- 5    # Number of time points
S <- 7    # Number of schools
p <- 2    # Number of covariates

# True values for the parameters
alpha_true <- c(2, -0.5)
xi_true <- c(1.2, -0.8, 0.5, 0.2, 1.8, -1.4, -0.1)
rho_true <- 0.8
sigma_true <- 1

# Simulating data
set.seed(42)

# Covariates: Let's assume they're drawn from a standard normal distribution
C <- matrix(rnorm(N * p), N, p)

# School assignments: Randomly assign students to schools over time
Sigma <- matrix(sample(1:S, N*T, replace = TRUE), N, T)

# Simulating outcome y based on the model
# Simulating outcome y based on the model
mu <- matrix(0, N, T)

for (i in 1:N) {
  for (t in 1:T) {
    mu[i,t] = rho_true^(T-t) * xi_true[Sigma[i,t]]
  }
}


y <- apply(mu, 1, sum) + (C %*% alpha_true) + rnorm(N, 0, sigma_true)
y <- as.vector(y)


# Combine data for Stan
stan_data <- list(N = N, T = T, S = S, Sigma = Sigma, C = C, y = y, p=p)

# Read in the Stan model from the file
model_code <- file.path(getwd(), "random_effect_model.stan")

fit <- stan(
  file = model_code,
  data = stan_data,
  chains = 2,
  iter = 2000,
  warmup = 500
)

# Check the results
print(fit)

# You can further analyze the results using diagnostics, posterior predictive checks, etc.
