set.seed(123)
n <- 100
x <- rnorm(n)
y <- 2 + 3*x + rnorm(n)

library(rstan)

# compile model
mod <- stan_model("development/stan_test.stan")

# fit model
fit <- sampling(mod, data = list(N = n, x = x, y = y))

# print results
print(fit)
