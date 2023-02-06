#' @examples 
simulateFRLM_binary <- function (n = 100, alpha = 0, delta = 3, omega = function(t) {4 * t^3}, gpStd = 1){
  functionalCovariates <- replicate(n, generateGaussianProcess(sig2 = gpStd))
  predictors <- sapply(functionalCovariates, function(X) {
    Xomega_Product <- function(t) X(t) * omega(t)
    myIntegrate(Xomega_Product)
  })
  eta <- alpha + delta * predictors
  prob <- 1/( 1+ exp(-eta) )
  y <- rbinom(n = n, size = 1, prob = prob )
  return(list(y = y, functionalCovariates = functionalCovariates, 
              delta = delta, omega = omega, alpha = alpha, errorStd = errorStd))
}

#' LOGISTIC FUNCTIONAL LIFECOURSE MODEL 
#' 
#' Estimate a logistic functional regression.
#' # y_i ~ Bernouilli(mu_i) 
#' # g(mu_i) = delta int X_i(t) * omega(t) dt + C_i' alpha
#' where omega(t) is a weight function, omega(t) > 0 and in omega(t) dt = 1.
#' 
#' @param y  The output vector (should be 0 and 1 values).
#' @param X  A n times M matrix containing the n functions evaluated at points specified in grid.
#' @param grid A vector of the points at which the functions are evaluated.
#' @param covariates A matrix of covariates.
#' @param beta_par Hyperparameter of the dirichlet prior of beta.
#' @param alpha_par Variance yperparameter of the prior of alpha.
#' @param L  Number of B-splines basis to approximate omega.
#' @param ... Optional parameters to pass to the stan function.
#' @return An list containing estimates for alpha, delta, beta, gamma, and sigma2.
#' @export
#' @examples 
#' library(fRLM)
#' library(dplyr)
#' library(devtools)
#' library(rstan)
#' # parallel computing for the bayesian fit:
#' # library("rstan") # observe startup messages
#' # options(mc.cores = parallel::detectCores())
#' # rstan_options(auto_write = TRUE)
#' n <- 1000
#' L <- 5
#' Nrange <- 3:5
#' Nvec <- sample( Nrange, n, replace = TRUE )
#' # Parameters:
#' omega <- function(t) cos(2 * pi * t) + 1
#' alpha <- 0.5
#' delta <- 1
#' grid <- seq(0,1, l = 150 )
#' # 1) Simulations:
#' sim <- simulateFRLM_binary( n, alpha = alpha, delta = delta, omega = omega )
#' tobs <- lapply( Nvec, function(x) sort( runif(x) ) )
#' Xobs <- lapply( 1:n, function(i) round( sim$functionalCovariates[[i]]( tobs[[i]] ), 3 ) )
#' 
#' # 2) Fit gaussian processes:
#' gpfitList <- lapply( 1:n, function(i) gpFit( tobs[[i]], Xobs[[i]], nugget = 0.001 ) )
#' pred <- t( sapply( gpfitList, predict, tnew = grid ) )
#' plot(gpfitList[[5]])
#' 
#' # 3) fit:
#' fit <- funcReg_bayes_binary(sim$y, pred, warmup = 2000, iter = 4000 )
#' plot(fit)
#' lines(omega(grid) ~ grid)
funcReg_bayes_binary <- function (y, X, L = 5, covariates = NULL, alpha_par = 1, beta_par = 1, grid = seq(0, 1, l = ncol(X)), ...){
  stopifnot(nrow(X) == length(y))
  C <- cbind(rep(1, length(y)), covariates)
  basis <- getBasis(L, grid)
  Int_XtimesBasis <- X %*% basis/ncol(X)
  dat <- list(n = length(y), L = L, d = ncol(C), y = y, eta = Int_XtimesBasis, 
              alpha_par = alpha_par, beta_par = beta_par, C = C)
  fileName <- "binary_fRLM.stan"
  # stanFile <- system.file("stan", fileName, package = "fRLM")
  stanFile <- "./binary_fRLM.stan"
  fit <- rstan::stan(file = stanFile, data = dat, ...)
  alpha <- rstan::extract(fit, "alpha")[[1]]
  beta <- rstan::extract(fit, "beta")[[1]]
  delta <- rstan::extract(fit, "delta")[[1]]
  out <- list(fit = fit, delta = delta, alpha = alpha, beta = beta, basis = basis, L = L)
  class(out) <- "funcRegBayes"
  return(out)
}

#' @export
simulateFRLM_tv_binary <- function (n = 100, alpha = 0, delta = 3, omega = function(t) {4 * t^3}, gpStd = 1 ){
  functionalCovariates <- replicate(n, generateGaussianProcess(sig2 = gpStd))
  TimeOfMeasurement <- runif(n, 0.5, 1 )
  predictors <- sapply( 1:n, function(i) {
    Xomega_Product <- function(t) functionalCovariates[[i]](t) * omega(t)
    myIntegrate2(Xomega_Product, to = TimeOfMeasurement[i] )
  })
  eta <- alpha + delta * predictors
  prob <- 1/( 1+ exp(-eta) )
  y <- rbinom(n = n, size = 1, prob = prob )
  return(list(y = y, functionalCovariates = functionalCovariates, 
              delta = delta, omega = omega, alpha = alpha, TimeOfMeasurement = TimeOfMeasurement ) )
}



myIntegrate2 <- function(fun, from = 0, to = 1 ){
  grid = seq( from, to, l = 200 )
  mean( fun( grid ) )
}


#' LOGISTIC FUNCTIONAL LIFECOURSE MODEL with time varying outcome
#' 
#' Estimate a logistic functional regression where final outcomes are not measured at the same times.
#' # y_i ~ Bernouilli(mu_i) 
#' # g(mu_i) = delta int_0^T_i X_i(t) * omega(t) dt + C_i' alpha
#' where T_i is the time point at which the outcome is measured.
#' 
#' @param y  The output vector (should be 0 and 1 values).
#' @param X  A n times M matrix containing the n functionsd evaluated at points specified in grid.
#' @param TimeOfMeasurement A vector of the same length as y, indicated the time points at which the outcomes are measured.
#' @param L  Number of B-splines basis to approximate omega.
#' @param grid A vector of the points at which the functions are evaluated.
#' @param covariates A matrix of covariates.
#' @param beta_par Hyperparameter of the dirichlet prior of beta.
#' @param alpha_par Variance yperparameter of the prior of alpha.
#' @param ... Optional parameters to pass to the stan function.
#' @return An list containing estimates for alpha, delta, beta, gamma, and sigma2.
#' @export
#' @examples 
#' @examples 
#' library(fRLM)
#' library(dplyr)
#' library(devtools)
#' library(rstan)
#' # parallel computing for the bayesian fit:
#' # library("rstan") # observe startup messages
#' # options(mc.cores = parallel::detectCores())
#' # rstan_options(auto_write = TRUE)
#' n <- 5000
#' L <- 7
#' Nrange <- 8:10
#' Nvec <- sample( Nrange, n, replace = TRUE )
#' 
#' # Parameters:
#' omega <- function(t) cos(2*pi*t) + 1
#' alpha <- -2
#' delta <- 1
#' grid <- seq(0,1, l = 150 )
#' 
#' # 1) Simulations:
#' sim <- simulateFRLM_tv_binary( n, alpha = alpha, delta = delta, omega = omega )
#' tobs <- lapply( 1:n, function(i){
#'   u <- runif( Nvec[i], min = 0.3, max = 0.7)
#'   cumsum(u) / sum(u) * sim$lastExpo[i]
#' })
#' # tobs <- c( sort( runif( Nvec[i] - 1, max = sim$lastExpo[i] ) ), sim$lastExpo[i] ) 
#' Xobs <- lapply( 1:n, function(i) round( sim$functionalCovariates[[i]]( tobs[[i]] ), 3 ) )
#' 
#' # 2) Fit gaussian processes:
#' gpfitList <- lapply( 1:n, function(i) gpFit( tobs[[i]], Xobs[[i]], nugget = 0.001 ) )
#' pred <- t( sapply( gpfitList, predict, tnew = grid ) )
#' plot(gpfitList[[ sample(1:n, 1)]])
#' 
#' # 3) fit:
#' # fit <- funcReg_bayes_binary(sim$y, pred, warmup = 2000, iter = 4000 )
#' # plot(fit)
#' # lines(omega(grid) ~ grid)
#' 
#' fit2 <- funcReg_bayes_tv_binary(sim$y, pred, sim$lastExpo, warmup = 2000, iter = 4000 )
#' plot(fit2)
#' lines(omega(grid) ~ grid)
#' 
funcReg_bayes_tv_binary <- function (y, X, TimeOfMeasurement, L = 5, covariates = NULL, alpha_par = 1, beta_par = 1, grid = seq(0, 1, l = ncol(X)), ...){
  stopifnot(nrow(X) == length(y))
  C <- cbind(rep(1, length(y)), covariates)
  basis <- getBasis(L, grid)
  # Int_XtimesBasis <- X[] %*% basis/ncol(X)
  Int_XtimesBasis <- t( sapply(1:n, function(i){
    integration_index <- which( grid < TimeOfMeasurement[i] )
    X[i, integration_index ] %*% basis[integration_index, ] / length(integration_index)
  }) )
  dat <- list(n = length(y), L = L, d = ncol(C), y = y, eta = Int_XtimesBasis, 
              alpha_par = alpha_par, beta_par = beta_par, C = C)
  fileName <- "binary_fRLM.stan"
  stanFile <- system.file("stan", fileName, package = "fRLM")
  fit <- rstan::stan(file = stanFile, data = dat, ...)
  alpha <- rstan::extract(fit, "alpha")[[1]]
  beta <- rstan::extract(fit, "beta")[[1]]
  delta <- rstan::extract(fit, "delta")[[1]]
  out <- list(fit = fit, delta = delta, alpha = alpha, beta = beta, basis = basis, L = L)
  class(out) <- "funcRegBayes"
  return(out)
}







