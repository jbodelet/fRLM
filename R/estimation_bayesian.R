#' FUNCTIONAL LIFECOURSE MODEL 
#' 
#' Estimate a functional regression,
#' # y_i = delta int X_i(t) * omega(t) dt + C_i' alpha + eps_i
#' where omega(t) is a weight function, omega(t) > 0 and in omega(t) dt = 1.
#' 
#' @param y  The output vector.
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
#' # parallel computing for the bayesian fit:
#' # library("rstan") # observe startup messages
#' # options(mc.cores = parallel::detectCores())
#' # rstan_options(auto_write = TRUE)
#' 
#' # 1) Setup data:
#' head(toy) # toy example dataset
#' toy2 <- toy
#' toy2[, 4] <- round( ( toy[, 4] - 10 ) / 30, 2 ) # i. rescale age:
#' y <- toy2 %>% group_by(id) %>% summarise(outcome = mean(outcome) ) %>% pull(outcome) # extract outcome:
#' t_obs <- toy2 %>% group_by(id) %>% summarise(age = list(age)) %>%  pull(age) # extract age
#' exposure <- toy2 %>% group_by(id) %>% summarise(exposure = list(exposure)) %>%  pull(exposure) # extract exposure
#' 
#' # 2) Fit gaussian processes:
#' gpfitList <- lapply( 1:length(y), function(i) gpFit( t_obs[[i]], exposure[[i]] ) )
#' grid <- seq(0,1, l = 150 )
#' condMu <- t( sapply( gpfitList, predict, tnew = grid ) )
#' 
#' # 3) Bayesian Estimation:
#' fitBayes <- funcReg_bayes( y, condMu, warmup = 500, iter = 1000, chains = 2 )
#' plot(fitBayes)
funcReg_bayes <- function( y, X, L = 5, covariates = NULL, alpha_par = 1, beta_par = 1, grid = seq(0,1, l = ncol(X) ), ... ){
  stopifnot( nrow(X) == length(y) )
  C <- cbind( rep(1, length(y) ), covariates )  # scalar predictors
  basis <- getBasis( L, grid )
  Int_XtimesBasis <- X %*% basis / ncol(X)
  # Stan:
  dat <- list( n = length(y), L = L, d = ncol(C), y = y, eta = Int_XtimesBasis, alpha_par = alpha_par, beta_par = beta_par, C = C )
  fileName <- "plugin.stan"
  stanFile <- system.file("stan", fileName, package = "fRLM")
  fit <- rstan::stan( file = stanFile, data = dat, ... )
  # output:
  alpha <- rstan::extract(fit, 'alpha')[[1]]
  beta  <- rstan::extract(fit, 'beta')[[1]]
  delta <- rstan::extract(fit, 'delta')[[1]]
  sigma <- rstan::extract(fit, 'sigma')[[1]]
  out <- list( fit = fit, delta = delta, alpha = alpha, beta = beta, sigma = sigma, basis = basis, L = L )
  class(out) <- "funcRegBayes"
  return(out)
}



#' @export
predict.funcRegBayes <- function( object, newdata = seq(0,1, l = 150), returnALL = FALSE ){
  L <- object$L
  basis <- getBasis( L, grid = newdata )
  omega_all <- basis %*% t(object$beta)
  omega <- rowMeans(omega_all)
  omega_ci <- apply( omega_all, 1, quantile, probs = c(0.025, 0.975) )
  out <- list( omega = omega, omega_ci = omega_ci, newdata = newdata )
  if(returnALL){
    out$omega_all <- omega_all    
  }
  return( out )
}


#' @export
plot.funcRegBayes <- function(object, ylim = range( pred$omega_ci ), ...){
  grid <- seq(0,1, l = 150)
  pred <- predict(object)
  plot( pred$omega ~ grid, lwd = 2, type = "l", col = "blue", ylim = ylim, ... )
  lines( pred$omega_ci[1, ] ~ grid, col = "lightblue" )
  lines( pred$omega_ci[2, ] ~ grid, col = "lightblue" )
}

getBasis <- function( L, grid ){
  knots <- seq(0, 1, length.out = L - 2 )
  return( t( bsplinePsd::dbspline( grid, knots = knots ) ) )  # basis for omega
}





#===========================
# TESTS:
#===========================

if(0){
  library(fRLM)
  library(dplyr)
  # parallel computing for the bayesian fit:
  # library("rstan") # observe startup messages
  # options(mc.cores = parallel::detectCores())
  # rstan_options(auto_write = TRUE)
  
  # 1) Setup data:
  head(toy) # toy example dataset
  toy2 <- toy
  toy2[, 4] <- round( ( toy[, 4] - 10 ) / 30, 2 ) # i. rescale age:
  y <- toy2 %>% group_by(id) %>% summarise(outcome = mean(outcome) ) %>% pull(outcome) # extract outcome:
  t_obs <- toy2 %>% group_by(id) %>% summarise(age = list(age)) %>%  pull(age) # extract age
  exposure <- toy2 %>% group_by(id) %>% summarise(exposure = list(exposure)) %>%  pull(exposure) # extract exposure
  
  
  # 2) Fit gaussian processes:
  gpfitList <- lapply( 1:length(y), function(i) gpFit( t_obs[[i]], exposure[[i]] ) )
  grid <- seq(0,1, l = 150 )
  condMu <- t( sapply( gpfitList, predict, tnew = grid ) )
  
  # 3) Bayesian Estimation:
  fitBayes <- funcReg_bayes( y, condMu, warmup = 500, iter = 1000, chains = 2 )
  plot(fitBayes)
}
