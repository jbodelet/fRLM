#' DISCRETE RELEVANT LIFE COURSE MODEL
#' 
#' Estimate Madathil's regression model
#' y_i = delta X_ij * weights_j + C_i' alpha + eps_i
#' 
#' @param y  The output vector.
#' @param x  A n times M matrix containing the n functions evaluated at points specified in grid.
#' @param Covariates A matrix of covariates.
#' @param robust If TRUE, use a t-distribution as a likelihood.
#' @param ... Optional parameters to pass to the stan function.
#' @return An list containing estimates for alpha, delta, beta, gamma, sigma2, and omega.
#' @examples 
#' n <- 100
#' J <- 5
#' x <- matrix( rnorm(n*J), ncol = J )
#' eps <- rnorm(n)
#' delta <- 2
#' weights <- runif(J); 
#' weights <- weights / sum(weights)
#' y <- as.numeric( delta * ( x %*% weights) + eps )
#' fitRLM(y, x, warmup = 2000, iter = 4000, chains = 4)
#' @export
fitRLM <- function(y, x, covariates = NULL, robust = FALSE, ...){
  C <- cbind( rep(1, length(y) ), covariates )  # scalar predictors
  dat <- list( n = length(y), J = ncol(x), d = ncol(C), y = y, x = x, C = C )
  if(robust){
    fileName <- "discreteModel_robust.stan"
  }else{
    fileName <- "discreteModel.stan"
  }
  stanFile <- system.file("stan", fileName, package = "fRLM")
  fit <- rstan::stan( file = stanFile, data = dat, ...)
}


#' @export
fitRLM_missingIncome <- function(y, x, z_obs, u, covariates = NULL, robust = FALSE, ...){
  n <- length(y)
  n_obs <- sum( !is.na(z_obs) )
  C <- cbind( rep(1, length(y) ), covariates )  # scalar predictors
  # data:
  obs_index <- (1:n)[!is.na( z_obs )]
  miss_index <- (1:n)[is.na( z_obs )]
  new_index <- c(obs_index, miss_index)
  dat <- list(n = n, n_obs = n_obs, J = J, K = K, d = d, y = y[new_index], z_obs = z_obs[new_index][1:n_obs], 
              x = x[new_index, ], u = u[new_index, ], C = as.matrix(C[new_index,] ) )
  if(robust){
    fileName <- "discreteModel_robust_missing_income.stan"
  }else{
    fileName <- "discreteModel_missing_income.stan"
  }
  stanFile <- system.file("stan", fileName, package = "fRLM")
  fit <- rstan::stan( file = stanFile, data = dat, ...)
}



