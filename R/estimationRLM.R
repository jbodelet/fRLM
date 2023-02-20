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
  dat <- list(n = n, n_obs = n_obs, J = ncol(x), K = ncol(u), d = ncol(C), y = y[new_index], z_obs = z_obs[new_index][1:n_obs], 
              x = x[new_index, ], u = u[new_index, ], C = as.matrix(C[new_index,] ) )
  if(robust){
    fileName <- "discreteModel_robust_missing_income.stan"
  }else{
    fileName <- "discreteModel_missing_income.stan"
  }
  stanFile <- system.file("stan", fileName, package = "fRLM")
  fit <- rstan::stan( file = stanFile, data = dat, ...)
}


#' Multivariate RELEVANT LIFE COURSE MODEL
#' 
#' Estimate the multivariate version of Madathil's regression model:
#' y_ij = delta_j X_it * weights_tj + C_i' alpha_j + eps_ij
#' 
#' @param y  The output matrix.
#' @param x  A n times M matrix containing the n functions evaluated at points specified in grid.
#' @param Covariates A matrix of covariates.
#' @param robust If TRUE, use a t-distribution as a likelihood.
#' @param ... Optional parameters to pass to the stan function.
#' @return An list containing estimates for alpha, delta, beta, gamma, sigma2, and omega.
#' @examples 
#' n <- 100
#' p <- 4
#' J <- 5
#' x <- matrix( rnorm(n*J), ncol = J )
#' delta <- c(2, 1, rep(0, p-2) )
#' weights <- rbind( c(0.1, 0.3, 0.6), c(0.8, 0.2, 0), matrix(0, nrow = p-2, ncol = J) )
#' sigma <- c(1, 0.5, rep(1, p-2))
#' x <- matrix( rnorm(n*J), ncol = J )
#' eps <- matrix( rnorm(n * p), ncol = p ) %*% diag(sigma)
#' y <- x %*% t(weights) %*% diag(delta) + eps
#' fitMV_RLM(y, x, warmup = 2000, iter = 4000, chains = 4)
#' @export
fitMV_RLM <- function(y, x, covariates = NULL, robust = FALSE, ...){
  C <- cbind( rep(1, nrow(y) ), covariates )  # scalar predictors
  dat <- list(n = nrow(y), p = ncol(y), J = ncol(x), d = ncol(C), y = y, x = x, C = C)
  if(robust){
    fileName <- "MV_RLM_robust.stan"
  }else{
    fileName <- "MV_RLM.stan"
  }
  stanFile <- system.file("stan", fileName, package = "fRLM")
  fit <- rstan::stan( file = stanFile, data = dat, ...)
}


#' @export
fitMV_RLM_missingIncome <- function(y, x, z_obs, u, covariates = NULL, robust = FALSE, ...){
  n <- nrow(y)
  p <- ncol(y)
  n_obs <- sum( !is.na(z_obs) )
  C <- cbind( rep(1, n ), covariates )  # scalar predictors
  # data:
  obs_index <- (1:n)[!is.na( z_obs )]
  miss_index <- (1:n)[is.na( z_obs )]
  new_index <- c(obs_index, miss_index)
  dat <- list(n = n, p = p, n_obs = n_obs, J = ncol(x), K = ncol(u), 
              d = ncol(C), y = y[new_index,], z_obs = z_obs[new_index][1:n_obs], 
              x = x[new_index, ], u = u[new_index, ], C = as.matrix(C[new_index, 
              ]))
  if(robust){
    fileName <- "MV_RLM_missingIncome_robust.stan"
  }else{
    fileName <- "MV_RLM_missingIncome.stan"
  }
  stanFile <- system.file("stan", fileName, package = "fRLM")
  fit <- rstan::stan( file = stanFile, data = dat, ...)
}
