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
#' @param L  Number of B-splines basis to approximate omega.
#' @return An list containing estimates for alpha, delta, beta, gamma, sigma2, and omega.
#' @export
#' @examples 
#' library(fRLM)
#' library(dplyr)
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
#' # 3) Frequentist Estimation:
#' fit <- funcReg( y, condMu, grid, L = 5 )
#' bootfit <- bootFunReg( y, condMu, grid, L = 5 )
#' plot(bootfit)
funcReg <- function( y, X, grid = seq(0,1, l = ncol(X) ), covariates = NULL, L = 6 ){
  stopifnot( nrow(X) == length(y) )
  stopifnot( ncol(X) == length(grid) )
  C <- cbind( rep(1, length(y) ), covariates )  # scalar predictors
  basis <- getBasis( L, grid )
  Int_XtimesBasis <- X %*% basis / ncol(X)
  # Algorithms:
  coeff <- constrainedOLS( y, cbind( Int_XtimesBasis, C ), constrainedId = 1:L )
  # Output:
  beta2 <- coeff[1:L]
  alpha <- coeff[-(1:L) ]
  predictor <- Int_XtimesBasis %*% beta2  + C %*% alpha
  sigma2 <- mean( ( y - predictor )^2 )
  out <- list( delta = sum(beta2), functionalCoeff = beta2 / sum(beta2), covariatesCoeff = alpha, sigma2 = sigma2,
               y = y, X = X, C = C, L = L )
  class(out) <- "funcReg"
  return( out )
}


getBasis <- function( L, grid ){
  knots <- seq(0, 1, length.out = L - 2 )
  return( t( bsplinePsd::dbspline( grid, knots = knots ) ) )  # basis for omega
}

#===================
# METHODS:
#==================

#' @export
predict.funcReg <- function( object, newdata = seq(0,1, l = 150), type = "omega"){
  basis <- getBasis( object$L, grid = newdata )
  return( basis %*% object$functionalCoeff )
}


#' @export
plot.funcReg <- function(object){
  grid <- seq(0,1, l = 150)
  plot( predict(object ) ~ grid, lwd = 2, type = "l" )
}




#===========================
# BOOTSTRAP:
#===========================




#' BOOTSRTAP FOR THE FUNCTIONAL LIFECOURSE MODEL 
#' 
#' Bootsrtap of the functional regression,
#' # y_i = delta int X_i(t) * omega(t) dt + C_i' alpha + eps_i
#' where omega(t) is a weight function, omega(t) > 0 and in omega(t) dt = 1.
#' 
#' @param y  The output vector.
#' @param X  A n times M matrix containing the n functions evaluated at points specified in grid.
#' @param grid A vector of the points at which the functions are evaluated.
#' @param covariates A matrix of covariates.
#' @param nBoot Number of bootstrap samples.
#' @param L  Number of B-splines basis to approximate omega.
#' @return An list containing estimates for alpha, delta, beta, gamma, sigma2, and omega.
#' @export
bootFunReg <- function( y, X, grid, covariates = NULL, L = 6, nBoot = 100 ){
  statistic <- function( y, indices ){
    out <- funcReg( y[indices], X[indices, ], grid, L = L )
    c( out$delta, out$functionalCoeff, out$covariatesCoeff )
  }
  bootSamples <- boot::boot( data = y, statistic, R = nBoot )$t
  bootSamples_mean <- colMeans( bootSamples )
  bootSamples_var  <- var( bootSamples )
  nbCov <- length( bootSamples_mean ) - 2 - L
  names( bootSamples_mean ) <- rownames( bootSamples_var ) <- colnames(bootSamples_var) <- 
    c( "delta", paste0("beta", 1:L), paste0("alpha", 0:nbCov) )
  out <- list( bootSamples_mean = bootSamples_mean, bootSamples_var = bootSamples_var, bootSamples = bootSamples, L = L )
  class(out) <- "bootFunReg"
  return( out )
}


#' @export
predict.bootFunReg <- function( object, newdata = seq(0,1, l = 150), type = "omega" ){
  L <- object$L
  basis <- getBasis( L, grid = newdata )
  beta_index <- 1:L + 1
  omega_mean <- basis %*% object$bootSamples_mean[beta_index ]
  omega_sd <- sqrt( diag( basis %*% object$bootSamples_var[beta_index,  beta_index ] %*% t( basis ) ) )
  return( list( omega_mean = omega_mean, omega_sd = omega_sd, grid = grid ) )
}


#' @export
plot.bootFunReg <- function(object){
  grid <- seq(0,1, l = 150)
  pred <- predict(object, grid)
  upperBand <- pred$omega_mean + 1.96 * pred$omega_sd
  lowerBand <- pred$omega_mean - 1.96 * pred$omega_sd
  ylim <- c( min(lowerBand), max(upperBand) )
  plot( pred$omega_mean ~ grid, lwd = 2, type = "l", col = "blue", ylim = ylim )
  lines( upperBand ~ grid, col = "lightblue" )
  lines( lowerBand ~ grid, col = "lightblue" )
}





#===========================
# CONSTRAINED ESTIMATION:
#===========================


constrainedOLS <- function(y, X, constrainedId ){
  # y = X %*% beta + eps, s.t. sign(beta_j) > 0 or < 0, for j in constrainedID
  Q <- t( X ) %*% X
  a <- - 2 * t( X) %*% y
  return( constrainedOptimFindSign( Q, a, constrainedId = constrainedId ) )
}


constrainedOptimFindSign <- function(Q, a, constrainedId = 1:nrow(Q) ){
  # min_x x'Qx +a'x + d
  # s.t.  sign( x[j]) =  sign( x[k]), for all j,k in constrainedId
  # it solves for negative, then for positive, 
  # and return the best solution of the two, by checking the reached minimum.
  out_lb <- constrainedOptim( Q, a, constrainedId, signIsPositive = TRUE )
  out_ub <- constrainedOptim( Q, a, constrainedId, signIsPositive = FALSE )
  if( out_lb$mse < out_ub$mse ){
    out <- out_lb
  }else{
    out <- out_ub
  }
  out$x <- round( as.numeric(out$x), 3 )  
}



constrainedOptim <- function( Q, a, constrainedId, signIsPositive = TRUE ){
  # solve:  x' Q x + x' a   
  # with constraint on the signs of x: 
  # positive (signIsPositive = FALSE) or negative (signIsPositive = TRUE)
  fun2min <- optiSolve::quadfun( Q = Q, a = a )
  bounds <- rep(0, length(constrainedId) )
  ub <- lb <- NULL
  if(signIsPositive){
    lb <- optiSolve::lbcon( bounds, id = constrainedId )
  }else{
    ub <- optiSolve::ubcon( bounds, id = constrainedId )
  }
  cop <- optiSolve::cop( f = fun2min, lb = lb, ub = ub )
  xhat <- optiSolve::solvecop( cop, quiet = TRUE )$x
  mse <- c( t(xhat) %*% Q %*% xhat + t( xhat ) %*% a )
  return( list(xhat = xhat, mse = mse ) )
}




#===========================
# TESTS:
#===========================



if(0){
  library(fRLM)
  library(dplyr)
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
  
  
  # 3) Frequentist Estimation:
  fit <- funcReg( y, condMu, grid, L = 5 )
  bootfit <- bootFunReg( y, condMu, grid, L = 5 )
  plot(bootfit)
}
