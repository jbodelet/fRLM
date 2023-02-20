#' Gaussian process simulation
#' 
#' Generate a Gaussian process on the interval [0,1].
#' @param mu Mean function.
#' @param kernelFunc Covariance kernel function (i.e. a positive function of two arguments). Default is the exponential kernel.
#' @param sig2 Variance of the gaussian process.
#' @param gridLength Number of points for simulations.
#' @return A function, i.e. a realization of a Gaussian process.
#' @export
#' @examples 
#'   gp <- generateGaussianProcess()
#'   tobs <- sort(runif(5))
#'   Xobs <- gp(tobs)
#'   gpHat <- gpFit(tobs, Xobs)
#'   plot(gpHat)
generateGaussianProcess <- function( mu = NULL, kernelFunc = NULL, sig2 = 1, gridLength = 100 ){
  if( is.null( mu ) ){ mu   <- zeroFunction }
  if( is.null(kernelFunc) ){ kernelFunc <- function(t,s) expoKernelMat(t,s, beta = 1, sigma2 = sig2 ) }
  grid <- seq( 0, 1, l = gridLength)
  VarMat <- kernelFunc(grid, grid)
  gp <- MASS::mvrnorm( mu = mu(grid), Sigma = VarMat )  # much faster!
  return( approxfun(grid, gp ) )
}




#' Gaussian process regression
#' 
#' Fit a Gaussian process regression with penalyzed likelihood.
#' @param tobs Observed time points.
#' @param Xobs Observed points of the gaussian process.
#' @param lambda penalty.
#' @param optimizationRange A vector containing the end-points of the interval to be searched.
#' @param power Power parameter for the exponential kernel.
#' @param nugget Nugget parameter for stable inversion of the correlation matrix.
#' @return  An list containing estimates
#' @export
#' @examples 
#'   gp <- generateGaussianProcess()
#'   tobs <- sort(runif(5))
#'   Xobs <- gp(tobs)
#'   gpHat <- gpFit(tobs, Xobs)
#'   plot(gpHat)
gpFit <- function( tobs, Xobs, lambda = 0.01, optimizationRange = c(-2, 4), power = 1.95, nugget = 0 ){
  loss <- Vectorize( function(beta){
    N <- length(tobs)
    corrMat <- getCorrMat( beta, tobs, nugget)
    muHat <- getMuHat(corrMat, tobs )
    sigmaHat <- getSigmaHat( Xobs, muHat, corrMat )
    log( det(corrMat) ) + N * log( N * sigmaHat )  
  })
  penalyzedLoss <- function(beta) loss(beta) + lambda * 10^beta
  # Optimization:
  opt <- optimize( penalyzedLoss, optimizationRange )
  beta <- opt$min
  corrMat <- getCorrMat( beta, tobs, nugget)
  mu <- getMuHat(corrMat, tobs)
  sig2 <- getSigmaHat( Xobs, mu, corrMat )
  # Output:
  out <- list( tobs = tobs, Xobs = Xobs, beta = beta, sig2 = sig2, mu = mu,
               nugget = nugget, power = power, lambda = lambda, 
               correlation_param = list(power = power ) 
  )
  class(out) <- "gpFit"
  return( out )
}




#==================
# METHODS:
#==================




#' @export
predict.gpFit <- function( object, tnew = object$tobs ){
  N <- length(object$Xobs)
  M <- length(tnew)
  tobs <- object$tobs
  corrMatObs <- getCorrMat(object$beta, tobs, object$nugget)
  mu <- getMuHat(corrMatObs, tobs)
  corrMatNewWithObs <- expoKernelMat( tnew, tobs, object$beta, 1 )
  Xnew <- rep(mu,M) + corrMatNewWithObs %*% solve(corrMatObs) %*% ( object$Xobs - rep(mu,N) )
  return( c(Xnew) )
}

#' @export
plot.gpFit <- function(object){
  grid <- seq(0,1, l = 150)
  pred <- predict(object, grid)
  kernel <- getCondKernelFunc(object)
  std <- sqrt( diag(kernel(grid, grid) ) )
  upperBand <- pred + 1.96 * std
  lowerBand <- pred - 1.96 * std
  ylim <- c( min(lowerBand), max(upperBand) )
  plot( pred ~ grid, type = "l", lwd = 2, ylim = ylim )
  points( object$Xobs ~ object$tobs, col = "red", lwd = 4 )
  lines( upperBand ~ grid, col = "lightblue")
  lines( lowerBand ~ grid, col = "lightblue")
  invisible(std)
}


#' @export
lines.gpFit <- function(object){
  grid <- seq(0,1, l = 150)
  pred <- predict( object, grid )
  lines( pred ~ grid, lwd = 2, type = "l", col= "blue" )
  points( object$Xobs ~ object$tobs, col = "red", lwd = 4 )
}

evalGP <- function(gp, l = 150 ){
  grid <- seq(0,1, l = l )
  return( list( X = gp(grid), tt = grid ) )
}


#=====================
# Internal functions:
#=====================


getCondMuFunc <- function( gpfit, gridLength = 150){
  grid <- seq(0,1, l = gridLength )
  pred <- predict( gpfit, grid )
}




getCondKernelFunc <- function( gpfit ){
  tobs <- gpfit$tobs
  kern <- function(t,s) expoKernelMat( t, s, gpfit$beta, sigma2 = gpfit$sig2 )
  invCovMat <- solve( expoKernelMat( tobs, tobs, gpfit$beta, sigma2 = gpfit$sig2 ) )
  return( function(t, s) pmax( kern(t, s) - kern(t, tobs) %*% invCovMat %*% kern(tobs, s ), 0 ) )  # sometimes, gets negative, when t is close to 0.
}


expoKernelMat <- function(xvec, yvec, beta, sigma2, power = 1.96 ){
  do.call( rbind, lapply( xvec, function(x) expoKernel( x, yvec, beta, sigma2, power ) ) )
}

expoKernel <- function( x, y, beta, sigma2 = 1, power = 1.96 ){ sigma2 * exp( -10^beta * abs( x - y )^power ) } 

zeroFunction <- Vectorize( function(x) 0  )

# corr is not kernel!!
getCorrMat  <- function(beta, tobs, nugget = 0 ) expoKernelMat( tobs, tobs, beta, 1 ) + nugget * diag(length(tobs))

getMuHat    <- function( corrMat, tobs ){ sum( solve( corrMat ) %*% tobs ) / sum( solve( corrMat ) ) }

getSigmaHat <- function( Xobs, muHat, corrMat ) as.numeric( t( Xobs - muHat ) %*% solve( corrMat ) %*% (Xobs - muHat ) / length(Xobs) )




#==================
# Tests:
#==================


if(0){
  # # Test for generateGaussianProcess and gpFit:
  gp <- generateGaussianProcess()
  tobs <- sort(runif(5))
  Xobs <- gp(tobs)
  object <- gpFit(tobs, Xobs)
  plot(object)
  
  # Redo it:
  # # Test for getCondKernelFunc:
  # grid <- seq(0,1, l = 150)
  # pred <- predict(gpHat, grid)
  # condMuFunc <- approxfun(grid, condMu )
  # condKernel <- getCondKernelFunc(gpHat)
  # simulatedGP <- replicate( 200, generateGaussianProcess( mu = condMuFunc, kernel = condKernel ) )
  # gpMat <- sapply( simulatedGP, function(x) evalGP( x )$X )
  # 
  # 
  # plot(gp)
  # grid <- seq(0,1, l=150)
  # apply( gpMat, 2, lines, x = grid, col ="gray" )
  # lines(gp(grid)~ grid, lwd = 2, col = "black")
  # lines(condMu(grid) ~ grid, lwd = 2, col = "blue")
  # points( Xobs ~ tobs, col = "red", lwd = 4 )
}



