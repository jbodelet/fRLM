#' @export
simulateFRLM <- function( sampleSize = 100, alpha = 0, delta = 3, omega = function(t){ 4 * t^3 }, errorStd = 1, gpStd = 1 ){
  # simulate a Functional Relevant Life course Model:
  # y_i = alpha + delta int_0^1 omega(t) X_i(t) dt + eps_i,  i=1,...,sampleSize
  errors <- rnorm( sampleSize, sd = errorStd )
  functionalCovariates <- replicate( sampleSize, generateGaussianProcess( sig2 = gpStd) ) 
  predictors <- sapply( functionalCovariates, function(X){ 
    Xomega_Product <- function(t) X(t) * omega(t)
    myIntegrate(Xomega_Product)
    })
  y <- alpha + delta * predictors + errors
  return( list( y = y, functionalCovariates = functionalCovariates, delta = delta, omega = omega, alpha = alpha, errorStd = errorStd ) )
}


myIntegrate <- function(fun, grid = seq(0,1, l = 200 ) ){
  mean( fun( grid ) )
}



