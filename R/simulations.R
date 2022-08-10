#' @export
simulateFuncReg <- function( n = 100, alpha = 0, delta = 3, omega = NULL, errorStd = 1, gpStd = 1 ){
  # simulate a functional regression model:
  # y_i = alpha + delta int_0^1 omega(t) X_i(t) dt + eps_i,  i=1,...,n
  if(is.null(omega)){
    omega <- function(t){ 4 * t^3 }    
  }
  errors <- rnorm(n, sd = errorStd )
  gpList <- replicate( n, generateGaussianProcess( sig2 = gpStd) ) 
  functionalPredictors <- lapply( gpList, function(gp){ function(x) gp(x) * omega(x) } )
  scalarPredictors <- sapply( functionalPredictors, function(fun){
    integrate( fun, lower = 0, upper = 1, subdivisions = 200, rel.tol = .Machine$double.eps^0.15 )$value 
  })
  y <- alpha + delta * scalarPredictors + errors
  return( list( y = y, gpList = gpList, delta = delta, omega = omega, alpha = alpha, errorStd = errorStd ) )
}

#' @export
getTimePoints <- function( n, N = 5, method = "sumOfUniforms" ){
  if(method == "random"){
    tobs <- replicate( n, cumsum( rdirichlet( N + 1, 2 ) )[1:N], simplify = FALSE )  # dirichlets
  }
  if(method == "sortedUniforms" ){
    tobs <- replicate( n, sort( runif(N) ), simplify = FALSE )
  }
  if(method == "rep10"){ # repeat the same 10 observed timePoint vectors, to speed up computation in the functional PCA.
    tobs <- rep( replicate( 10, cumsum( rdirichlet( N + 1, 2 ) )[1:N], simplify = FALSE ), n / 10 )
  }
  if(method == "rep10Random"){ # same but index is randomized
    tobs_uniq <- replicate( 10, cumsum( rdirichlet( N + 1, 2 ) )[1:N], simplify = FALSE )
    tobs <- rep( tobs_uniq, n / 10 )
    index <- c( 1:10, sample(1:10, size = n - 10, replace = TRUE ) )
    tobs <- tobs[index]
  }
  if(method == "sumOfUniforms"){
    tobs <- replicate(n,{
      u <- runif(N+1);
      cumsum(u[-1]) / sum(u)
    }, simplify = FALSE )
  }
  tobs <- lapply(tobs, round, digits = 2 )
  return( tobs )
}



#' @export
rdirichlet <- function( n = 2, a = 1 ){ # a can be either a scalar or a vector of length n.
  if( length(a) == 1 ){ a <- rep(a, n ) }
  x <- sapply( a, function(s) rgamma(1, shape = s ) )
  return( x / sum(x) )
}



# 
# library(fRLM)
# 
# 
# # 1) Simulations:
# n <- 200
# sim <- simulateFuncReg(n = n, errorStd = 1, gpStd = 1 )
# timePoints <- getTimePoints( n, N = 5, method = "rep10Random" )
# tobs <- timePoints$tobs
# XobsList <- lapply( 1:n, function(i) sim$gpList[[i]]( tobs[[i]] ) )
# 
# # 2) Fit gaussian processes:
# gpfitList <- lapply( 1:n, function(i) gpFit( tobs[[i]], XobsList[[i]] ) )
# grid <- seq(0,1, l = 150 )
# condMu <- t( sapply( gpfitList, predict, tnew = grid ) )


