#' Goldsmith method
#' @param y: a vector of numeric values, the outcome, with vector elements corresponding to the units.
#' @param tobs: a list of time of observations, with list elements corresponding to the units
#' @param xobs: a list of exposures corresponding to the time observeations, with list elements corresponding to the units
#' @param covariates: a matrix of covariates, with rows corresponding to the units.
lpfr_stan <- function(y, tobs, xobs, L = 4, K = 5, covariates = NULL,
                           grid = seq(0,1, l = 150 ), compiled_file = NULL, ...){
  stopifnot(  K >= L ) # condition should be met (see 2010 Goldsmith paper)
  padding <- function(list_of_vec) t( sapply( list_of_vec, function(x) c( x, rep(0, Nmax - length(x) ) ) ) )
  Nvec <- sapply(tobs, length)
  Nmax <- max(Nvec)
  tobs_mat <- padding(tobs)
  xobs_mat <- padding(xobs)
  # bases:
  phi <- getDensitySplines(L, grid)
  psi <- get_splines(grid, K)$psi
  phi_mat <- apply( tobs_mat, 1, function(tt) get_splines(tt, K )$psi, simplify = F)
  J <- t(psi) %*% phi / length(grid)
  # fit:
  C <- cbind( rep(1, length(y) ), covariates )  # scalar predictors
  dat <- list( n = length(y), L = L, K = K, d = ncol(C), Nmax = Nmax, y = y, C = C, phi_mat = phi_mat,
               J = J, Nvec = Nvec, tobs = tobs_mat, xobs = xobs_mat )

  fit <- rstan::stan( file = compiled_file, data = dat, ... )
  fit_list <- list(
    fit = rstan::extract(fit),
    fit_object = fit,
    psi = psi,
    Xhat = t( sapply( apply(out$xi, 2, function(x) x %*% t(psi), simplify = F), colMeans ) )
  )
  return(fit)
}



get_splines <- function( grid, K, range = c(0,1), degree = 3 ){
  nknots <- K - 2
  knots <- orthogonalsplinebasis::expand.knots( seq( range[1], range[2], length.out = nknots ) )
  basis  <-  orthogonalsplinebasis::SplineBasis(knots = knots, order = degree + 1 )
  d_basis <- orthogonalsplinebasis::deriv( basis)
  psi <- orthogonalsplinebasis::evaluate( basis, grid )
  d2psi <- orthogonalsplinebasis::OuterProdSecondDerivative(basis)
  return( list( psi = psi, basis = basis, d_basis = d_basis, d2psi = d2psi ) )
}

