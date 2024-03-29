# Goldsmith method
#'@export
lpfr <- function(y, tobs, xobs, L = 4, K = 5, covariates = NULL,
                           grid = seq(0,1, l = 150 ), compiled_file = NULL, ...){
  stopifnot(  K >= L ) # condition should be met (see 2010 Goldsmith paper)

  Nvec <- sapply(tobs, length)
  Nmax <- max(Nvec)

  tobs_padded <- padding(tobs)
  tobs_mask <- tobs_padded$mask
  tobs_padded <- tobs_padded$padded

  xobs_padded <-
  xobs_mat <- padding(xobs)
  # bases:
  phi <- getBasis(L, grid) # density splines
  psi <- get_splines(grid, K)$psi
  phi_mat <- apply( tobs_mat, 1, function(tt) get_splines(tt, K )$psi, simplify = F)
  J <- t(psi) %*% phi / length(grid)
  # fit:
  C <- cbind( rep(1, length(y) ), covariates )  # scalar predictors
  dat <- list( n = length(y), L = L, K = K, d = ncol(C), Nmax = Nmax, y = y, C = C, phi_mat = phi_mat,
               J = J, Nvec = Nvec, tobs = tobs_mat, xobs = xobs_mat )
  if(is.null(compiled_file)){
    fit <- rstan::stan( file = "./stan/lpfr.stan", data = dat, ... )
  }else{
    fit <- rstan::sampling(compiled_file, data = dat, ...)
  }
  out <- c( fit = fit, rstan::extract(fit ), L = L )
  out$psi <- psi
  out$Xhat <- t( sapply( apply(out$xi, 2, function(x) x %*% t(psi), simplify = F), colMeans ) )
  class(out) <- "funcRegBayes"
  return(out)
}


