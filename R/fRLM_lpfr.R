#' fRLM Function
#'
#' Wrapper function that estimates a Functional Regression Linear Model (fRLM) with a given set of predictors.
#'
#' @param data A tidy dataframe containing the data.
#' @param id A string representing the identifier variable in the dataframe.
#' @param time A string representing the time variable in the dataframe, at what time the outcome or exposures were measured.
#' @param exposures A character vector of strings representing the exposure variables in the dataframe. For now limited to length 1.
#' @param outcome A string representing the outcome variable in the dataframe.
#' @param family A string representing the family. Either "gaussian" or "binary".
#' @param bounded_exposures A boolean vector specifying whether the exposure should be bounded between 0 and 1. For now of length 1.
#' @param controls A character vector of strings representing the control variables in the dataframe. Default is NULL.
#' @param L An integer representing the number of basis functions. Default is 4.
#' @param K An integer representing the number of knots for each basis function. Default is K=5.
#' @param ... Additional arguments passed to the `rstan::stan` function.
#'
#' @return A list containing the fitted model and other relevant information.
#'
#' @details
#' The `outcome`and `controls` must be the same (repeated) for all time stamps for each unit. See the examples below.
#'
#'
#' @examples
#'
#' library(fRLM)
#' library(dplyr)
#'
#' data(toy)
#'
#' # Gaussian outcome
#' fit <- fRLM_lpfr(
#'   data = toy,
#'   id="id",  # string, name of the subject identifier
#'   exposures = "exposure", # string, the variable name of the exposures (for now, one)
#'   outcome="outcome", # string, the variable name of the outcome
#'   time = "age", # string, the variable name of the time at which measures were taken
#'   warmup= 500,
#'   iter=1000,
#'   chains=4
#' )
#'
#' plot(fit)
#' w = predict(fit)
#'
#'
#' # Binary outcome and bounded exposure (between 0 and 1)
#' toy <- toy %>% group_by(id) %>% mutate(outcome_binary = rbinom(1,1,1/(1+exp(-outcome)))) %>% ungroup()
#'
#' fit_binary <- fRLM_lpfr(
#'   data = toy,
#'   id="id",  # string, name of the subject identifier
#'   exposures = "exposure", # string, the variable name of the exposures (for now, one)
#'   outcome="outcome_binary", # string, the variable name of the outcome
#'   family="binary",
#'   time = "age",
#'   bounded_exposures = TRUE,
#'   warmup= 500,
#'   iter=1000,
#'   chains=4
#' )
#'
#' plot(fit_binary)
#' w_binary = predict(fit_binary)
#' @export

fRLM_lpfr <- function(data, id, time, exposures, outcome, family="gaussian", bounded_exposures=FALSE, controls=NULL, L=4, K=5, grid= seq(0,1,l=150), ...) {
  # Create the time scaler. This function scales the data to standardized scale (on [0,1]) and back to the original scale
  timeScaler <- fRLM::timeScaler_ff(data %>% pull(!!sym(time)), min=0.06, max=0.93)
  # Standardize time
  data <- data %>% mutate(!!time := timeScaler(data[[time]]))

  # Extract quantities
  # ------------------

  # Extract y
  y_with_id <- data %>% group_by(!!sym(id)) %>% summarise(!!outcome := mean(!!sym(outcome)))
  y <- y_with_id %>% dplyr::pull(!!sym(outcome))
  # Extract tobs (xobs, resp.) as a list of times for each unit (and corresponding exposure, resp.)
  data_split <- split(data, data[id])
  tobs <- lapply(data_split, function(di) di[time] %>% unlist %>% as.vector)
  xobs <- lapply(data_split, function(di) di[exposures] %>% unlist %>% as.vector)

  # Padding, adding observation at t=0 with values 0 TODO: improve the padding, or remove them from the likelihood
  tobs <- padding(tobs)
  xobs <- padding(xobs)


  stopifnot(all.equal(tobs$mask, xobs$mask))

  # Create the Nvec list TODO: use the mask instead
  Nvec <- apply(tobs$mask, 1, function(x)sum(x==0))
  Nmax <- max(Nvec)

  # Extract the covariates
  C <- data %>% group_by(!!sym(id)) %>% summarise(across(all_of(controls), ~mean(.))) %>% ungroup() %>%
    # select the controls
    select(all_of(controls)) %>%
    mutate(intercept=1, .before=1) %>% # add intercept (don't! already is in Julien's code!)
    as.matrix

  # Bases:
  phi <- getBasis(L, grid) # density splines TODO: write a single function get_splines that can also have density splines
  psi <- get_splines(grid, K)$psi
  phi_mat <- apply( tobs$mat, 1, function(tt) get_splines(tt, K )$psi, simplify = F)
  J <- t(psi) %*% phi / length(grid)

  # Stan
  # -------------------------

  if (family == "gaussian") {
    file_name <- "lpfr.stan"
  } else if (family == "binary" & bounded_exposures) {
    file_name <- "lpfr_binary.stan"
    y = as.integer(y)
  } else {
    stop("Combination of family and bounded_exposure not yet implemented.")
  }

  stan_file <- system.file("stan", file_name, package = "fRLM")

  # Creat the data for stan
  dat <- list(
    n = length(y),
    L = L,
    K = K,
    d = ncol(C),
    Nmax = Nmax,
    y = y,
    C = C,
    phi_mat = phi_mat,
    J = J,
    Nvec = Nvec,
    tobs = tobs$mat,
    xobs = xobs$mat )

  fit <- rstan::stan( file = stan_file, data = dat, ... )

  out <- c( fit = fit, rstan::extract(fit ), L = L )
  out$psi <- psi
  out$Xhat <- t( sapply( apply(out$xi, 2, function(x) x %*% t(psi), simplify = F), colMeans ) )
  out$grid <- timeScaler(grid, original_scale=TRUE)
  out$timeScaler <- timeScaler
  class(out) <- "funcRegBayes"
  return(out)
}

# padding with mask
# mask is 1 if padded, 0 othewerise
#' @export
padding <- function(list_of_vec) {
  # Maximum number of observations across all subjects
  Nmax <- lapply(list_of_vec, function(x) length(x)) %>% unlist %>% max
  # Function to pad a single vector and create its mask
  pad_and_mask <- function(x) {
    pad_length <- Nmax - length(x)
    padded <- c(x, rep(0, pad_length))
    mask <- c(rep(0, length(x)), rep(1, pad_length))
    list(padded = padded, mask = mask)
  }

  # Apply the function to each vector in the list
  results <- lapply(list_of_vec, pad_and_mask)

  # Combine the padded vectors and masks
  padded_matrix <- do.call(rbind, lapply(results, `[[`, "padded"))
  mask_matrix <- do.call(rbind, lapply(results, `[[`, "mask"))
  list(mat = padded_matrix, mask = mask_matrix)
}


#' Generate Spline Basis and Derivatives
#'
#' This function generates a spline basis and its derivatives using a specified grid, number of knots,
#' and degree. It uses the `orthogonalsplinebasis` package for creating and manipulating spline basis.
#'
#' @param grid A numeric vector that specifies the points at which the splines and their derivatives
#'             are to be evaluated.
#' @param K An integer specifying the total number of knots to be used in the spline basis.
#'          It is to be noted that the actual number of knots used inside the function will be `K - 2`.
#' @param range A numeric vector of length 2 specifying the range over which the knots are distributed.
#'              Defaults to c(0, 1).
#' @param degree An integer specifying the degree of the spline. Defaults to 3.
#'
#' @return A list containing the following elements:
#'         - `psi`: The evaluated spline basis at the specified grid points.
#'         - `basis`: The spline basis object as created by `SplineBasis`.
#'         - `d_basis`: The first derivative of the spline basis.
#'         - `d2psi`: The second derivative of the spline basis evaluated at the grid points.
#'
#' @examples
#' grid <- seq(0, 1, length.out = 100)
#' splines <- get_splines(grid, K = 5)
#' plot(grid, splines$psi[, 1], type = 'l')
#'
#' @importFrom orthogonalsplinebasis SplineBasis evaluate deriv OuterProdSecondDerivative
#' @export
get_splines <- function( grid, K, range = c(0,1), degree = 3 ){
  # Calculate the number of knots as K - 2
  nknots <- K - 2

  # Generate knots within the specified range
  knots <- orthogonalsplinebasis::expand.knots( seq( range[1], range[2], length.out = nknots ) )

  # Create spline basis with the specified order (degree + 1)
  basis <- orthogonalsplinebasis::SplineBasis(knots = knots, order = degree + 1 )

  # Compute the first derivative of the spline basis
  d_basis <- orthogonalsplinebasis::deriv( basis)

  # Evaluate the spline basis at the specified grid points
  psi <- orthogonalsplinebasis::evaluate( basis, grid )

  # Compute the second derivative of the spline basis
  d2psi <- orthogonalsplinebasis::OuterProdSecondDerivative(basis)

  # Return the spline basis and its derivatives
  return( list( psi = psi, basis = basis, d_basis = d_basis, d2psi = d2psi ) )
}
