#' fRLM Function
#'
#' Wrapper function that estimates a Functional Regression Linear Model (fRLM) with a given set of predictors.
#'
#' @param data A tidy dataframe containing the data.
#' @param id A string representing the identifier variable in the dataframe.
#' @param time A string representing the time variable in the dataframe.
#' @param exposures A character vector of strings representing the exposure variables in the dataframe.
#' @param outcome A string representing the outcome variable in the dataframe.
#' @param controls A character vector of strings representing the control variables in the dataframe. Default is NULL.
#' @param grouping A data-frame with the same variable name for respondent identification, and as many columns as there are exposures to the random effect, in chronological order. It is understood that the time-increment is constant between exposures.
#' @param L An integer representing the number of basis functions. Default is 5.
#' @param alpha_par A numeric value for the alpha hyperparameter. Default is 1.
#' @param beta_par A numeric value for the beta hyperparameter. Default is 1.
#' @param ... Additional arguments passed to the `rstan::stan` function.
#'
#' @return A list containing the fitted model and other relevant information.
#'
#' @details
#' The function begins by scaling the time variable using a time scaler function created with `timeScaler_ff`.
#' Next, the function defines a grid of 150 points between 0 and 1 and fits a Gaussian process for each exposure
#' in the `exposures` vector using the `gpFit` function. The `gpFit` function is applied to each subject
#' and the predicted values are computed on the grid for each subject.
#' The resulting predicted values are then multiplied by a set of basis functions computed on the grid using the `getBasis` function,
#' and averaged across subjects.
#' The `data` is then aggregated by `id` and the means of the `outcome` and `controls` variables are computed for each subject.
#' The `data` is then passed to the `rstan::stan` function to fit a Bayesian regression model.
#' The fitted model and other relevant information is then returned as a list.
#'
#'
#' @examples
#'
#' # Minimal example usage of the fRLM function
#' library(fRLM)
#' data(toy)
#' output <- fRLM(data=toy,
#'   id = "id",
#'   time="age",
#'   exposures="exposure",
#'   outcome="outcome",
#'   warmup = 500, iter = 1000, chains = 2) # this is passed to stan
#'
#' # Example with more exposures and random effects
#'
#'
#' # Example usage of the fRLM function
#'
#' library(fRLM)
#' set.seed(1234)
#' data(toy)
#'
#' toy2 <- toy
#' # Add an unrelated exposure
#' toy2$exposure2 <- round(rnorm(nrow(toy)),2)
#' # Add a grouping
#' grouping <- unique(toy2$id) %>% as_tibble() %>% rename(id=value)
#' grouping <- grouping %>% mutate(group = 1 + (runif(n()) <0.5)*1)
#'
#' # Modify the outcome to add random effects
#' add_to_outcome <- grouping %>% mutate(to_add=ifelse(group==1, -0.5, +0.5)) %>% dplyr::select(id, to_add)
#' toy2 <- toy2 %>% left_join(add_to_outcome) %>% mutate(outcome = outcome + to_add) %>% select(-to_add)
#'
#'
#' output <- fRLM(data=toy2,
#'                id = "id",
#'                time="age",
#'                exposures=c("exposure", "exposure2"),
#'                grouping = grouping,
#'                outcome="outcome",
#'                warmup = 1000, iter = 2000, chains = 2) # this is passed to stan
#'
#'
#' samples <- rstan::extract(output$fit, permuted = TRUE)
#' # Create confidence intervals for delta (effect size of exposures)
#' apply(output$delta, 2, function(x) quantile(x, c(0.025, 0.975)))
#'
#' plot(output)

#' @export

fRLM <- function(data, id, time, exposures, outcome, controls=NULL, grouping=NULL, L=5, alpha_par=1, beta_par=1, ...) {
  # Create the time scaler. This function scales the data down and back to the original scale
  timeScaler <- timeScaler_ff(data %>% pull(!!sym(time)), min=0.06, max=0.93)
  # Standardize time
  data <- data %>% mutate(!!time := timeScaler(data[[time]]))

  # Define the grid
  grid <- seq(0,1, l = 150 )
  grid_original_scale <- timeScaler(grid, original_scale = TRUE)

  # Fit the gaussian processes for each exposure
  condMu <- list()
  for (exposure in exposures) {
    # extract outcome for each exposure
    is_NA <- is.na(data[[exposure]])
    subdata <- data %>% dplyr::filter(!is_NA)
    t_obs <- subdata %>% group_by(!!sym(id)) %>% summarise(t_obs = list(!!sym(time))) %>% pull(t_obs)
    expo  <- subdata %>% group_by(!!sym(id)) %>% summarise(expo = list(!!sym(exposure))) %>% pull(expo)
    gpfitList <- lapply( 1:length(t_obs), function(i) gpFit( t_obs[[i]], expo[[i]] ) )
    condMu[[exposure]] <- t( sapply( gpfitList, predict, tnew = grid ) )
  }

  # Extract quantities
  # ------------------
  # Extract y
  # y_with_id is the outcome variable along with the id. it is used thoughout the code to ensure that each line of matrices / vectors correspond to the same respondent_id across all different stan input
  y_with_id <- data %>% group_by(!!sym(id)) %>% summarise(!!outcome := mean(!!sym(outcome)))
  respondent_id <- y_with_id %>% dplyr::select(!!sym(id)) # tibble with the id: they should consistently appear in that order for all relevant stan inputs
  y <- y_with_id %>% dplyr::pull(!!sym(outcome)) # TODO keep it a tibble?

  # Declare the dimensions
  dim <- list(
    "n" = length(y),
    "L" = L,
    "p" = length(exposures),
    "d" = length(controls) + 1
  )


  # Compile the array of exposures
  basis <- getBasis( L, grid )
  eta_list <- lapply(condMu, function(X) X %*% basis / ncol(X))
  # get them in a p x n x L array
  eta <- abind::abind(eta_list, along=0)

  # Compile the controls
  # We take the mean by id (it should be time-varying anyway)
  C <- data %>% group_by(!!sym(id)) %>% summarise(across(all_of(controls), ~mean(.))) %>% ungroup() %>%
    # Making sure the ordering is consistent across
    right_join(respondent_id) %>%
    # select the controls
    select(all_of(controls)) %>%
    # add intercept
    mutate(intercept=1, .before=1) %>% as.matrix

  # Stan: Bayesian regression
  # -------------------------

  data_stan <- c(dim, list(
    y = y,
    eta = eta,
    C = C,
    alpha_par = alpha_par,
    beta_par = beta_par
  ))


  if(is.null(grouping)) {
    fileName <- "fRLM_additive.stan"
    stanFile <- system.file("stan", fileName, package = "fRLM")
    fit <- rstan::stan( file = stanFile, data = data_stan, ...)
    # output:
    alpha <- rstan::extract(fit, 'alpha')[[1]]
    beta  <- rstan::extract(fit, 'beta')[[1]]
    delta <- rstan::extract(fit, 'delta')[[1]]
    sigma <- rstan::extract(fit, 'sigma')[[1]]

    out <- list(
      fit = fit,
      condMu = condMu,
      delta = delta,
      alpha = alpha,
      beta = beta,
      sigma = sigma,
      grid=grid,
      grid_original_scale=grid_original_scale,
      basis = basis,
      L = L,
      data=data,
      timeScaler=timeScaler,
      dim=dim,
      exposures=exposures)

  } else {
    # Create the grouping matrix and data relevant to groupings
    groups_matrix <- respondent_id %>% left_join(grouping, by = join_by(!!sym(id))) %>% select(-!!sym(id)) %>% as.matrix()

    # add dimensions
    grouping_T <- dim(groups_matrix)[2]
    grouping_S <- length(unique(grouping$group))

    data_stan <- c(data_stan, list(T = grouping_T, S= grouping_S, Sigma=groups_matrix))


    fileName <- "fRLM_additive_random-effects.stan"
    stanFile <- system.file("stan", fileName, package = "fRLM")
    fit <- rstan::stan( file = stanFile, data = data_stan, ...)
    # output:
    alpha <- rstan::extract(fit, 'alpha')[[1]]
    beta  <- rstan::extract(fit, 'beta')[[1]]
    delta <- rstan::extract(fit, 'delta')[[1]]
    sigma <- rstan::extract(fit, 'sigma')[[1]]
    sigma_xi <- rstan::extract(fit, "sigma")[[1]]
    xi    <- rstan::extract(fit, "xi")[[1]]
    rho   <- rstan::extract(fit, "rho")[[1]]
    out <- list(
      fit = fit,
      condMu = condMu,
      delta = delta,
      alpha = alpha,
      beta = beta,
      sigma = sigma,
      grid=grid,
      grid_original_scale=grid_original_scale,
      basis = basis,
      L = L,
      data=data,
      timeScaler=timeScaler,
      dim=dim,
      exposures=exposures,
      xi=xi,
      rho=rho,
      sigma_xi=sigma_xi)

  }

  class(out) <- "additive_fRLM"

  out$w <- predict(out)

  return(out)
}

#' Time Scaler Function Factory
#'
#' Returns a function to scale a time dataset down and back up.
#' The returned function can either standardize the time data to the specified
#' min and max values or scale it back to its original scale.
#'
#' @param time A numeric vector representing the time data.
#' @param min A numeric value representing the minimum value to which the time data should be scaled. Default is 0.
#' @param max A numeric value representing the maximum value to which the time data should be scaled. Default is 1.
#'
#' @return A function that can be used to scale the time data either to the specified range or back to its original scale.
#'
#' @details
#' The function first checks that the input parameters are of the correct type and that min is less than max.
#' It then computes the minimum and maximum values of the input time data and defines a new function, `timeScaler`,
#' that can be used to scale the time data either to the specified range or back to its original scale.
#' The `timeScaler` function takes two parameters: `time`, a numeric vector representing the time data to be scaled,
#' and `original_scale`, a logical value indicating whether the time data should be scaled back to its original scale.
#' If `original_scale` is TRUE, the `timeScaler` function will scale the `time` data back to its original scale.
#' If `original_scale` is FALSE, the `timeScaler` function will standardize the `time` data to the specified `min` and `max` values.
#'
#' @note
#' 1. Error Handling: The function includes error handling to check that the `time`, `min`, and `max` parameters are numeric
#'    and that `min` is less than `max`.
#' 2. Scaling: The `timeScaler` function returned by `timeScaler_ff` scales the `time` data to the range 0-1 before
#'    scaling it either to the specified `min` and `max` values or back to its original scale.
#'
#' @examples
#' # Example usage of the timeScaler_ff function
#' time = 1:10
#' timeScaler <- timeScaler_ff(time, min=0, max = 1)
#' time_scaled <- timeScaler(time)
#' time_unscaled <- timeScaler(time_scaled, original_scale=TRUE)
#' print(paste("Time scaled: ", time_scaled))
#' print(paste("Time unscaled: ", time_unscaled))
#'
#' @export
timeScaler_ff <- function(time, min=0, max=1) {
  # Check input parameters
  if(!is.numeric(time)) stop("time must be a numeric vector")
  if(!is.numeric(min) | !is.numeric(max)) stop("min and max must be numeric")
  if(min >= max) stop("min must be less than max")

  # Find the minimum and maximum values of the input time data
  obs_min = min(time)
  obs_max = max(time)

  # Define the timeScaler function
  timeScaler <- function(time, original_scale=FALSE) {
    if (original_scale) {
      # Scale the time data to the range 0-1
      time <- (time - min ) / (max - min)
      # Scale the time data back to its original scale
      time <- time * (obs_max - obs_min) + obs_min
    }  else {
      # Scale the time data to the range 0-1
      time <- (time - obs_min ) / (obs_max - obs_min)
      # Standardize the time data to the specified min and max values
      time <- time * (max - min) + min
    }
    return(time)
  }

  return(timeScaler)
}

# multiply the output
multiply_slices <- function(A, B) {
  # Check dimensions of A and B
  if (dim(A)[2] != dim(B)[1] || dim(A)[3] != dim(B)[2]) {
    stop("Dimensions of A and B are not compatible for multiplication.")
  }

  # Extract the dimensions for result C
  dim1 <- dim(A)[2]
  dim2 <- dim(A)[1]
  dim3 <- dim(B)[1]

  # Initialize C
  C <- array(0, dim = c(dim1, dim2, dim3))

  # Perform matrix multiplication for each slice
  for (i in 1:dim1) {
    C[i,,] <- A[,i,] %*% t(B)
  }

  return(C)
}





# S3 methods

#' @export
predict.additive_fRLM <- function( object, returnALL = FALSE ){
  out <- list()
  for (i in seq_along(object$exposures)) {
    L <- object$L
    basis <- object$basis
    beta <- object$beta[,i,]
    omega_all <- basis %*% t(beta)
    omega <- rowMeans(omega_all)
    omega_median <- apply(omega_all, 1, median)
    omega_ci <- apply( omega_all, 1, quantile, probs = c(0.025, 0.975) )
    out_i <- list( omega = omega, omega_ci = omega_ci, omega_median = omega_median)
    if(returnALL){
      out_i$omega_all <- omega_all
    }
    out[[object$exposures[i]]] <- out_i
  }
  return( out )
}


#' @export
plot.additive_fRLM <- function(object, ...){
  grid <- object$grid_original_scale
  pred <- predict(object)
  y_max <- max(unlist(lapply(pred, function(x) max(unlist(x)))))
  y_min <- min(unlist(lapply(pred, function(x) min(unlist(x)))))
  if (length(pred) > 1) {
    par(ask = TRUE)
  }
  for (i in 1:length(pred)) {
    plot_name = paste0("Relative importance of life course exposure to ", object$exposures[i], ".")
    plot( pred[[i]]$omega ~ grid, lwd = 2, type = "l", col = "blue", ylim=c(y_min, y_max), main=plot_name, xlab="age", ylab="",... )
    lines( pred[[i]]$omega_ci[1, ] ~ grid, col = "lightblue" )
    lines( pred[[i]]$omega_ci[2, ] ~ grid, col = "lightblue" )
  }
  if (length(pred) > 1) {
    par(ask = FALSE)
  }
}


