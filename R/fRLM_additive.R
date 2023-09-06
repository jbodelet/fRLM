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
#' # Example usage of the fRLM function
#' library(fRLM)
#' data(toy)
#' output <- fRLM(data=toy,
#'   id = "id",
#'   time="age",
#'   exposures="exposure",
#'   outcome="outcome",
#'   warmup = 500, iter = 1000, chains = 2) # this is passed to stan
#'
#' @export

fRLM <- function(data, id, time, exposures, outcome, controls=NULL, L=5, alpha_par=1, beta_par=1, ...) {
  # Create the time scaler. This function scales the data down and back to the original scale
  timeScaler <- timeScaler_ff(data %>% pull(!!sym(time)), min=0.05, max=0.95)
  # Standardize time
  data <- data %>% mutate(!!time := timeScaler(data[[time]]))

  # Define the grid
  grid <- seq(0,1, l = 150 )

  # Fit the gaussian processes for each exposure
  condMu <- list()
  for (exposure in exposures) {
    # extract outcome for each exposure
    subdata <- data %>% filter(!is.na({{exposure}}))
    t_obs <- subdata %>% group_by(!!sym(id)) %>% summarise(t_obs = list(!!sym(time))) %>% pull(t_obs)
    expo  <- subdata %>% group_by(!!sym(id)) %>% summarise(expo = list(!!sym(exposure))) %>% pull(expo)
    gpfitList <- lapply( 1:length(t_obs), function(i) gpFit( t_obs[[i]], expo[[i]] ) )
    condMu[[exposure]] <- t( sapply( gpfitList, predict, tnew = grid ) )
  }

  # Extract quantities
  # ------------------
  # Extract y
  y <- data %>% group_by(!!sym(id)) %>% summarise(outcome = mean(!!sym(outcome))) %>% pull(outcome)

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
  C <- data %>% group_by(!!sym(id)) %>% summarise(across(all_of(controls), ~mean(.))) %>% ungroup() %>%
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
  fileName <- "fRLM_additive.stan"
  stanFile <- system.file("stan", fileName, package = "fRLM")
  fit <- rstan::stan( file = stanFile, data = data_stan, ...)
  # output:
  alpha <- rstan::extract(fit, 'alpha')[[1]]
  beta  <- rstan::extract(fit, 'beta')[[1]]
  delta <- rstan::extract(fit, 'delta')[[1]]
  sigma <- rstan::extract(fit, 'sigma')[[1]]
  out <- list( fit = fit, condMu = condMu, delta = delta, alpha = alpha, beta = beta, sigma = sigma, basis = basis, L = L, data=data, timeScaler=timeScaler)
  class(out) <- "funcRegBayes"

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


