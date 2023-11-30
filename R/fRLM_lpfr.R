#' fRLM Function
#'
#' Wrapper function that estimates a Functional Regression Linear Model (fRLM) with a given set of predictors.
#'
#' @param data A tidy dataframe containing the data.
#' @param id A string representing the identifier variable in the dataframe.
#' @param time A string representing the time variable in the dataframe, at what time the outcome or exposures were measured.
#' @param exposures A character vector of strings representing the exposure variables in the dataframe.
#' @param outcome A string representing the outcome variable in the dataframe.
#' @param family A string representing the family. Either "gaussian" or "binary".
#' @param controls A character vector of strings representing the control variables in the dataframe. Default is NULL.
#' @param L An integer representing the number of basis functions. Default is 4.
#' @param K An integer representing the number of knots for each basis function. Default is K=5.
#' @param ... Additional arguments passed to the `rstan::stan` function.
#'
#' @return A list containing the fitted model and other relevant information.
#'
#' @details
#' The `outcome`and `controls` must be the same (repeated) for all time stamps for each unit. See the examples below.
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
#' output <- fRLM_lpfr(data=toy2,
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

fRLM_lpfr <- function(data, id, time, exposures, outcome, family="gaussian", controls=NULL, L=4, K=5, grid= seq(0,1,l=150), ...) {
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

  # Extract the covariates
  C <- data %>% group_by(!!sym(id)) %>% summarise(across(all_of(controls), ~mean(.))) %>% ungroup() %>%
    # select the controls
    select(all_of(controls)) %>%
    # mutate(intercept=1, .before=1) %>% # add intercept (don't! already is in Julien's code!)
    as.matrix

  # Feed quantities to Julien's function
  # -------------------------

  if (family == "gaussian") {
    file_name <- "lpfr.stan"
  } else if (family == "binary") {
    file_name <- "lpfr_binary.stan"
  }

  stan_file <- system.file("stan", file_name, package = "fRLM")

  fit_list <- lpfr(y = y, tobs=tobs, xobs=xobs, covariates=C, compiled_file =stan_file, L=L, K=K, grid=grid, ...)

  out <- c(fit_list,
    list(
      L=L,
      K=K,
      grid=grid,
      grid_original_scale = timeScaler(grid, original_scale=TRUE),
      timeScaler=timeScaler
    )
  )
  return(out)
}
