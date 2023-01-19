# Packages:
# devtools::install_github("Jaeoc/lmhyp")
# install.packages("combinat")
# install.packages("ggformula")
# install.packages("matrixStats") # for mcmc fco



#' Bayesian test for the three main hypotheses,
#' 
#' Compute the posterior probabilities of the main different models: Accumulation, Sensitivive, Critical.
#' It computes the posterior distribution of phi (the range of the weights ), and compare it to a threshold.
#' 
#' @param post_w Posterior samples of the weights.
#' @param l  lower treshold for phi  (under which we accept the accumulation model).
#' @param u  upper treshold for phi  (above which we accept the critical model).
#' @return   The posterior probabilities of each three hypotheses.
#' @export
#' @examples
#' library(latFunReg)
#' library(dplyr)
#' library(GPfit)
#' 
#' # 1) Set-up data:
#' head(toy) # toy example dataset
#' toy2 <- toy
#' toy2[, 4] <- round( ( toy[, 4] - 10 ) / 30, 2 ) # i. rescale age:
#' y <- toy2 %>% group_by(id) %>% summarise(outcome = mean(outcome) ) %>% pull(outcome) # extract outcome:
#' t_obs <- toy2 %>% group_by(id) %>% summarise(age = list(age)) %>%  pull(age) # extract age
#' exposure <- toy2 %>% group_by(id) %>% summarise(exposure = list(exposure)) %>%  pull(exposure) # extract exposure
#' 
#' # 2) Fit Gaussian processes
#' gpfit <- lapply( 1:length(y), function(i) GPfit::GP_fit( t_obs[[i]], exposure[[i]] ) )
#' 
#' # 3) Karhunen loeve expansion:
#' kl <- KL_expansion( gpfit, t_obs )
#' 
#' # 4) Estimation (em algorithm):
#' est2 <- fReg_bayes(y, kl, L = 6 )
#' #==============
#' # Testing:
#' #==============
#' library(lmhyp)
#' library(stringi)
#' library(stringr)
#' library(purrr)
#' library(tidyr)
#' library(tibble)
#' 
#' # composite test: Accumulation, sensitive, critical
#' n_period <- 5
#' paritions <- cbind( 1:n_period - 1, 1:n_period ) / n_period
#' colnames( paritions) <- c("l", "u")
#' post_w <- get_post_weights(est2, paritions )
#' composite_test( post_w )
#' 
#' # Finest credible intervals:
#' get_FCI( post_w )
composite_test <- function( post_w, l = 0.15, u = 0.85 ){
  n_samp <- nrow(post_w)
  post_w <- as_tibble(post_w)
  aug  <- bind_cols( post_w) %>% 
    mutate(post_range = apply(post_w, 1, max) - apply(post_w, 1, min))  
  aug <- mutate(aug, class = cut(post_range, c(0, l, u, 1), labels = c("a", "s", "c")))
  # "MODEL PROBABILITY": CREDIBLE IF EXCEEDS SAY 0.9
  return( summary(aug$class)/nrow(aug) )
}


#' Finest Credible Intervals (FCI)
#' 
#' Compute the Finest Credible Intervals (FCI) for the sensitive models.
#' 
#' @param post_w Posterior samples of the weights.
#' @return   The rankings and their different probabilities.
#' @export
#' @examples
#' library(latFunReg)
#' library(dplyr)
#' library(GPfit)
#' 
#' # 1) Set-up data:
#' head(toy) # toy example dataset
#' toy2 <- toy
#' toy2[, 4] <- round( ( toy[, 4] - 10 ) / 30, 2 ) # i. rescale age:
#' y <- toy2 %>% group_by(id) %>% summarise(outcome = mean(outcome) ) %>% pull(outcome) # extract outcome:
#' t_obs <- toy2 %>% group_by(id) %>% summarise(age = list(age)) %>%  pull(age) # extract age
#' exposure <- toy2 %>% group_by(id) %>% summarise(exposure = list(exposure)) %>%  pull(exposure) # extract exposure
#' 
#' # 2) Fit Gaussian processes
#' gpfit <- lapply( 1:length(y), function(i) GPfit::GP_fit( t_obs[[i]], exposure[[i]] ) )
#' 
#' # 3) Karhunen loeve expansion:
#' kl <- KL_expansion( gpfit, t_obs )
#' 
#' # 4) Estimation (em algorithm):
#' est2 <- fReg_bayes(y, kl, L = 6 )
#' #==============
#' # Testing:
#' #==============
#' library(lmhyp)
#' library(stringi)
#' library(stringr)
#' library(purrr)
#' library(tidyr)
#' library(tibble)
#' 
#' # composite test: Accumulation, sensitive, critical
#' n_period <- 5
#' paritions <- cbind( 1:n_period - 1, 1:n_period ) / n_period
#' colnames( paritions) <- c("l", "u")
#' post_w <- get_post_weights(est2, paritions )
#' composite_test( post_w )
#' 
#' # Finest credible intervals:
#' get_FCI( post_w )
get_FCI <- function( post_w ){
  n_samp <- nrow(post_w)
  post_w <- as_tibble(post_w)
  # SENSITIVE MODELS:
  ordering <- as.matrix( apply( post_w, 1, function(u) paste0( order(u), collapse = "<" ) ) )
  colnames(ordering) <- "H"
  post <- as_tibble( ordering ) %>% group_by(H) %>% 
    count %>% mutate(Hp = n/n_samp, M = "M0") %>% 
    dplyr::select(-n, Hp, H, M) %>% ungroup %>% list
  # ALGORITHM 1
  max_ranking <- get_fco_det(post) # initialize at max probability full ranking
  if(nrow(max_ranking) != 1){ #  why do I have several max rankings sometimes?
    warning("warning: several max_rankings! Pick the first one.")
    max_ranking <- max_ranking[1, ]    
  }
  # print( paste( "maxRand is", nrow(max_ranking) ) )  # should I print this message?
  fco_recurse <- find_local_fco(partial_rank = max_ranking$H, dat = post[[1]] )
  CDF  <- tidy_output( max_ranking, fco_recurse )
  # ALGORITHM 2 (more complex: to be depricated?)
  recursion <- purrr::map(list(post_w), get_discrete_posterior, type = "greedy") 
  CDF <- recursion %>% 
    map(pluck, "out_wux") %>% # over all params
    map(filter, what == "max" )
  return( CDF )
}


#' Get posterior distribution of the weights corresponding to omega
#' 
#' Compute the Finest Credible Intervals (FCI) for the sensitive models.
#' 
#' @param fitBayes An object created by fReg_bayes.
#' @param paritions A matrix with first column being the lower bounds and second column the upper bounds of the periods.
#'              The number of rows provides the number of periods.
#' @param M Number of points to approximate the integrals.            
#' @return   The posterior distributions of the weights.
#' @export
#' @examples
#' library(latFunReg)
#' library(dplyr)
#' library(GPfit)
#' 
#' # 1) Set-up data:
#' head(toy) # toy example dataset
#' toy2 <- toy
#' toy2[, 4] <- round( ( toy[, 4] - 10 ) / 30, 2 ) # i. rescale age:
#' y <- toy2 %>% group_by(id) %>% summarise(outcome = mean(outcome) ) %>% pull(outcome) # extract outcome:
#' t_obs <- toy2 %>% group_by(id) %>% summarise(age = list(age)) %>%  pull(age) # extract age
#' exposure <- toy2 %>% group_by(id) %>% summarise(exposure = list(exposure)) %>%  pull(exposure) # extract exposure
#' 
#' # 2) Fit Gaussian processes
#' gpfit <- lapply( 1:length(y), function(i) GPfit::GP_fit( t_obs[[i]], exposure[[i]] ) )
#' 
#' # 3) Karhunen loeve expansion:
#' kl <- KL_expansion( gpfit, t_obs )
#' 
#' # 4) Estimation (em algorithm):
#' est2 <- fReg_bayes(y, kl, L = 6 )
#' n_period <- 5
#' paritions <- cbind( 1:n_period - 1, 1:n_period ) / n_period
#' colnames( paritions) <- c("l", "u")
#' post_w <- get_post_weights(est2, paritions )
#' composite_test( post_w )
#' get_FCI( post_w )
get_post_weights <- function( fitBayes, paritions ){
  omega <- predict(fitBayes, returnALL = TRUE )$omega_all
  M <- nrow( omega )
  Mgrid <- seq(0,1, l = M )
  post_omega <- apply( omega, 2, function(u) approxfun(Mgrid, y = u ) )
  post_w <- t( sapply( post_omega, function(om) apply( paritions, 1, function(x) integrate( om, lower = x[1], upper = x[2] )$value ) ) )
  colnames( post_w ) <- paste0("w", 1:ncol(post_w) )
  post_w <- round( post_w, 4)
}
