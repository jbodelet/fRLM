
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fRLM

<!-- badges: start -->
<!-- badges: end -->

The functional Relevant Life course Model (fRLM) is a bayesian approach
to test life course hypotheses over the continuous time.

## Installation

You can install the development version of fRLM from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jbodelet/fRLM")
```

## Example

``` r
library(fRLM)
library(dplyr)
# parallel computing for the bayesian fit:
# library("rstan") # observe startup messages
# options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = TRUE)
```

### 1) Data setup:

``` r
head(toy) # toy example dataset
toy2 <- toy
toy2[, 4] <- round( ( toy[, 4] - 10 ) / 30, 2 ) # i. rescale age:
y <- toy2 %>% group_by(id) %>% summarise(outcome = mean(outcome) ) %>% pull(outcome) # extract outcome:
t_obs <- toy2 %>% group_by(id) %>% summarise(age = list(age)) %>%  pull(age) # extract age
exposure <- toy2 %>% group_by(id) %>% summarise(exposure = list(exposure)) %>%  pull(exposure) # extract exposure
```

### 2) Fit gaussian processes:

``` r
gpfitList <- lapply( 1:length(y), function(i) gpFit( t_obs[[i]], exposure[[i]] ) )
plot(gpfitList[[1]])
```

<img src="man/figures/README-Fit gaussian processes-1.png" width="100%" />

``` r
grid <- seq(0,1, l = 150 )
pred <- t( sapply( gpfitList, predict, tnew = grid ) )
```

### 3) Frequentist Estimation:

``` r
fit <- funcReg( y, pred, grid = grid, L = 8 )
bootfit <- bootFunReg( y, pred, grid = grid, L = 8 )
plot(bootfit, ylim = c(0, 3))
```

<img src="man/figures/README-Frequentist Estimation-1.png" width="100%" />

### 5) Bayesian Estimation:

``` r
fitBayes <- funcReg_bayes( y, pred, warmup = 250, iter = 500, chains = 2 )
plot(fitBayes)
```

<img src="man/figures/README-Bayesian Estimation-1.png" width="100%" />

### Testing:

``` r
is_lmhyp_available <- require("lmhyp")
if(!is_lmhyp_available){
  devtools::install_github("Jaeoc/lmhyp")  
}
library(lmhyp)
library(stringi)
library(stringr)
library(purrr)
library(tidyr)
library(tibble)

# composite test: Accumulation, sensitive, critical
n_period <- 3
paritions <- cbind( 1:n_period - 1, 1:n_period ) / n_period
colnames( paritions) <- c("l", "u")
post_w <- get_post_weights( fitBayes, paritions )
composite_test( post_w )

# Finest credible intervals:
get_FCI( post_w )
```
