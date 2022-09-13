
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fRLM

<!-- badges: start -->
<!-- badges: end -->

The goal of fRLM is to …

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
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
# parallel computing for the bayesian fit:
# library("rstan") # observe startup messages
# options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = TRUE)
```

### 1) Data setup:

``` r
head(toy) # toy example dataset
#>   id outcome exposure  age
#> 1  1    2.43     0.28 14.5
#> 2  1    2.43     0.37 17.4
#> 3  1    2.43     0.46 19.9
#> 4  1    2.43     0.62 23.5
#> 5  1    2.43     1.03 33.3
#> 6  2    0.13     0.18 11.9
toy2 <- toy
toy2[, 4] <- round( ( toy[, 4] - 10 ) / 30, 2 ) # i. rescale age:
y <- toy2 %>% group_by(id) %>% summarise(outcome = mean(outcome) ) %>% pull(outcome) # extract outcome:
t_obs <- toy2 %>% group_by(id) %>% summarise(age = list(age)) %>%  pull(age) # extract age
exposure <- toy2 %>% group_by(id) %>% summarise(exposure = list(exposure)) %>%  pull(exposure) # extract exposure
```

\###2) Fit gaussian processes:

``` r
gpfitList <- lapply( 1:length(y), function(i) gpFit( t_obs[[i]], exposure[[i]] ) )
plot(gpfitList[[1]])
```

<img src="man/figures/README-Fit gaussian processes-1.png" width="100%" />

``` r
grid <- seq(0,1, l = 150 )
pred <- t( sapply( gpfitList, predict, tnew = grid ) )
```

\###3) Frequentist Estimation:

``` r
fit <- funcReg( y, pred, grid = grid, L = 8 )
bootfit <- bootFunReg( y, pred, grid = grid, L = 8 )
plot(bootfit, ylim = c(0, 3))
```

<img src="man/figures/README-Frequentist Estimation-1.png" width="100%" />

\###5) Bayesian Estimation:

``` r
fitBayes <- funcReg_bayes( y, pred, warmup = 250, iter = 500, chains = 2 )
#> 
#> SAMPLING FOR MODEL 'plugin' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.001 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 10 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 500 [  0%]  (Warmup)
#> Chain 1: Iteration:  50 / 500 [ 10%]  (Warmup)
#> Chain 1: Iteration: 100 / 500 [ 20%]  (Warmup)
#> Chain 1: Iteration: 150 / 500 [ 30%]  (Warmup)
#> Chain 1: Iteration: 200 / 500 [ 40%]  (Warmup)
#> Chain 1: Iteration: 250 / 500 [ 50%]  (Warmup)
#> Chain 1: Iteration: 251 / 500 [ 50%]  (Sampling)
#> Chain 1: Iteration: 300 / 500 [ 60%]  (Sampling)
#> Chain 1: Iteration: 350 / 500 [ 70%]  (Sampling)
#> Chain 1: Iteration: 400 / 500 [ 80%]  (Sampling)
#> Chain 1: Iteration: 450 / 500 [ 90%]  (Sampling)
#> Chain 1: Iteration: 500 / 500 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 1.587 seconds (Warm-up)
#> Chain 1:                0.59 seconds (Sampling)
#> Chain 1:                2.177 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'plugin' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:   1 / 500 [  0%]  (Warmup)
#> Chain 2: Iteration:  50 / 500 [ 10%]  (Warmup)
#> Chain 2: Iteration: 100 / 500 [ 20%]  (Warmup)
#> Chain 2: Iteration: 150 / 500 [ 30%]  (Warmup)
#> Chain 2: Iteration: 200 / 500 [ 40%]  (Warmup)
#> Chain 2: Iteration: 250 / 500 [ 50%]  (Warmup)
#> Chain 2: Iteration: 251 / 500 [ 50%]  (Sampling)
#> Chain 2: Iteration: 300 / 500 [ 60%]  (Sampling)
#> Chain 2: Iteration: 350 / 500 [ 70%]  (Sampling)
#> Chain 2: Iteration: 400 / 500 [ 80%]  (Sampling)
#> Chain 2: Iteration: 450 / 500 [ 90%]  (Sampling)
#> Chain 2: Iteration: 500 / 500 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 1.303 seconds (Warm-up)
#> Chain 2:                0.427 seconds (Sampling)
#> Chain 2:                1.73 seconds (Total)
#> Chain 2:
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
plot(fitBayes)
```

<img src="man/figures/README-Bayesian Estimation-1.png" width="100%" />

\###Testing:

``` r
is_lmhyp_available <- require("lmhyp")
#> Loading required package: lmhyp
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
#> a s c 
#> 0 1 0
# Finest credible intervals:
get_FCI( post_w )
#> [[1]]
#> # A tibble: 2 x 3
#>   what  hyp   value
#>   <chr> <chr> <dbl>
#> 1 max   3|2|1 0.986
#> 2 max   3|2,1 1
```
