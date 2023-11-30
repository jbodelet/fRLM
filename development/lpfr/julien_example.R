library(fRLM)
library(dplyr)
# parallel computing for the bayesian fit:
# library("rstan") # observe startup messages
# options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = TRUE)

# 1) Setup data:
head(toy) # toy example dataset
toy2 <- toy
toy2[, 4] <- round( ( toy[, 4] - 10 ) / 30, 2 ) # i. rescale age:
y <- toy2 %>% group_by(id) %>% summarise(outcome = mean(outcome) ) %>% pull(outcome) # extract outcome:
t_obs <- toy2 %>% group_by(id) %>% summarise(age = list(age)) %>%  pull(age) # extract age
exposure <- toy2 %>% group_by(id) %>% summarise(exposure = list(exposure)) %>%  pull(exposure) # extract exposure

# 2) Fit gaussian processes:
gpfitList <- lapply( 1:length(y), function(i) gpFit( t_obs[[i]], exposure[[i]] ) )
grid <- seq(0,1, l = 150 )
condMu <- t( sapply( gpfitList, predict, tnew = grid ) )

# 3) Bayesian Estimation:
fitBayes <- funcReg_bayes( y, condMu, warmup = 500, iter = 1000, chains = 2 )
plot(fitBayes)


# 4) splines:
lp <- lpfr( y, t_obs, xobs = exposure )
plot(lp)
