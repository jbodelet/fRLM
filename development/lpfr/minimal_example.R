library(fRLM)
library(dplyr)

data(toy)

# Gaussian outcome
fit <- fRLM_lpfr(data = toy,
          id="id",  # string, name of the subject identifier
          exposures = "exposure", # string, the variable name of the exposures (for now, one)
          outcome="outcome", # string, the variable name of the outcome
          time = "age" # string, the variable name of the time at which measures were taken
          )

plot(fit)
w = predict(fit)

# Binary outcome and bounded exposure (between 0 and 1)
toy <- toy %>% group_by(id) %>% mutate(outcome_binary = rbinom(1,1,1/(1+exp(-outcome)))) %>% ungroup()

fit_binary <- fRLM_lpfr(data = toy,
          id="id",  # string, name of the subject identifier
          exposures = "exposure", # string, the variable name of the exposures (for now, one)
          outcome="outcome_binary", # string, the variable name of the outcome
          family="binary",
          time = "age",
          bounded_exposures = TRUE,
          # warmup=5,
          # iter=10,
          # chains=1
          )

plot(fit_binary)
w_binary = predict(fit_binary)
