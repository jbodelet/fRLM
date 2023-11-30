library(fRLM)
library(dplyr)

# 1) Setup data:
data(toy)
toy2 <- toy
toy2[, 4] <- round( ( toy[, 4] - 10 ) / 30, 2 ) # i. rescale age:
y <- toy2 %>% group_by(id) %>% summarise(outcome = mean(outcome) ) %>% pull(outcome) # extract outcome:
t_obs <- toy2 %>% group_by(id) %>% summarise(age = list(age)) %>%  pull(age) # extract age
exposure <- toy2 %>% group_by(id) %>% summarise(exposure = list(exposure)) %>%  pull(exposure) # extract exposure

# 2) run lpfr
lp <- lpfr( y, t_obs, xobs = exposure, compiled_file = "./inst/stan/lpfr.stan")

# fit <- fRLM_lpfr(data = toy,
#           id="id",  # string, name of the subject identifier
#           exposures = "exposure", # string, the variable name of the exposures (for now, one)
#           outcome="outcome", # string, the variable name of the outcome
#           time = "age" # string, the variable name of the time at which measures were taken
#           )
