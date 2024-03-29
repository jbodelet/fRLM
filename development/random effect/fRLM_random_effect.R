# Example usage of the fRLM function
devtools::load_all()
library(fRLM)
library(tidyverse)
toy <- fRLM::toy
set.seed(1234)

# Add an unrelated exposure
toy$exposure2 <- round(rnorm(nrow(toy)),2)
# Add a grouping
grouping <- unique(toy$id) %>% as_tibble() %>% rename(id=value)
grouping <- grouping %>% mutate(group = 1 + (runif(n()) <0.5)*1)

# Modify the outcome to add random effects
add_to_outcome <- grouping %>% mutate(to_add=ifelse(group==1, -0.5, +0.5)) %>% dplyr::select(id, to_add)
toy <- toy %>% left_join(add_to_outcome) %>% mutate(outcome = outcome + to_add) %>% select(-to_add)


toy <- toy %>% rename(y=outcome)

toy$exposure[c(10, 25, 38, 52)] <- NA
toy$exposure2[c(18, 33, 38, 52)] <- NA

output <- fRLM(data=toy,
               id = "id",
               time="age",
               exposures=c("exposure", "exposure2"),
               grouping = grouping,
               outcome="y",
               warmup = 1000, iter = 2000, chains = 2) # this is passed to stan


samples <- rstan::extract(output$fit, permuted = TRUE)
# Create confidence intervals for delta (effect size of exposures)
apply(output$delta, 2, function(x) quantile(x, c(0.025, 0.975)))

plot(output)
