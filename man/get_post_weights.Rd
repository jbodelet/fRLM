% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testing.R
\name{get_post_weights}
\alias{get_post_weights}
\title{Get posterior distribution of the weights corresponding to omega}
\usage{
get_post_weights(fitBayes, paritions)
}
\arguments{
\item{fitBayes}{An object created by fReg_bayes.}

\item{paritions}{A matrix with first column being the lower bounds and second column the upper bounds of the periods.
The number of rows provides the number of periods.}

\item{M}{Number of points to approximate the integrals.}
}
\value{
The posterior distributions of the weights.
}
\description{
Compute the Finest Credible Intervals (FCI) for the sensitive models.
}
\examples{
library(latFunReg)
library(dplyr)
library(GPfit)

# 1) Set-up data:
head(toy) # toy example dataset
toy2 <- toy
toy2[, 4] <- round( ( toy[, 4] - 10 ) / 30, 2 ) # i. rescale age:
y <- toy2 \%>\% group_by(id) \%>\% summarise(outcome = mean(outcome) ) \%>\% pull(outcome) # extract outcome:
t_obs <- toy2 \%>\% group_by(id) \%>\% summarise(age = list(age)) \%>\%  pull(age) # extract age
exposure <- toy2 \%>\% group_by(id) \%>\% summarise(exposure = list(exposure)) \%>\%  pull(exposure) # extract exposure

# 2) Fit Gaussian processes
gpfit <- lapply( 1:length(y), function(i) GPfit::GP_fit( t_obs[[i]], exposure[[i]] ) )

# 3) Karhunen loeve expansion:
kl <- KL_expansion( gpfit, t_obs )

# 4) Estimation (em algorithm):
est2 <- fReg_bayes(y, kl, L = 6 )
n_period <- 5
paritions <- cbind( 1:n_period - 1, 1:n_period ) / n_period
colnames( paritions) <- c("l", "u")
post_w <- get_post_weights(est2, paritions )
composite_test( post_w )
get_FCI( post_w )
}
