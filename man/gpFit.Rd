% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gaussianProcesses.R
\name{gpFit}
\alias{gpFit}
\title{Gaussian process regression}
\usage{
gpFit(
  tobs,
  Xobs,
  lambda = 0.01,
  optimizationRange = c(-2, 4),
  power = 1.95,
  nugget = 0
)
}
\arguments{
\item{tobs}{Observed time points.}

\item{Xobs}{Observed points of the gaussian process.}

\item{lambda}{penalty.}

\item{optimizationRange}{A vector containing the end-points of the interval to be searched.}

\item{power}{Power parameter for the exponential kernel.}

\item{nugget}{Nugget parameter for stable inversion of the correlation matrix.}
}
\value{
An list containing estimates
}
\description{
Fit a Gaussian process regression with penalyzed likelihood.
}
\examples{
  gp <- generateGaussianProcess()
  tobs <- sort(runif(5))
  Xobs <- gp(tobs)
  gpHat <- gpFit(tobs, Xobs)
  plot(gpHat)
}
