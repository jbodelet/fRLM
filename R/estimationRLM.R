fitRLM <- function(y, x, ...){
  dat <- list( n = length(y), J = ncol(x), y = y, x = x )
  fileName <- "discreteModel.stan"
  stanFile <- system.file("stan", fileName, package = "fRLM")
  fit <- rstan::stan( file = stanFile, data = dat, ...)
  edit_file(stanFile)
}
