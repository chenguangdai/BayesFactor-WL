### Importance sampling for the Log-Gaussian Cox process
rm(list = ls())
library(MASS)
library(debiasedhmc)
library(spatstat)

### define the number of grids
ngrid <- 10

### configuration
source("configuration.R")

### importance sampling
start_time <- Sys.time()
num_samples <- 10^6
xparticles <- surrogate$rinit(num_samples)
logw <- rep(0, num_samples)
logw <- target$logdensity(xparticles) - surrogate$logdensity(xparticles)
maxlw <- max(logw)
lognormconst <- log(mean(exp(logw - maxlw))) + maxlw
end_time <- Sys.time()
running_time <- end_time - start_time


