### importance sampling for Bayesian Lasso
rm(list = ls())
library(debiasedhmc)

### SNR and lambda
SNR <- 3
lambda <- 5

### configuration
source("configuration.R")

### importance sampling
start_time <- Sys.time()
num_samples <- 5*10^5
xparticles <- proposal$rinit(num_samples)
logw <- target$logdensity(xparticles) - proposal$logdensity(xparticles)
maxlw <- max(logw)
lognormconst <- log(mean(exp(logw - maxlw))) + maxlw
end_time <- Sys.time()
running_time <- end_time - start_time
