### WL mixture method for Bayesian Lasso
rm(list = ls())
library(MLWL)
library(geoR)
library(statmod)
library(mvtnorm)

### SNR and lambda
SNR <- 3
lambda <- 5

### configuration
source("configuration.R")
num_iterations <- 10^4
num_burnin <- 2000
learning_rate <- 1
flatness_criterion <- 0.2
x0 <- surrogate_mu

### the WL mixture method
start_time <- Sys.time()
result <- MLWL(num_iterations, num_burnin, target$kernel, surrogate$kernel, target$logdensity, 
               surrogate$logdensity, learning_rate, flatness_criterion, x0, suppress = F)
end_time <- Sys.time()
lognormconst <- result$log_marginal_likelihood
running_time <- end_time - start_time
