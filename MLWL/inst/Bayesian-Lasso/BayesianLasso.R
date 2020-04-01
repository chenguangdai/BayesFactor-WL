##### Normalizing constant estimation for Bayesian Lasso
rm(list = ls())
library(MLWL)
library(geoR)
library(statmod)
library(mvtnorm)

##### set a correct path
setwd(paste(getwd(), '/inst/Bayesian-Lasso', sep = ""))

##### load in simulated data
SNR <- 0.1
lambda <- 20
load(paste(getwd(), "/data/XY/BayesianLasso_SNR_", SNR, ".RData", sep = ""))

##### configuration
source("configuration.R")
num_iterations <- 10^4
num_burnin <- 2000
learning_rate <- 1
flatness_criterion <- 0.2
x0 <- surrogate_mu

##### run the Wang-Landau algorithm
start_time <- Sys.time()
log_marginal_likelihood <- MLWL(num_iterations, num_burnin, target_kernel, surrogate_kernel,
                                target_logdensity, surrogate_logdensity, learning_rate, flatness_criterion, x0)
end_time <- Sys.time()
running_time <- end_time - start_time
print(paste("The log marginal likelihood is", log_marginal_likelihood))
print(paste("The running time is", running_time))








