### the WL mixture method for Logistic regression
rm(list = ls())
library(MLWL)
library(debiasedhmc)

### lambda
lambda <- 0.01

### configuration
source("configuration.R")

num_iterations <- 5000
num_burnin <- num_iterations/2
learning_rate <- 1
flatness_criterion <- 0.2
x0 <- surrogate_mean

### the WL mixture method
start_time <- Sys.time()
result <- MLWL(num_iterations, num_burnin, target$kernel, surrogate$kernel, target$logdensity, 
               surrogate$logdensity, learning_rate, flatness_criterion, x0, suppress = F)
lognormconst <- result$log_marginal_likelihood
end_time <- Sys.time()
running_time <- end_time - start_time

