### marginal likelihood estimation for the Log-Gaussian-Cox process using the WL mixture method
rm(list = ls())
library(MLWL)
library(debiasedhmc)
library(spatstat)

### define the number of grids
ngrid <- 10

##### configuration
source("configuration.R")
num_iterations <- 10^5
num_burnin <- num_iterations/2
learning_rate <- 1
flatness_criterion <- 0.2
x0 <- surrogate_mu

##### run the WL mixture algorithm
start_time <- Sys.time()
log_marginal_likelihood <- MLWL(num_iterations, num_burnin, target_kernel, surrogate_kernel,
                                target_logdensity, surrogate_logdensity, learning_rate, flatness_criterion, x0)
end_time <- Sys.time()
running_time <- end_time - start_time
data_save <- list(lognormconst = log_marginal_likelihood, running_time = running_time)





