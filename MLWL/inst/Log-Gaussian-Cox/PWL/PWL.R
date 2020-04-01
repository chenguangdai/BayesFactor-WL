### Parallel Wang-Landau algorithm for the Log-Gaussian Cox process
rm(list = ls())
library(MLWL)
library(debiasedhmc)
library(spatstat)

### number of grids
ngrid <- 10

### configuration
source("configuration.R")

### algorithmic settings of PWL
num_iterations <- 1.5*10^3
num_burnin <- num_iterations/2
learning_rate <- 1
flatness_criterion <- 0.2

### PWL
start_time <- Sys.time()
lognormconst_ratio <- rep(NA, length(temperature) - 1)
for(iter in 1:(length(temperature) - 1)){
  ### temperature
  gamma_target <- temperature[iter + 1]
  gamma_surrogate <- temperature[iter]
  
  ### define the target
  target <- list()
  target$logdensity <- function(x) return(prior$logdensity(x) + gamma_target*loglikelihood(x))
  target$gradlogdensity <- function(x) return(prior$gradlogdensity(x) + gamma_target*gradloglikelihood(x))
  
  ### define the surrogate
  surrogate <- list()
  surrogate$logdensity <- function(x) return(prior$logdensity(x) + gamma_surrogate*loglikelihood(x))
  surrogate$gradlogdensity <- function(x) return(prior$gradlogdensity(x) + gamma_surrogate*gradloglikelihood(x))
  
  ### define the target kernel
  target_kernel <- construct_kernel(target)
  
  ### define the surrogate kernel
  if(iter == 0){
    surrogate_kernel <- function(x){return(prior$rinit(1))}
  }else{
    surrogate_kernel <- construct_kernel(surrogate)
  }
  
  ### run the WL mixture method
  result <- MLWL(num_iterations, num_burnin, target_kernel, surrogate_kernel,
                                    target$logdensity, surrogate$logdensity, learning_rate, flatness_criterion, x0)
  lognormconst_ratio[iter] <- result$log_marginal_likelihood
}
end_time <- Sys.time()

lognormconst <- sum(lognormconst_ratio)
running_time <- end_time - start_time



