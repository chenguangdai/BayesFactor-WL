setwd('..')

### load in data
source("load_data.R")

### the target distribution
target <- list()
target$logdensity <- function(x) debiasedhmc:::logistic_logtarget_c(x, Y, X, lambda)
target$gradlogdensity <- function(x) debiasedhmc:::logistic_gradlogtarget_c(x, Y, X, lambda)

### the surrogate distribution
surrogate <- list()
load(paste("data/posterior_moments_lambda_", lambda, ".RData", sep = ""))
surrogate_mean <- data_save$mean
surrogate_cov <- data_save$cov
surrogate_precision <- solve(surrogate_cov)
surrogate_precision_chol <- t(chol(surrogate_precision))
surrogate$logdensity <- function(x) return(fast_dmvnorm_chol_inverse(matrix(x, nrow = 1), surrogate_mean, surrogate_precision_chol))
surrogate$rinit <- function() return(fast_rmvnorm(n = 1, mean = surrogate_mean, covariance = surrogate_cov))

### HMC
stepsize <- 0.03
nsteps <- 10
### leapfrog
leapfrog <- function(x, v, gradient){
  v <- v + stepsize * gradient(x) / 2
  for (step in 1:nsteps){
    x <- x + stepsize * v
    if (step != nsteps){
      v <- v + stepsize * gradient(x)
    }
  }
  v <- v + stepsize * gradient(x) / 2
  return(list(x = x, v = v))
}

target$kernel <- function(x){
  x_current <- x
  v_current <- matrix(rnorm(dimension), nrow = 1)
  ### leapfrog steps
  leapfrog_result <- leapfrog(x_current, v_current, target$gradlogdensity)
  v_propose <- - leapfrog_result$v
  x_propose <- leapfrog_result$x
  loglik_current <- target$logdensity(matrix(x_current, nrow = 1))
  loglik_propose <- target$logdensity(matrix(x_propose, nrow = 1))
  ### metropolis acceptance step
  logacceptratio <- loglik_propose - loglik_current + sum(v_current^2)/2 - sum(v_propose^2)/2
  if(log(runif(1)) < logacceptratio){x = x_propose}
  return(x)
}

### directly sampling from the surrogate distribution
surrogate$kernel <- function(x) return(surrogate$rinit())

setwd('./WL')
