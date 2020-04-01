##### load in C functions.
sourceCpp("mvnorm.cpp", verbose = F)

##### define the target distribution
##### prior distribution
prior_mean <- rep(parameter_mu, dimension)
prior_cov <- matrix(nrow = dimension, ncol = dimension)
for(m in 1:dimension){
  for(n in 1:dimension){
    index_m <- c(floor((m - 1)/ngrid) + 1, ((m - 1)%%ngrid) + 1)
    index_n <- c(floor((n - 1)/ngrid) + 1, ((n - 1)%%ngrid) + 1)
    prior_cov[m, n] <- parameter_sigmasq * exp(-sqrt(sum((index_m - index_n)^2))/(ngrid * parameter_beta))
  }
}
prior_precision <- solve(prior_cov)
prior_precision_chol <- t(chol(prior_precision))
priorlogdensity <- function(x) return(dmvnorm_cholesky_inverse(matrix(x, nrow = 1), prior_mean, prior_precision_chol))
gradpriorlogdensity <- function(x) return(eigenMapMatMult(matrix(parameter_mu - x, nrow = 1),  prior_precision))
##### likelihood
loglikelihood <- function(x) return(coxprocess_loglikelihood(matrix(x, nrow = 1), data_counts, parameter_area))
gradloglikelihood <- function(x) return(matrix(data_counts, nrow = 1) - parameter_area * exp(matrix(x, nrow = 1)))
##### target distribution
target_logdensity <- function(x) return(priorlogdensity(x) + loglikelihood(x))
target_gradlogdensity <- function(x) return(gradpriorlogdensity(x) + gradloglikelihood(x))

##### An HMC kernel invariant to the target distribution
stepsize <- 0.25
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
### HMC
target_kernel <- function(x){
  x_current <- x
  v_current <- matrix(rnorm(ngrid^2), nrow = 1)
  ### leapfrog steps
  leapfrog_result <- leapfrog(x_current, v_current, target_gradlogdensity)
  v_propose <- - leapfrog_result$v
  x_propose <- leapfrog_result$x
  loglik_current <- target_logdensity(x_current)
  loglik_propose <- target_logdensity(x_propose)
  ### metropolis acceptance step
  logacceptratio <- loglik_propose - loglik_current + sum(v_current^2)/2 - sum(v_propose^2)/2
  if(log(runif(1)) < logacceptratio){x <- x_propose}
  return(x)
}

### define the surrogate distribution
load(paste("~/mode/mode", ngrid, ".RData", sep = ""))
surrogate_mu <- data_save$surrogate_mu
if(ngrid == 10){
  surrogate_sigma <- rep(1.0, dimension)
}else if(ngrid == 20){
  surrogate_sigma <- rep(1.2, dimension)
}else{
  surrogate_sigma <- rep(1.3, dimension)
}
surrogate_logdensity <- function(x){
  return(-dimension/2*log(2*pi) - sum(log(surrogate_sigma)) - 0.5 * sum((x - surrogate_mu)^2/surrogate_sigma^2))
} 

##### direct sampling from the surrogate distribution
surrogate_kernel <- function(x){
  return(rnorm(dimension, mean = surrogate_mu, sd = surrogate_sigma))
}




