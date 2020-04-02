setwd('..')

### load in C functions
sourceCpp("code/mvnorm.cpp", verbose = F)

### load in data and set up parameters
source("code/load_data.R")

### define the target distribution
### prior distribution
prior <- list()
source("code/prior.R")
prior$logdensity <- function(x) return(dmvnorm_cholesky_inverse(x, prior_mean, prior_precision_chol))

### likelihood
source("code/likelihood.R")

### target
target <- list()
target$logdensity <- function(x) return(prior$logdensity(x) + loglikelihood_matrix(x))

### define the surrogate distribution
load(paste("mode/", ngrid, ".RData", sep = ""))
surrogate_mu <- data_save$surrogate_mu
if(ngrid == 10){
  surrogate_sigma <- rep(1.0, dimension)
}else if(ngrid == 20){
  surrogate_sigma <- rep(1.2, dimension)
}else{
  surrogate_sigma <- rep(1.3, dimension)
}
surrogate <- list()
surrogate_precision <- diag(1/surrogate_sigma^2)
surrogate_precision_chol <- t(chol(surrogate_precision))
surrogate$logdensity <- function(x) return(dmvnorm_cholesky_inverse(x, surrogate_mu, surrogate_precision_chol))
surrogate$rinit <- function(n){
  xparticles <- matrix(0, nrow = n, ncol = dimension)
  for(j in 1:dimension){
    xparticles[, j] <- rnorm(n, mean = surrogate_mu[j], sd = surrogate_sigma[1])
  }
  return(xparticles)
} 

setwd('./IS')
