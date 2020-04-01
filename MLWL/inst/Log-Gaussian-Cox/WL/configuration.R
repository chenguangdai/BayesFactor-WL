setwd('..')
### load in data and set up parameters
source("code/load_Data.R")

### load in C functions.
sourceCpp("code/mvnorm.cpp", verbose = F)

### define the target distribution
### prior distribution
source("code/prior.R")

### likelihood
source("code/likelihood.R")

### target distribution
target <- list()
target$logdensity <- function(x) return(prior$logdensity(x) + loglikelihood(x))
target$gradlogdensity <- function(x) return(prior$gradlogdensity(x) + gradloglikelihood(x))

### an HMC kernel invariant to the target distribution
source("code/HMC.R")
target$kernel <- construct_kernel(target)

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
surrogate$logdensity <- function(x){
  return(-dimension/2*log(2*pi) - sum(log(surrogate_sigma)) - 0.5 * sum((x - surrogate_mu)^2/surrogate_sigma^2))
} 

### direct sampling from the surrogate distribution
surrogate$kernel <- function(x){
  return(rnorm(dimension, mean = surrogate_mu, sd = surrogate_sigma))
}

setwd("./WL")


