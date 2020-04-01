### The stepping stone method for the Log-Gaussian Cox process
rm(list = ls())
library(debiasedhmc)
library(spatstat)

### number of grids
ngrid <- 10

### configuration
source("configuration.R")
### load the mode
load(paste("/Users/chenguang/Desktop/working papers/MTM/revision/code/log-cox-process/Stepping Stone/mode/mode", ngrid, ".RData", sep = ""))
x0 <- data_save$surrogate_mu
### load the temeprature sequence
load(paste("/Users/chenguang/Desktop/working papers/MTM/revision/code/log-cox-process/Stepping Stone/sequence/", ngrid, ".RData", sep = ""))
temperature <- data_save$sequence
nparticles <- 1.5*10^3

### Stepping stone
start_time <- Sys.time()
lognormconst_ratio <- rep(NA, length(temperature) - 1)
for(iter in 1:(length(temperature) - 1)){
  ### temperature
  gamma_target <- temperature[iter + 1]
  gamma_proposal <- temperature[iter]
  
  ### generate samples
  if(iter == 1){
    xparticles <- prior$rinit(n = nparticles)
  }else{
    ### define the proposal used in importance sampling
    proposal <- list()
    proposal$logdensity <- function(x) return(prior$logdensity(x) + gamma_proposal*loglikelihood(x))
    proposal$gradlogdensity <- function(x) return(prior$gradlogdensity(x) + gamma_proposal*gradloglikelihood(x))
    
    ### define the proposal kernel
    proposal_kernel <- construct_kernel(proposal)
    
    ### generate samples using MCMC
    xparticles <- matrix(0, nrow = nparticles, ncol = dimension)
    xparticles[1, ] <- x0
    for(i in 2:nparticles){xparticles[i, ] <- proposal_kernel(matrix(xparticles[i - 1, ], nrow = 1))}
  }
  
  ### calculate the log normalizing constant ratio
  logw <- (gamma_target - gamma_proposal) * loglikelihood(xparticles)
  maxlw <- max(logw)
  lognormconst_ratio[iter] <- log(mean(exp(logw - maxlw))) + maxlw
  
  print(iter/length(temperature))
}

end_time <- Sys.time()


lognormconst <- sum(lognormconst_ratio)
running_time <- end_time - start_time


