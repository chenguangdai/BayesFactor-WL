### MTM-RJMCMC on pollution data
rm(list = ls())
library(mvnfast)
library(LaplacesDemon)
library(pscl)

### the standard deviation of the auxiliary distribution
sigma_auxiliary <- 1

### algorithmic setting of MTM
### number of tries
m <- 5
### standard deviation of jumping distance 
sigma_r <- 1 

### configuration
source("configuration.R")

### MCMC
num_iter <- 5*10^4
beta <- matrix(NA, nrow = num_iter, ncol = p)
model_dimension <- rep(0, p)
margial_prob <- rep(0, p)

### random initialization
initial_variable_index <- sample(x = 1:p, size = p, replace = F)
beta[1, initial_variable_index] <- lm(Y ~ as.matrix(X[, initial_variable_index], nrow = n) - 1)$coefficients
sigma[1] <- update_sigma(beta[1, ])

start_time <- Sys.time()
for(iter in 2:num_iter){
  ### condition on sigma, update beta
  beta[iter, ] <- MTM_RJMCMC(beta[iter - 1, ], sigma[iter - 1])
  
  ### condition on beta, update sigma
  sigma[iter] <- update_sigma(beta[iter, ])
  
  ### burn in the first 10%
  if(iter > 0.1*num_iter){
    selected_variable_index <- which(!is.na(beta[iter, ]))
    margial_prob[selected_variable_index] <- margial_prob[selected_variable_index] + 1
  }
}
end_time <- Sys.time()
margial_prob <- margial_prob/sum(margial_prob)
running_time <- as.numeric(end_time - start_time)