### MTM-RJMCMC on pollution data
rm(list = ls())
library(mvnfast)
library(LaplacesDemon)
library(pscl)

### sparsity
g <- exp(10)

### the standard deviation of the proposal distribution for beta
### also the step size in MH for beta
sigma_proposal <- 0.5

### configuration
source("configuration.R")

### MCMC
num_iter <- 1.5*10^5
beta <- matrix(NA, nrow = num_iter, ncol = p)
sigma <- rep(1, num_iter)
marginal_prob <- rep(0, p)

### initialization
beta[1, ] <- lm(Y ~ as.matrix(X[, ], nrow = n) - 1)$coef
sigma[1] <- update_sigma(beta[1, ])

start_time <- Sys.time()
for(iter in 2:num_iter){
  ### condition on sigma, update beta
  beta[iter, ] <- BD_RJMCMC(beta[iter - 1, ], sigma[iter - 1])
  
  ### condition on beta, update sigma
  sigma[iter] <- update_sigma(beta[iter, ])
  
  ### burn in the first 10%
  if(iter > 0.1*num_iter){
    selected_variable_index <- which(!is.na(beta[iter, ]))
    marginal_prob[selected_variable_index] <- marginal_prob[selected_variable_index] + 1
  }
}
end_time <- Sys.time()
marginal_prob <- marginal_prob/sum(marginal_prob)
running_time <- as.numeric(end_time - start_time)

