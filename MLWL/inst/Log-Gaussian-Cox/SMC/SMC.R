### SMC for the Log-Gaussian Cox process
rm(list = ls())
library(MASS)
library(debiasedhmc)
library(spatstat)
library(smcUtils)
library(mgcv)

### define the number of grids
ngrid <- 10

### configuration
source("configuration.R")

### SMC configuration
start_time <- Sys.time()
mcmc_kernel <- list()
mcmc_kernel$parameters <- list(stepsize = 0.25, nsteps = 10)
mcmc_kernel$nmcmc <- 10
nparticles <- 500
cess_criterion <- 0.9
smc_output <- run_smc(target, proposal, nparticles, mcmc_kernel = mcmc_kernel, cess_criterion)
end_time <- Sys.time()
running_time <- end_time - start_time
data_save <- list(lognormconst = smc_output$log_ratio_normconst, running_time = running_time)

