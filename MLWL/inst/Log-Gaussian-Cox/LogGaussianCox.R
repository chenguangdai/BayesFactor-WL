##### Normalizing constant estimation for the Log-Gaussian-Cox process
rm(list = ls())
library(MLWL)
library(debiasedhmc)
library(spatstat)

##### setup the correct poth
setwd(paste(getwd(), '/inst/Log-Gaussian-Cox', sep = ""))

##### load in data and set up parameters
data(finpines)
data_x <- (finpines$x + 5) / 10
data_y <- (finpines$y + 8) / 10
ngrid <- 30
grid <- seq(from = 0, to = 1, length.out = ngrid + 1)
dimension <- ngrid^2
data_counts <- rep(0, dimension)
for (i in 1:ngrid){
  for (j in 1:ngrid){
    logical_y <- (data_x > grid[i]) * (data_x < grid[i+1])
    logical_x <- (data_y > grid[j]) * (data_y < grid[j+1])
    data_counts[(i-1)*ngrid + j] <- sum(logical_y * logical_x)
  }
}
parameter_sigmasq <- 1.91
parameter_beta <- 1/33
parameter_area <- 1 / dimension
parameter_mu <- log(126) - 0.5 * parameter_sigmasq

##### configuration
source("configuration.R")
num_iterations <- 10^5
num_burnin <- num_iterations/2
learning_rate <- 1
flatness_criterion <- 0.1
x0 <- surrogate_mu

##### run the Wang-Landau algorithm
start_time <- Sys.time()
log_marginal_likelihood <- MLWL(num_iterations, num_burnin, target_kernel, surrogate_kernel,
                                target_logdensity, surrogate_logdensity, learning_rate, flatness_criterion, x0)
end_time <- Sys.time()
running_time <- end_time - start_time
print(paste("The log marginal likelihood is", log_marginal_likelihood))
print(paste("The running time is", running_time))





