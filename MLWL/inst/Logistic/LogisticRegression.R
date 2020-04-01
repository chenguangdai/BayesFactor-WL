##### Normalizing constant estimation for Logistic regression
rm(list = ls())
library(MLWL)
library(debiasedhmc)

##### set a correct path
setwd(paste(getwd(), '/inst/Logistic', sep = ""))

##### load in data
data <- read.table(paste(getwd(), "/data/germancredit.txt", sep = ""), header = F, sep = '')
lambda <- 1

##### configuration
source("configuration.R")
num_iterations <- 5000
num_burnin <- num_iterations/2
learning_rate <- 1
flatness_criterion <- 0.2
x0 <- surrogate_mean

##### run the Wang-Landau algorithm
start_time <- Sys.time()
log_marginal_likelihood <- MLWL(num_iterations, num_burnin, target_kernel, surrogate_kernel,
                                target_logdensity, surrogate_logdensity, learning_rate, flatness_criterion, x0)
end_time <- Sys.time()
running_time <- end_time - start_time
print(paste("The log marginal likelihood is", log_marginal_likelihood))
print(paste("The running time is", running_time))

