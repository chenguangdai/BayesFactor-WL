### model configuration
library(MASS)
library(Rcpp)
library(RcppEigen)
library(spatstat)
library(smcUtils)
library(mgcv)

### source the code
source("smc.R")
source("hmc.R")
source("compute_cess.R")
source("find_next_temp.R")
sourceCpp("coxprocess.cpp")

### model parameters
data(finpines)
data_x <- (finpines$x + 5) / 10 
data_y <- (finpines$y + 8) / 10
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

### prior distribution
prior <- list()
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
prior$logdensity <- function(x) return(dmvnorm_cholesky_inverse(x, prior_mean, prior_precision_chol))
prior$gradlogdensity <- function(x) return(eigenMapMatMult(parameter_mu - x,  prior_precision))
prior$rinit <- function(n) return(return(rmvn(n = n, mu = prior_mean, V = prior_cov)))

### likelihood
loglikelihood <- function(x) return(coxprocess_loglikelihood(x, data_counts, parameter_area))
gradloglikelihood <- function(x) return(matrix(rep(data_counts, nparticles), nrow = nparticles, byrow = T) - parameter_area * exp(x))

### target distribution
target <- list()
target$logdensity <- function(x) return(prior$logdensity(x) + loglikelihood(x))
target$gradlogdensity <- function(x) return(prior$gradlogdensity(x) + gradloglikelihood(x))

### using prior as proposal
proposal <- prior

### using independent normal as proposal
# load(paste(getwd(), "/surrogate_mean/", ngrid, ".RData", sep = ""))
# proposal_mean <- data_save
# if(ngrid == 10){
#   proposal_sigma <- rep(1.0, dimension)
# }else if(ngrid == 20){
#   proposal_sigma <- rep(1.2, dimension)
# }else{
#   proposal_sigma <- rep(1.3, dimension)
# }
# proposal_precision <- solve(diag(proposal_sigma))
# proposal_precision_chol <- t(chol(proposal_precision))
# proposal <- list()
# proposal$logdensity <- function(x) return(dmvnorm_cholesky_inverse(x, proposal_mean, proposal_precision_chol))
# proposal$gradlogdensity <- function(x) return(eigenMapMatMult(parameter_mu - x,  proposal_precision))
# proposal$rinit <- function(n) return(return(rmvn(n = n, mu = proposal_mean, V = diag(proposal_sigma))))





