setwd('..')

### load in data
load(paste("data/XY/BayesianLasso_SNR_", SNR, ".RData", sep = ""))
n <- length(y)
p <- ncol(X)
dimension <- 2*p + 1

### pre-calculate quantity
XtX <- t(X)%*%X
Xty <- c(t(X)%*%y)
yty <- sum(y^2)

### prior
prior <- list()
prior$logdensity <- function(x){
  beta <- x[, 1:p]
  eta <- x[, (p + 1):(2*p)]
  xi <- x[, 2*p + 1]
  ### log density of beta
  transformed_beta <- sweep(beta/exp(eta/2), 1, exp(xi/2), '/')
  logjacobian <- -p*xi/2 - rowSums(eta)/2
  logdensity_beta <- dmvnorm_cholesky_inverse(transformed_beta, rep(0, p), diag(p)) + logjacobian
  ### log density of eta
  logdensity_eta <- rowSums(2*log(lambda) - log(2) - lambda^2/2*exp(eta) + eta)
  ### log density of xi
  logdensity_xi <- 0
  return(logdensity_beta + logdensity_eta + logdensity_xi)
}

### likelihood
loglikelihood <- function(x){
  num_samples <- nrow(x)
  beta <- x[, 1:p]
  eta <- x[, (p + 1):(2*p)]
  xi <- x[, 2*p + 1]
  part1 <- -n/2*log(2*pi) - n/2*xi 
  part2 <- -0.5*(colSums((matrix(rep(y, num_samples), ncol = num_samples) - X%*%t(beta))^2)*exp(-xi))
  return(part1 + part2)
} 

### the target distribution
target <- list()
target$logdensity <- function(x) return(prior$logdensity(x) + loglikelihood(x))

### the proposal distribution
load(paste("data/VI/VI_SNR_", SNR, "_lambda_", lambda, ".RData", sep = ""))
proposal_mean <- c(data_save$surrogate_m, data_save$surrogate_phi, data_save$surrogate_u)
proposal_cov <- diag(c(data_save$surrogate_ssq, data_save$surrogate_zetasq, data_save$surrogate_vsq))
proposal <- list()
proposal_precision <- solve(proposal_cov)
proposal_precision_chol <- t(chol(proposal_precision))
proposal$logdensity <- function(x) return(dmvnorm_cholesky_inverse(x, proposal_mean, proposal_precision_chol))
proposal$rinit <- function(n) return(fast_rmvnorm(n = n, mean = proposal_mean, covariance = proposal_cov))

setwd('./IS')

### source C functions
sourceCpp("mvnorm.cpp")


