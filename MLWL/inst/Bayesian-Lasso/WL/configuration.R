setwd('..')

### load in data
load("data/XY/BayesianLasso_SNR_", SNR, ".RData", sep = ""))

### pre-calculate quantity
XtX <- t(X)%*%X
Xty <- c(t(X)%*%y)
yty <- sum(y^2)
n <- length(y)
p <- ncol(X)
dimension <- 2*p + 1

### the target distribution
target <- list()
target$logdensity <- function(x){
  beta <- x[1:p]
  eta <- x[(p + 1):(2*p)]
  xi <- x[2*p + 1]
  part1 <- -(n + p)/2*log(2*pi) + (2*p)*log(lambda) - p*log(2) - (n + p)/2*xi + 0.5*sum(eta) - lambda^2/2*sum(exp(eta))
  part2 <- -0.5*exp(-xi)*(yty - 2*sum(beta*Xty) + t(beta)%*%XtX%*%beta + sum(beta^2 * exp(-eta)))
  return(part1 + part2)
}

### the surrogate distribution
surrogate <- list()
load("data/VI/VI_SNR_", SNR, "_lambda_", lambda, ".RData", sep = ""))
surrogate_mu <- c(data_save$surrogate_m, data_save$surrogate_phi, data_save$surrogate_u)
surrogate_sigma <- sqrt(c(data_save$surrogate_ssq, data_save$surrogate_zetasq, data_save$surrogate_vsq))
surrogate$logdensity <- function(x){
  return(-dimension/2*log(2*pi) - sum(log(surrogate_sigma)) - 0.5 * sum((x - surrogate_mu)^2/surrogate_sigma^2))
} 

### a Gibbs kernel invariant to the target distribution
target$kernel <- function(x){
  beta <- x[1:p]
  tausq <- exp(x[(p + 1):(2*p)])
  sigmasq <- exp(x[2*p + 1])
  ### given the rest, update beta
  A_tau_inverse <- solve(XtX + diag(tausq^(-1)))
  beta <- as.vector(rmvnorm(1, mean = A_tau_inverse%*%Xty, sigma = sigmasq * A_tau_inverse))
  ### given the rest, update tausq
  tausq <- (rinvgauss(p, lambda * sqrt(sigmasq)/abs(beta), lambda^2))^(-1)
  ### given the rest, update sigmasq
  sigmasq <- rinvchisq(1, n + p, scale = (sum((y - X%*%beta)^2) + sum(beta^2/tausq))/(n + p))
  return(c(beta, log(tausq), log(sigmasq)))
}

### direct sampling from the surrogate distribution
surrogate$kernel <- function(x){
  return(rnorm(dimension, mean = surrogate_mu, sd = surrogate_sigma))
}

