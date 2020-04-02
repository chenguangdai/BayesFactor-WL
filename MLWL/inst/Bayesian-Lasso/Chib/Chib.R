### Chib's method for Bayesian lasso
rm(list = ls())
library(MASS)
library(statmod)
library(geoR)
library(mvtnorm)

### SNR and lambda
SNR <- 3
lambda <- 5

### load in data
setwd('..')
load(paste("data/XY/BayesianLasso_SNR_", SNR, ".RData", sep = ""))
setwd('./Chib')

### pre-calculate quantity
n <- length(y)
p <- ncol(X)
dimension <- 2*p + 1
XtX <- t(X)%*%X
Xty <- c(t(X)%*%y)
yty <- sum(y^2)

### Gibbs sampler
gibbs_sampler <- function(beta, tausq, sigmasq){  
  ### given the rest, update beta
  A_tau_inverse <- solve(XtX + diag(tausq^(-1)))
  beta <- as.vector(rmvnorm(1, mean = A_tau_inverse%*%Xty, sigma = sigmasq * A_tau_inverse))
  ### given the rest, update tausq
  tausq <- (rinvgauss(p, lambda * sqrt(sigmasq)/abs(beta), lambda^2))^(-1)
  ### given the rest, update sigmasq
  sigmasq <- rinvchisq(1, n + p, scale = (sum((y - X%*%beta)^2) + sum(beta^2/tausq))/(n + p))
  return(list(beta = beta, tausq = tausq, sigmasq = sigmasq))
}

### get full posterior samples
start_time <- Sys.time()
num_iter <- 5000
beta <- matrix(0, nrow = num_iter, ncol = p)
tausq <- matrix(1, nrow = num_iter, ncol = p)
sigmasq <- rep(1, num_iter)
for(iter in 2:num_iter){
  gibbs_result <- gibbs_sampler(beta[iter - 1, ], tausq[iter - 1, ], sigmasq[iter - 1])
  beta[iter, ] <- gibbs_result$beta
  tausq[iter, ] <- gibbs_result$tausq
  sigmasq[iter] <- gibbs_result$sigmasq
}

### burn in the first 10%
burnin <- c(1:(num_iter*0.1))
beta <- beta[-burnin, ]
tausq <- tausq[-burnin, ]
sigmasq <- sigmasq[-burnin]

### estimate the (normalized) posterior density
beta_star <- colMeans(beta)
tausq_star <- colMeans(tausq)
sigmasq_star <- mean(sigmasq)

### step 1: estimate p(beta_star|y)
num_samples <- nrow(beta)
logdensity <- rep(0, num_samples)
for(iter in 1:num_samples){
  A_tau_inverse <- solve(XtX + diag(tausq[iter, ]^(-1)))
  logdensity[iter] <- dmvnorm(beta_star, mean = A_tau_inverse%*%Xty, sigma = sigmasq[iter] * A_tau_inverse, log = T)
}
maxlw <- max(logdensity)
logpart1 <- log(mean(exp(logdensity - maxlw))) + maxlw

### step 2: estimate p(tausq_star|y, beta_star)
### simulate sigmasq given y and beta_star
sigmasq_temp <- rep(1, num_iter)
tausq_temp <- matrix(1, nrow = num_iter, ncol = p)
for(iter in 2:num_iter){
  ### given the rest, update tausq
  tausq_temp[iter, ] <- (rinvgauss(p, lambda * sqrt(sigmasq_temp[iter - 1])/abs(beta_star), lambda^2))^(-1)
  ### given the rest, update sigmasq
  sigmasq_temp[iter] <- rinvchisq(1, n + p, scale = (sum((y - X%*%beta_star)^2) + sum(beta_star^2/tausq_temp[iter, ]))/(n + p))
}
tausq_temp <- tausq_temp[-c(1:(num_iter)*0.1), ]
sigmasq_temp <- sigmasq_temp[-c(1:(num_iter)*0.1)]

### estimate p(tausq_star|y, beta_star)
logdensity <- rep(0, num_samples)
for(iter in 1:num_samples){
  logdensity[iter] <- sum(dinvgauss(1/tausq_star, lambda * sqrt(sigmasq_temp[iter])/abs(beta_star), lambda^2, log = T)) - 2*sum(log(tausq_star))
}
maxlw <- max(logdensity)
logpart2 <- log(mean(exp(logdensity - maxlw))) + maxlw

### step 3: estimate p(sigmasq_star|beta_star, tausq_star, y)
nu <- n + p
tau0sq <- (sum((y - X%*%beta_star)^2) + sum(beta_star^2/tausq_star))/(n + p)
logpart3 <- nu/2*log(tau0sq*nu/2) - lgamma(nu/2) - nu*tau0sq/(2*sigmasq_star) - (1 + nu/2)*log(sigmasq_star)

### log normalizing constant
logposterior <- logpart1 + logpart2 + logpart3
logprior <- -log(sigmasq_star) + sum(dexp(tausq_star, rate = lambda^2/2, log = T)) + dmvnorm(beta_star, mean = rep(0, p), sigma = sigmasq_star*diag(tausq_star), log = T)
loglikelihood <- sum(dnorm(y, mean = X%*%beta_star, sd = sqrt(sigmasq_star), log = T))
lognormconst <- loglikelihood + logprior - logposterior
end_time <- Sys.time()
running_time <- end_time - start_time
