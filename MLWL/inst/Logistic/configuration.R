##### data
X <- data[,1:24]
X <- scale(X)
Y <- data[,25]
Y <- Y - 1
n <- nrow(X)
p <- ncol(X)

Xinter <- matrix(nrow = n, ncol = p*(p-1) / 2)
index <- 1
for (j in 1:(p-1)){
  for (jprime in (j+1):p){
    Xinter[,index] <- X[,j] * X[,jprime]
    index <- index + 1
  }
}
Xinter <- scale(Xinter)
X <- cbind(X, Xinter)
p <- ncol(X)
dimension <- p + 2

##### The target distribution
target_logdensity <- function(x) debiasedhmc:::logistic_logtarget_c(x, Y, X, lambda)
target_gradlogdensity <- function(x) debiasedhmc:::logistic_gradlogtarget_c(x, Y, X, lambda)

##### The surrogate distribution
load(paste(getwd(), "/data/posterior_moments_lambda_", lambda, ".RData", sep = ""))
surrogate_mean <- data_save$mean
surrogate_cov <- data_save$cov
surrogate_precision <- solve(surrogate_cov)
surrogate_precision_chol <- t(chol(surrogate_precision))
surrogate_logdensity <- function(x) return(fast_dmvnorm_chol_inverse(matrix(x, nrow = 1), surrogate_mean, surrogate_precision_chol))
surrogate_rinit <- function() return(fast_rmvnorm(n = 1, mean = surrogate_mean, covariance = surrogate_cov))

##### A HMC kernel invariant to the target distribution
stepsize <- 0.03
nsteps <- 10
### leapfrog
leapfrog <- function(x, v, gradient){
  v <- v + stepsize * gradient(x) / 2
  for (step in 1:nsteps){
    x <- x + stepsize * v
    if (step != nsteps){
      v <- v + stepsize * gradient(x)
    }
  }
  v <- v + stepsize * gradient(x) / 2
  return(list(x = x, v = v))
}
### HMC
target_kernel <- function(x){
  x_current <- x
  v_current <- matrix(rnorm(dimension), nrow = 1)
  ### leapfrog steps
  leapfrog_result <- leapfrog(x_current, v_current, target_gradlogdensity)
  v_propose <- - leapfrog_result$v
  x_propose <- leapfrog_result$x
  loglik_current <- target_logdensity(matrix(x_current, nrow = 1))
  loglik_propose <- target_logdensity(matrix(x_propose, nrow = 1))
  ### metropolis acceptance step
  logacceptratio <- loglik_propose - loglik_current + sum(v_current^2)/2 - sum(v_propose^2)/2
  if(log(runif(1)) < logacceptratio){x = x_propose}
  return(x)
}

### Directly sampling from the surrogate distribution
surrogate_kernel <- function(x) return(surrogate_rinit())
