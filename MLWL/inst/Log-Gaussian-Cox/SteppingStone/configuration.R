##### load in C functions.
sourceCpp("mvnorm.cpp", verbose = F)

##### load in data and set up parameters
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

##### prior distribution
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
prior <- list()
prior$logdensity <- function(x) return(dmvnorm_cholesky_inverse(x, prior_mean, prior_precision_chol))
prior$gradlogdensity <- function(x) return(eigenMapMatMult(parameter_mu - x,  prior_precision))
prior$rinit <- function(n) return(return(fast_rmvnorm(n, mean = prior_mean, covariance = prior_cov)))

##### likelihood
loglikelihood <- function(x) return(coxprocess_loglikelihood(x, data_counts, parameter_area))
gradloglikelihood <- function(x) return(matrix(rep(data_counts, nrow(x)), nrow = nrow(x), byrow = T) - parameter_area * exp(x))

### leapfrog
stepsize <- 0.25
nsteps <- 10
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
construct_kernel <- function(target){
  logdensity <- target$logdensity
  gradlogdensity <- target$gradlogdensity
  
  kernel <- function(x){
    x_current <- x
    v_current <- matrix(rnorm(ngrid^2), nrow = 1)
    ### leapfrog steps
    leapfrog_result <- leapfrog(x_current, v_current, gradlogdensity)
    v_propose <- - leapfrog_result$v
    x_propose <- leapfrog_result$x
    loglik_current <- logdensity(x_current)
    loglik_propose <- logdensity(x_propose)
    ### metropolis acceptance step
    logacceptratio <- loglik_propose - loglik_current + sum(v_current^2)/2 - sum(v_propose^2)/2
    if(log(runif(1)) < logacceptratio){x <- x_propose}
    return(x)
  }
  
  return(kernel)
}
