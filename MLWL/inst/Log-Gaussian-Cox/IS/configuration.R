# parameters
data(finpines)
data_x <- (finpines$x + 5) / 10 # normalize data to unit square
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

# target
likelihood <- list()
likelihood$log <- function(x) return(coxprocess_loglikelihood(x, data_counts, parameter_area))
targetlogdensity <- function(x) return(prior$logdensity(x) + likelihood$log(x))


##### surrogate distribution
# load(paste("/Users/chenguang/Desktop/working papers/MTM/revision/code/log-cox-process/WL/surrogate_moments/ngrid_", 
#            ngrid, ".RData", sep = ""))
load(paste("/Users/chenguang/Desktop/working papers/MTM/revision/code/log-cox-process/IS/mode/mode", 
           ngrid, ".RData", sep = ""))
surrogate_mean <- data_save$surrogate_mu
surrogate_cov <- matrix(0, nrow = dimension, ncol = dimension)
# diag(surrogate_cov) <- mean(diag(surrogate_moments$surrogate_cov))
if(ngrid == 10){
  diag(surrogate_cov) <- 1.0^2
}else if(ngrid == 20){
  diag(surrogate_cov) <- 1.2^2
}else{
  diag(surrogate_cov) <- 1.3^2
}
surrogate_precision <- solve(surrogate_cov)
surrogate_precision_chol <- t(chol(surrogate_precision))
surrogatelogdensity <- function(x) return(dmvnorm_cholesky_inverse(x, surrogate_mean, surrogate_precision_chol))
surrogaterinit <- function(n){
  xparticles <- matrix(0, nrow = n, ncol = dimension)
  sigma <- sqrt(diag(surrogate_cov))
  for(j in 1:dimension){
    xparticles[, j] <- rnorm(n, mean = surrogate_mean[j], sd = sigma[j])
  }
  return(xparticles)
} 
