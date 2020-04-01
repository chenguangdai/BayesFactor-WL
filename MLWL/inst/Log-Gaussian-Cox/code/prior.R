### prior distribution
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
prior$logdensity <- function(x) return(dmvnorm_cholesky_inverse(matrix(x, nrow = 1), prior_mean, prior_precision_chol))
prior$gradlogdensity <- function(x) return(eigenMapMatMult(matrix(parameter_mu - x, nrow = 1),  prior_precision))
prior$rinit <- function(n) return(return(fast_rmvnorm(n, mean = prior_mean, covariance = prior_cov)))
