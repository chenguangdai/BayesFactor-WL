### conditioning on sigma, update beta
MH_gibbs <- function(beta, sigma){
  beta_current <- beta
  variable_index <- which(!is.na(beta))
  logposterior_current <- logposterior(beta_current, sigma)
  for(j in variable_index){
    beta_forward <- beta_current
    beta_forward[j] <- rnorm(1, mean = beta_current[j], sd = sigma_proposal)
    logposterior_forward <- logposterior(beta_forward, sigma)
    if(log(runif(1)) < (logposterior_forward - logposterior_current)){
      beta_current <- beta_forward
      logposterior_current <- logposterior_forward
    }
  }
  return(beta_current)
}

### conditioning on beta, update sigma
update_sigma <- function(beta){
  index <- which(!is.na(beta))
  fittedvalue <- as.matrix(X[, index], nrow = n)%*%beta[index]
  a <- (n + length(index))/2
  b <- 0.5*((1/g)*sum(fittedvalue^2) + sum((Y - fittedvalue)^2))
  sigma <- sqrt(rigamma(1, alpha = a, beta = b))
  return(sigma)
}
