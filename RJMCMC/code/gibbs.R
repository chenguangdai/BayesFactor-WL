### conditioning on sigma, update beta
update_beta <- function(beta, sigma){
  index <- which(!is.na(beta))
  Xgamma <- as.matrix(X[, index], nrow = n)
  XtXinv <- solve(t(Xgamma)%*%Xgamma)
  XtXinv <- 0.5*(XtXinv + t(XtXinv))
  XtY <- t(Xgamma)%*%Y
  beta[index] <- rmvn(n = 1, mu = as.vector(g/(g + 1)*XtXinv%*%XtY), Sigma = g/(g + 1)*sigma^2*XtXinv)
  return(beta)
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
