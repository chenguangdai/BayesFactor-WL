setwd('..')

### load in data
data <- read.table("data/pollution.txt", header = T, sep = "")
Y <- data[, 16]
Y <- scale(Y)
n <- length(Y)
X <- apply(data[, 1:15], 2, scale)
p <- 15
g <- exp(15)

### mode finding 
findmode <- function(variable_index){
  mode <- rep(NA, p)
  mode[variable_index] <- .lm.fit(as.matrix(X[, variable_index], nrow = n), Y)$coefficients
  return(mode)
}

### auxiliary distribution of beta
auxiliary <- list()
auxiliary$rinit <- function(n) return(rnorm(n, mean = 0, sd = sigma_auxiliary))
auxiliary$logdensity <- function(beta){
  return(sum(dnorm(beta, mean = 0, sd = sigma_auxiliary, log = TRUE)))
} 

### log posterior of beta conditioning on sigma
logposterior <- function(beta, sigma){
  index <- which(!is.na(beta))
  qgamma <- length(index)
  Xgamma <- as.matrix(X[, index], nrow = n)
  XtX <- t(Xgamma)%*%Xgamma
  fittedvalue <- Xgamma%*%beta[index]
  logp <- -qgamma/2*log(g) + 0.5*log(det(XtX)) - qgamma/2*log(sigma^2) - 0.5/sigma^2*((g + 1)/g*sum(fittedvalue^2) - 2*sum(fittedvalue*Y))
  return(logp)
}
