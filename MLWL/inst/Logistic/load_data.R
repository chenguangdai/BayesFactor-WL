### load in data
data <- read.table("data/germancredit.txt", header = F)
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
