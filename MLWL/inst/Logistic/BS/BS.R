rm(list = ls())
library(rstan)
library(bridgesampling)

##### stan model
LogisticRegression.stan <- '
data {
  int<lower = 0> n;       // number of observations
  int<lower = 0> p;       // number of variables
  real<lower = 0> lambda;
  matrix[n, p] X;         // design matrix X
  int y[n];               // observations y
}
parameters {
  real alpha;
  vector[p] beta; 
  real logssq;
}
transformed parameters {
  vector[n] X_beta;
  X_beta = X * beta + alpha;
}
model {
  // prior
  beta ~ normal(0, exp(0.5 * logssq));
  alpha ~ normal(0, exp(0.5 * logssq));
  exp(logssq) ~ exponential(lambda);
  target += logssq;

  // likelihood
  y ~ bernoulli_logit(X_beta);
}
'

##### data preparation 
germancredit <- read.table("/Users/chenguang/Desktop/working papers/MTM/revision/code/Logistic/data/germancredit.txt", header = F, sep = '')
X <- germancredit[, 1:24]
X <- scale(X)
y <- germancredit[, 25]
y <- y - 1
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
lambda <- 1

##### run Stan 
num_replicate <- 10
for(replicate in 1:num_replicate){
  warmup <- 1000
  iteration <- 3500
  stan.fit <- stan(model_code = LogisticRegression.stan, 
                   data = list(
                     n = n,
                     p = p,
                     lambda = lambda,
                     X = X,
                     y = y),
                   warmup = warmup, iter = iteration, chains = 1, control = list(adapt_delta = 0.8))
  stan.fit.extract = extract(stan.fit)
  bridge_result <- bridge_sampler(stan.fit, method = "normal", maxiter = 1000)
  print(replicate)
  data_save <- list(lognormconst = bridge_result$logml)
  save(data_save, file = paste("/Users/chenguang/Desktop/working papers/MTM/revision/code/Logistic/result/BS/", 
                               iteration - warmup, "/1.00", "/replicate_", 
                               replicate, ".RData", sep = ""))
}





