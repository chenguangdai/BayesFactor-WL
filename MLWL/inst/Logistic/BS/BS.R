### bridge sampling for Logistic regression
rm(list = ls())
library(rstan)
library(bridgesampling)

### stan model
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

### load in data
setwd('..')
source("load_data.R")
setwd('./BS')

### lambda
lambda <- 1

### run stan
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
lognormconst = bridge_result$logml
 
