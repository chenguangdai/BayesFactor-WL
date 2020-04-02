setwd('..')

### load in data and set up parameters
source("code/load_data.R")

### load in C functions.
sourceCpp("code/mvnorm.cpp", verbose = F)

### prior distribution
source("code/prior.R")
prior$logdensity <- function(x) return(dmvnorm_cholesky_inverse(x, prior_mean, prior_precision_chol))

### likelihood
source("code/likelihood.R")

### target distribution
target <- list()
target$logdensity <- function(x) return(prior$logdensity(x) + loglikelihood_matrix(x))
target$gradlogdensity <- function(x) return(prior$gradlogdensity(x) + gradloglikelihood_matrix(x))

### set prior as proposal
proposal <- prior

setwd('./SMC')

### source the code
source("adaptive_smc.R")
source("hmc.R")
source("compute_cess.R")
source("find_next_temp.R")

