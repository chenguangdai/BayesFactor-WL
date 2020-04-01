setwd('..')
### load in C functions.
sourceCpp("code/mvnorm.cpp", verbose = F)

### load in data and set up parameters
source("code/load_data.R")

### prior distribution
source("code/prior.R")

### likelihood
source("code/likelihood.R")

### HMC
source("code/HMC.R")
