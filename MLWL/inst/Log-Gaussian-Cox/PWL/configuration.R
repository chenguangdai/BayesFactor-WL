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

### load the temeprature sequence
load(paste("sequence/", ngrid, ".RData", sep = ""))
temperature <- data_save$sequence

### load the mode
load(paste("mode/", ngrid, ".RData", sep = ""))
x0 <- data_save$surrogate_mu

