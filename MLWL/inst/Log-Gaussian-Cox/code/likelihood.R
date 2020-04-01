### likelihood
loglikelihood <- function(x) return(coxprocess_loglikelihood(matrix(x, nrow = 1), data_counts, parameter_area))
gradloglikelihood <- function(x) return(matrix(data_counts, nrow = 1) - parameter_area * exp(matrix(x, nrow = 1)))
