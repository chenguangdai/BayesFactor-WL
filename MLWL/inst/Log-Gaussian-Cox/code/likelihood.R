### likelihood
loglikelihood <- function(x) return(coxprocess_loglikelihood(matrix(x, nrow = 1), data_counts, parameter_area))
gradloglikelihood <- function(x) return(matrix(data_counts, nrow = 1) - parameter_area * exp(matrix(x, nrow = 1)))

loglikelihood_matrix <- function(x) return(coxprocess_loglikelihood(x, data_counts, parameter_area))
gradloglikelihood_matrix <- function(x) return(matrix(rep(data_counts, nrow(x)), nrow = nrow(x), byrow = T) - parameter_area * exp(x))

