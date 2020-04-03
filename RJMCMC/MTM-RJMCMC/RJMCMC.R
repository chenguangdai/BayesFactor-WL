### augmented logposterior
logaugtarget <- function(beta, sigma, variable_index){
  beta_temp <- beta
  beta_temp[variable_index] <- NA
  return(logposterior(beta_temp, sigma) + logproposal(beta[variable_index]))
}
