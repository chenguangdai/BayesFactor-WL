### birth and death RJMCMC
BD_RJMCMC <- function(beta, sigma){
  beta_current <- beta_forward <- beta
  variable_index_current <- which(!is.na(beta))
  num_variables_current <- length(variable_index_current)

  ### decide to move within model or across models
  move_across_model <- rbern(1, prob = 0.5)
  if(move_across_model == 0){
    ### stay within the same dimension
    ### gibbs move
    beta_current <- update_beta(beta_current, sigma)
  }else{
    ### move across different dimensions
    moveup_index <- rbern(1, prob = 0.5)
    ### handle boundary cases
    if(num_variables_current == 1) moveup_index <- 1
    if(num_variables_current == p) moveup_index <- 0
    
    if(moveup_index == 1){
      ### move to a higher dimension
      ### select a variable to add in
      add_index <- sample(x = setdiff(1:p, variable_index_current), size = 1, replace = F)
      variable_index_forward <- c(variable_index_current, add_index)
      beta_forward[add_index] <- rproposal(1)
    }else{
      ### move to a lower dimension
      ### select a variable to remove
      remove_index <- sample(x = variable_index_current, size = 1, replace = F)
      variable_index_forward <- setdiff(variable_index_current, remove_index)
      beta_forward[remove_index] <- NA
    }
    
    ### calculate the acceptance probability
    logposterior_forward <- logposterior(beta_forward, sigma)
    logposterior_current <- logposterior(beta_current, sigma)
    logacceptance_prob <- logposterior_forward - logposterior_current
    if(num_variables_current == 1) logacceptance_prob <- logacceptance_prob - log(2)
    if(num_variables_current == p) logacceptance_prob <- logacceptance_prob - log(2)
    if(log(runif(1)) < logacceptance_prob){
      beta_current <- beta_forward
    }
  }
  return(beta_current)
}
