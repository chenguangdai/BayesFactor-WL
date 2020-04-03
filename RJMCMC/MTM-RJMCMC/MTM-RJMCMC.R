### MTM-RJMCMC
MTM_RJMCMC <- function(beta, sigma){
  beta_current <- beta_forward <- beta
  variable_index_current <- which(!is.na(beta))
  num_variables_current <- length(variable_index_current)
 
  ### choose to move within model or across models
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
      ### move to higher dimension
      ### select a variable to add in
      add_index <- sample(x = setdiff(1:p, variable_index_current), size = 1, replace = F)
      variable_index_forward <- c(variable_index_current, add_index)
      beta_forward[add_index] <- auxiliary$rinit(1)
      
      ### calculate the moving direction (between two modes)
      mode_forward <- findmode(variable_index_forward)
      mode_current <- findmode(variable_index_current)
      mode_current[add_index] <- 0
      e <- mode_forward - mode_current
      
      ### MTM move
      MTM_result <- MTM(x = beta_forward, e = e, sigma = sigma, moveup_index = moveup_index, variable_index = add_index)
      if(num_variables_current == 1) MTM_result$accept_prob <- MTM_result$accept_prob * 0.5
      if(log(runif(1)) < log(MTM_result$accept_prob)){
        beta_current <- MTM_result$x
      } 
    }else{
      ### move to lower dimension
      ### select a variable to remove
      remove_index <- sample(x = variable_index_current, size = 1, replace = F)
      variable_index_forward <- setdiff(variable_index_current, remove_index)
      
      ### calculate the moving direction (between two modes)
      mode_forward <- findmode(variable_index_forward)
      mode_current <- findmode(variable_index_current)
      mode_forward[remove_index] <- 0
      e <- mode_forward - mode_current
      
      ### MTM move
      MTM_result <- MTM(x = beta_current, e = e, sigma = sigma, moveup_index = moveup_index, variable_index = remove_index)
      if(num_variables_current == p) MTM_result$accept_prob <- MTM_result$accept_prob * 0.5
      if(log(runif(1)) < log(MTM_result$accept_prob)){
        beta_current <- MTM_result$x
      } 
    }
  }
  return(beta_current)
}

