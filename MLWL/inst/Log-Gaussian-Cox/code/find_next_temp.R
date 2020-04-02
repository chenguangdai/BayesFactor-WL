find_next_temp <- function(gamma_current, cess_criterion, logweights_IS, weights_current){
  ### initialize the proposed inverse temperature to be 1
  gamma_proposal <- 1
  
  ### compute the conditional effective sample size
  cess_result <- compute_cess(gamma_current, gamma_proposal, cess_criterion, logweights_IS, weights_current)
  
  if(cess_result$criterion_condition){
    return(list(cess = cess_result$ess, gamma_next = gamma_proposal, log_ratio_normconst = cess_result$log_ratio_normconst,
                normweights = cess_result$normweights))
  }else{
    ### update the proposed inverse temperature to be the middle point
    gamma_lower <- gamma_current
    gamma_upper <- 1
    gamma_proposal <- (gamma_lower + gamma_upper) / 2
    
    ### calculate the conditional effective sample size
    cess_result <- compute_cess(gamma_current, gamma_proposal, cess_criterion, logweights_IS, weights_current)
  }

  while(!cess_result$criterion_condition){
    cess_current <- cess_result$cess
    ### if the conditional effective sample size is larger than the given threshold
    if(cess_current > cess_criterion){
      ### increase the proposed inverse temperature by moving up the lower bound
      gamma_lower <- gamma_proposal
    }else{
      ### otherwise decrease the proposed inverse temperature by moving down the upper bound
      gamma_upper <- gamma_proposal
    }
    ### update the proposed inverse temperature to be the middle point
    gamma_proposal <- (gamma_lower + gamma_upper) / 2
    
    ### calculate the conditional effective sample size
    cess_result <- compute_cess(gamma_current, gamma_proposal, cess_criterion, logweights_IS, weights_current)
  }

  return(list(cess = cess_result$cess, gamma_next = gamma_proposal, log_ratio_normconst = cess_result$log_ratio_normconst,
              normweights = cess_result$normweights))
}



