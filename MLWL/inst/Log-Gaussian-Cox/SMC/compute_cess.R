compute_cess <- function(gamma_current, gamma_proposal, cess_criterion, logweights_IS, weights_current){
  ### calculate the incremental weights
  logweights_increment <- (gamma_proposal - gamma_current) * logweights_IS
  maxlogweights_increment <- max(logweights_increment)
  weights_increment <- exp(logweights_increment - maxlogweights_increment)
  
  ### calculate the conditional effective sample size
  weights_sum <- sum(weights_current * weights_increment)
  cess <- weights_sum^2/sum(weights_current * weights_increment^2)
  criterion_condition <- FALSE
  
  if((gamma_proposal == 1 & cess > cess_criterion) | abs(cess - cess_criterion) < 10^(-2)){
    ### calculate the normalized weights
    normweights <- weights_current * weights_increment/weights_sum
    
    ### calculate the increment of the log normalizing constant estimate
    log_ratio_normconst <- log(weights_sum) + maxlogweights_increment
    criterion_condition <- TRUE
    return(list(cess = cess, log_ratio_normconst = log_ratio_normconst, normweights = normweights, criterion_condition = criterion_condition))
  }else{
    return(list(cess = cess, criterion_condition = criterion_condition))
  }
}

