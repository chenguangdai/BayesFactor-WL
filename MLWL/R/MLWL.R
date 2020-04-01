#' Marginal likelihood estimation using the Wang-Landau algorithm
#'
#' @param num_iterations        Total number of iterations
#' @param num_burnin            Number of iterations for burn-in
#' @param target_kernel         A Markov kernel invariant to the target distribution
#' @param surrogate_kernel      A Markov kernel invariant to the surrogate distribution
#' @param target_logdensity     A function evaluating the log density of the unnormalized target distribution
#' @param surrogate_logdensity  A function evaluating the log density of the normalized surrogate distribution
#' @param learning_rate         Initialization of the learning_rate
#' @param flatness_criterion    The flatness criterion to adaptly decrease the learning rate. 
#' @param x0                    An initializing sample
#' @return The log normalizing constant (marginal likelihood) of the target distribution.
#' @export
MLWL <- function(num_iterations, num_burnin, target_kernel, surrogate_kernel, target_logdensity, 
                 surrogate_logdensity, learning_rate, flatness_criterion, x0, suppress = False){
  ### Initialization of the Wang-Landau algorithm
  target_logZ <- surrogate_logZ <- rep(0, num_iterations)
  target_hist <- surrogate_hist <- 0
  dimension <- length(x0)
  x <- matrix(0, nrow = num_iterations, ncol = dimension)
  x[1, ] <- x0
  progress_index <- 1
  
  ### The Wang-Landau algorithm
  for(iter in 2:num_iterations){
    ### Sample an indicator given the current x.
    logtarget <- target_logdensity(x[iter - 1, ]) - target_logZ[iter - 1]
    logsurrogate <- surrogate_logdensity(x[iter - 1, ]) - surrogate_logZ[iter - 1]
    maxlw <- max(logtarget, logsurrogate)
    ratio <- exp(logtarget - maxlw)/sum(exp(logtarget - maxlw) + exp(logsurrogate - maxlw))
    indicator <- ifelse(runif(1) < ratio, 1, 2)
    
    ### sample x give the current indicator
    if(indicator == 1){
      ### Local move around the target distribution
      x[iter, ] <- target_kernel(x[iter - 1, ])
    }else{
      ### Local move around the surrogate distribution
      x[iter, ] <- surrogate_kernel(x[iter - 1, ])
    }
    
    ### The Wang-Landau weight adjustment
    if(indicator == 1){
      target_logZ[iter] <- target_logZ[iter - 1] + log(1 + learning_rate)
      surrogate_logZ[iter] <- surrogate_logZ[iter - 1]
      target_hist <- target_hist + 1
    }else{
      target_logZ[iter] <- target_logZ[iter - 1]
      surrogate_logZ[iter] <- surrogate_logZ[iter - 1] + log(1 + learning_rate)
      surrogate_hist <- surrogate_hist + 1
    }
    normalize_result <- normalization(target_logZ[iter], surrogate_logZ[iter])
    target_logZ[iter] <- normalize_result$target_logZ
    surrogate_logZ[iter] <- normalize_result$surrogate_logZ
    
    ### Decrease the learning rate if the flatness criterion is satisfied
    flat_discrepancy <- abs(c(target_hist, surrogate_hist)/sum(target_hist + surrogate_hist) - rep(0.5, 2))
    if(max(flat_discrepancy) < (flatness_criterion / 2)){
      learning_rate <- 1/(1 + 1/learning_rate)
      target_hist <- surrogate_hist <- 0
    }
    
    if(iter == (progress_index/10*num_iterations) && !suppress){
      print(paste("Finish", progress_index * 10, "%! ", "The current learning rate is", round(learning_rate, 3)))
      progress_index <- progress_index + 1
    }
  }
  
  ### Return
  log_marginal_likelihood <- mean(target_logZ[-(1:num_burnin)] - surrogate_logZ[-(1:num_burnin)])
  return(list(log_marginal_likelihood = log_marginal_likelihood))
}


