### MTM move 
MTM <- function(x, e, sigma, moveup_index, variable_index){
  ### moveup_index = 1 represents adding one variable
  ### moveup_index = 0 represents removing one variable
  ### variable_index represents the variable added or removed
  
  ### propose multiple tries
  r <- rnorm(m, mean = 1, sd = sigma_r)
  y <- outer(rep(1, m), x) + outer(r, e)  
  
  ### calculate the weight function
  logtarget_y <- rep(0, m)
  for(j in 1:m){
    if(moveup_index == 1) logtarget_y[j] <- logposterior(y[j, ], sigma)
    if(moveup_index == 0) logtarget_y[j] <- logaugtarget(y[j, ], sigma, variable_index)
  }
  maxlogtarget_y <- max(logtarget_y)
  if(maxlogtarget_y == -Inf) return(list(x = x, accept_prob = 0))
  logsumtarget_y <- maxlogtarget_y + log(sum(exp(logtarget_y - maxlogtarget_y)))
  
  ### sample y_star
  index <- which(rmultinom(1, 1, prob = exp(logtarget_y - maxlogtarget_y)) == 1)
  ystar <- y[index, ]
  
  ### reverse sampling
  xrev <- outer(rep(1, m), ystar) - outer(r, e)
  logtarget_xrev <- rep(0, m)
  for(j in 1:m){
    if(moveup_index == 1) logtarget_xrev[j] <- logaugtarget(xrev[j, ], sigma, variable_index)
    if(moveup_index == 0) logtarget_xrev[j] <- logposterior(xrev[j, ], sigma)
  } 
  maxlogtarget_xrev <- max(logtarget_xrev)
  logsumtarget_xrev <- maxlogtarget_xrev + log(sum(exp(logtarget_xrev - maxlogtarget_xrev)))
  
  ### calculate the acceptance probability
  accept_prob <- min(1, exp(logsumtarget_y - logsumtarget_xrev))
  if(moveup_index == 0) ystar[variable_index] <- NA
  return(list(x = ystar, accept_prob = accept_prob))
}
