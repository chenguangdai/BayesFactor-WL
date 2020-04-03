### MTM move 
MTM <- function(x, e, sigma, moveup_index, variable_index){
  ### both x and e and length p where non-selected variables are indicated by NA
  ### moveup_index represents the direction of moving
  ### variable_index represents the variable added or removed
  r <- rnorm(m, mean = 1, sd = sigma_r)
  y <- outer(rep(1, m), x) + outer(r, e)   ### Y denote data and y denote proposal
  logtarget_y <- rep(0, m)
  for(j in 1:m){
    if(moveup_index == 1) logtarget_y[j] <- logposterior(y[j, ], sigma)
    if(moveup_index == 0) logtarget_y[j] <- logaugtarget(y[j, ], sigma, variable_index)
  }
  maxlogtarget_y <- max(logtarget_y)
  if(maxlogtarget_y == -Inf) return(list(x = x, accept_prob = 0))
  logsumtarget_y <- maxlogtarget_y + log(sum(exp(logtarget_y - maxlogtarget_y)))
  index <- which(rmultinom(1, 1, prob = exp(logtarget_y - maxlogtarget_y)) == 1)
  ystar <- y[index, ]
  xrev <- outer(rep(1, m), ystar) - outer(r, e)
  logtarget_xrev <- rep(0, m)
  for(j in 1:m){
    if(moveup_index == 1) logtarget_xrev[j] <- logaugtarget(xrev[j, ], sigma, variable_index)
    if(moveup_index == 0) logtarget_xrev[j] <- logposterior(xrev[j, ], sigma)
  } 
  maxlogtarget_xrev <- max(logtarget_xrev)
  logsumtarget_xrev <- maxlogtarget_xrev + log(sum(exp(logtarget_xrev - maxlogtarget_xrev)))
  accept_prob <- min(1, exp(logsumtarget_y - logsumtarget_xrev))
  if(moveup_index == 0) ystar[variable_index] <- NA
  return(list(x = ystar, accept_prob = accept_prob))
}
