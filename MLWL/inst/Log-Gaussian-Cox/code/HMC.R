### leapfrog
stepsize <- 0.25
nsteps <- 10
leapfrog <- function(x, v, gradient){
  v <- v + stepsize * gradient(x) / 2
  for (step in 1:nsteps){
    x <- x + stepsize * v
    if (step != nsteps){
      v <- v + stepsize * gradient(x)
    }
  }
  v <- v + stepsize * gradient(x) / 2
  return(list(x = x, v = v))
}

### HMC
construct_kernel <- function(target){
  logdensity <- target$logdensity
  gradlogdensity <- target$gradlogdensity
  
  kernel <- function(x){
    x_current <- x
    v_current <- matrix(rnorm(ngrid^2), nrow = 1)
    ### leapfrog steps
    leapfrog_result <- leapfrog(x_current, v_current, gradlogdensity)
    v_propose <- - leapfrog_result$v
    x_propose <- leapfrog_result$x
    loglik_current <- logdensity(x_current)
    loglik_propose <- logdensity(x_propose)
    ### metropolis acceptance step
    logacceptratio <- loglik_propose - loglik_current + sum(v_current^2)/2 - sum(v_propose^2)/2
    if(log(runif(1)) < logacceptratio){x <- x_propose}
    return(x)
  }
  
  return(kernel)
}
