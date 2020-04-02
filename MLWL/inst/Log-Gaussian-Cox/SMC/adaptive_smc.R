### SMC sampler
run_smc <- function(target, proposal, nparticles, mcmc_kernel, cess_criterion){
  ### sample from the initial distribution
  xparticles <- proposal$rinit(nparticles)
  logtarget <- target$logdensity(xparticles)
  logproposal <- proposal$logdensity(xparticles)
  
  ### initialize the inverse temperature
  gamma_current <- 0
  temperature <- gamma_current
  dimension <- ncol(xparticles)
  weights_current <- rep(1, nparticles)/nparticles
  cess <- numeric()
  
  ### initialize the log normalizing constant estimator
  log_ratio_normconst <- 0

  while(gamma_current < 1){
    ### find the next inverse temperature
    logweights_IS <- logtarget - logproposal
    next_temp_result <- find_next_temp(gamma_current, cess_criterion, logweights_IS, weights_current)
    cess <- c(cess, next_temp_result$cess)
    
    ### update the current temperature
    gamma_current <- next_temp_result$gamma_next
    temperature <- c(temperature, gamma_current)
    
    ### update the log normalizing constant estimator
    log_ratio_normconst <- log_ratio_normconst + next_temp_result$log_ratio_normconst
    weights_current <- next_temp_result$normweights

    ### systematic resampling if ESS drops below 0.5
    ess <- 1/sum(weights_current^2)/nparticles
    if(ess < 0.5){
      ancestors <- systematic.resample(weights_current)
      xparticles <- xparticles[ancestors, ]
      logtarget <- logtarget[ancestors]
      logproposal <- logproposal[ancestors]
      weights_current <- rep(1, nparticles)/nparticles 
    }
    
    ### MCMC rejuvenation move
    transition_kernel <- hmc(target, proposal, mcmc_kernel$parameters, gamma_current)$kernel
    for(i in 1:mcmc_kernel$nmcmc){
      rejuv_result <- transition_kernel(xparticles, logtarget, logproposal)
      xparticles <- rejuv_result$chain_state
      logtarget <- rejuv_result$logtarget
      logproposal <- rejuv_result$logproposal
    }
  }

  return(list(log_ratio_normconst = log_ratio_normconst, temperature = temperature, cess = cess))
}

