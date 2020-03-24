#  #'Metropolis-Hastings algorithm
  #' 
  #'@description Implement the MH algorithm for drawing samples from the exact ABC posterior density.

source("src/mcmc/MH_process.R")

MH <- function(nsamples, start, proposal_draw, proposal_density, target_density, burnin=0, thinning=1){

  ntotal <- (nsamples - 1) * thinning + burnin + 1
  samples <-  rep(NA, ntotal)
  accept_rate <- rep(0, ntotal)

  samples[1] <- start
  for (i in 2:ntotal){
    theta_new <- proposal_draw(samples[i - 1])
    accept_prob <- proposal_density(samples[i - 1], theta_new) * target_density(theta_new) / (
                      proposal_density(theta_new, samples[i - 1]) * target_density(samples[i - 1])       
                    )

    if (runif(1) < accept_prob){
      samples[i] <- theta_new
    } else{
      samples[i] <- samples[i - 1]
    }

    if (i %% 5e4 == 0){
      print(sprintf("[%i / %i]", i, ntotal))
    }
  }

  if (burnin == 0 & thinning == 1){
    return(samples)
  } else{
    return( MH_process(samples, nsamples, burnin, thinning) )
  }
}
