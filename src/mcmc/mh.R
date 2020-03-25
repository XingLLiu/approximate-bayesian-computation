#  #'Metropolis-Hastings algorithm
  #' 
  #'@description Implement the MH algorithm for drawing samples from the exact ABC posterior density.

library(rje)
library(winference)
source("src/mcmc/mh_process.R")

# mh <- function(nsamples, start, proposal_draw, proposal_density, target_density, burnin=0, thinning=1){
#   ntotal <- (nsamples - 1) * thinning + burnin + 1
#   samples <-  rep(NA, ntotal)
#   accept_num <- 1

#   samples[1] <- start
#   for (i in 2:ntotal){
#     theta_new <- proposal_draw(samples[i - 1])
#     accept_prob <- proposal_density(samples[i - 1], theta_new) * target_density(theta_new) / (
#                       proposal_density(theta_new, samples[i - 1]) * target_density(samples[i - 1])       
#                     )

#     if (runif(1) < accept_prob){
#       samples[i] <- theta_new
#       accept_num <- accept_num + 1
#     } else{
#       samples[i] <- samples[i - 1]
#     }

#     if (i %% 5e4 == 0){
#       print(sprintf("[%i / %i]", i, ntotal))
#     }
#   }

#   print(sprintf("acceptance rate: %.3f", accept_num / ntotal))
#   if (burnin == 0 & thinning == 1){
#     return(samples)
#   } else{
#     return( mh_process(samples, nsamples, burnin, thinning) )
#   }
# }



mh <- function (observations, target, tuning_parameters, savefile = NULL, verbose = TRUE) 
{
    posterior <- function(thetas) {
        logdens <- target$dprior(thetas, target$parameters)
        which.ok <- which(is.finite(logdens))
        if (length(which.ok) > 0) {
            theta.ok <- thetas[which.ok, , drop = FALSE]
            logdens[which.ok] <- logdens[which.ok] + 
                                  target$loglikelihood(theta.ok, observations, target$parameters)
        }
        return(logdens)
    }
    niterations <- tuning_parameters$niterations
    nchains <- tuning_parameters$nchains
    cov_proposal <- tuning_parameters$cov_proposal
    p <- ncol(tuning_parameters$init_chains)
    chains <- rep(list(matrix(nrow = niterations, ncol = p)), nchains)
    current_chains <- matrix(nrow = nchains, ncol = p)
    current_chains <- matrix(tuning_parameters$init_chains, nrow = nchains, ncol = p)
    for (ichain in 1:nchains) {
        chains[[ichain]][1, ] <- current_chains[ichain, ]
    }
    current_dtarget <- posterior(current_chains)
    naccepts <- 0
    for (iteration in 2:niterations) {
        # if ((iteration %% max(1, floor(niterations/100)) == 1) && 
        #     (verbose)) {
        #     cat("iteration ", iteration, "/", niterations, "\n")
        #     cat("average acceptance:", naccepts/(iteration * 
        #         nchains) * 100, "%\n")
        # }
        if (verbose) {
          printPercentage(iteration, niterations)
        }
        if (iteration > 250 && tuning_parameters$adaptation > 
            0 && (iteration%%tuning_parameters$adaptation) == 
            0) {
            mcmc_samples <- foreach(ichain = 1:nchains, .combine = rbind) %do% 
                {
                  matrix(chains[[ichain]][max(1, iteration - 
                    tuning_parameters$adaptation):(iteration - 
                    1), ], ncol = p)
                }
            cov_proposal <- cov(mcmc_samples)/p
        }
        proposals <- current_chains + fast_rmvnorm(nchains, rep(0, p), cov_proposal)
        proposal_dtarget <- posterior(proposals)
        acceptance_ratios <- (proposal_dtarget - current_dtarget)
        uniforms <- runif(n = nchains)
        accepts <- (log(uniforms) < acceptance_ratios)
        naccepts <- naccepts + sum(accepts)
        current_chains[accepts, ] <- proposals[accepts, ]
        if (is.null(dim(current_chains))) 
          current_chains <- matrix(current_chains, ncol = p)
          current_dtarget[accepts] <- proposal_dtarget[accepts]
        for (ichain in 1:nchains) {
            chains[[ichain]][iteration, ] <- current_chains[ichain, ]
        }
        if (!is.null(savefile) && iteration%%1000 == 1) {
            mh_results <- list(chains = chains,
                               naccepts = naccepts, 
                               cov_proposal = cov_proposal,
                               iteration = iteration
                              )
            save(mh_results, file = savefile)
        }
    }
    cat("average acceptance:", naccepts/(niterations * nchains) * 
        100, "%\n")
    return(list(chains = chains, naccepts = naccepts, cov_proposal = cov_proposal))
}