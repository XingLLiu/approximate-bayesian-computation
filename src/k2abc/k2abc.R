source("src/mmd/mmd.R")
source("src/k2abc/k2_posterior_sampler.R")

k2abc <- function(args) 
# args must contains
# alpha, beta: parameters for gamma prior
# y: observed sample
# epsilon: K2-ABC parameter
# kernel: a S4 object of class 'kernel' as in 'kernlab'
# nsample_out: number of posterior samples
{

    y <- args$y
    N <- args$nsample_out
    # initialize sample vector
    samples <- rep(NA, N)
    # initialize weights
    weights <- rep(NA, N)

    # sample from prior
    samples <- rgamma(n = N, shape = args$alpha, rate = args$beta)
    for (i in 1:N){
      # sample synthetic data
      z <- rexp(length(y), rate = samples[i])
      # compute weight
      weights[i] <- exp(- mmd_sq(args$kernel, y, z) / args$epsilon)

      if (i %% 100 == 0){
        cat(sprintf("[%i / %i]", i, N), "\n")
      }
    }

    # return normalized weights and prior samples
    weights <- weights / sum(weights)
    return(list("weights" = weights, 
                "prior_samples" = samples
                )
          )
}



