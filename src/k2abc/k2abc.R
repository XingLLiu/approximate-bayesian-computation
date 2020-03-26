library(Matrix)
library(rje)
source("src/mmd/mmdsq_c.R")
source("src/metrics.R")
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
    N <- args$nthetas
    # initialize sample vector
    samples <- rep(NA, N)
    # initialize weights
    weights <- rep(NA, N)
    if (is.null(args$kernel_epsilon)) {
      kernel_epsilon <- median(apply(y, 1, l1norm))
    } else {
      kernel_epsilon <- args$kernel_epsilon
    }

    # sample from prior
    samples <- args$rprior(n = N) # rgamma(n = N, shape = args$alpha, rate = args$beta)
    for (i in 1:N){
      # sample synthetic data
      z <- args$simulate(nrow(y), samples[i]) # rexp(length(y), rate = samples[i])
      # compute weight
      weights[i] <- exp(- mmdsq_c(y, z, kernel_epsilon) / (2 * args$epsilon^2))

      printPercentage(i, N)
    }

    # return normalized weights and prior samples
    weights <- weights / sum(weights)
    return(list("samples" = samples, "weights" = weights))
}

