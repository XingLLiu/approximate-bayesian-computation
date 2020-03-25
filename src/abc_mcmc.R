abc_mcmc_sampler <- function(N, epsilon, y, kernel="identity", sumstat="original", burnin=500, thinning=3, alpha, beta){
#'
#' thinning = [int] save 1 sample per "thinning" ones
#' 

    N_total <- (N - 1) * thinning + burnin + 1
    samples <- rep(0, N_total)
    accept_rate <- rep(0, N_total)
    y_summary <- sumstat(y, family = sumstat)

    # initialize chain
    samples[1] <- rgamma(n = 1, shape = alpha, rate = beta)

    # initialize trial index
    i <- 2

    while (i <= N_total) {

        # sample candidate parameter
        theta <- rnorm(n = 1, mean = samples[i - 1])

        # reject proposal if theta is outside of param. space
        if (theta < 0){
            samples[i] <- samples[i - 1]
            next
        }

        # sample synthetic data
        z <- rexp(length(y), rate = theta)

        # compute acceptance probability
        z_summary <- sumstat(z, family = sumstat)
        accept_prob_kernel <- discrep_kernel(z_summary, y_summary, epsilon = epsilon, kernel = kernel)
        
        mcmc_ratio <- dgamma(theta, shape = alpha, rate = beta) * dnorm(samples[i - 1], mean = theta) / ( 
                        dgamma(samples[i - 1], shape = alpha, rate = beta) * dnorm(theta, mean = samples[i - 1]) )

        accept_prob <- min(1, mcmc_ratio * accept_prob_kernel)

        # accept_prob <- min( 1, dgamma(theta, shape = alpha + length(y), rate = beta + sum(y)) / ( 
        #                 dgamma(samples[i - 1], shape = alpha + length(y), rate = beta + sum(y)) ) )


        # store sample
        if (runif(1) < accept_prob){
            samples[i] <- theta
        } else {
            samples[i] <- samples[i - 1]
        }

        # update trial index
        i <- i + 1

        if (i %% 500 == 0){
            cat(sprintf("[%s / %s]\n", i, N_total))
        }
    }
    
    return( samples[-(1:burnin)][seq(from = 1, length.out = N, by = thinning)] )
}