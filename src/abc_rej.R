abc_sample <- function(N, epsilon, y, alpha, beta, kernel="identity", summary_stat="original"){

    # initialize sample vector
    samples <- rep(0, N)
    # initialize acceptance rate vector
    accept_rate <- rep(0, N)
    # compute summary statistic of y
    y_summary <- summary_stat(y, family = summary_stat)

    for (i in 1:N) {
        # initialize criterion
        rho <- epsilon + 1
        accept_prob <- rho <= epsilon

        # initialize trial index
        trials <- 0
        while (runif(n = 1) > accept_prob){
            # sample candidate parameter
            theta <- rgamma(n = 1, shape = alpha, rate = beta)
            # sample synthetic data
            z <- rexp(length(y), rate = theta)
            # compute acceptance probability
            z_summary <- summary_stat(z, family = summary_stat)
            accept_prob <- discrep_kernel(z_summary, y_summary, epsilon, kernel)
            # update trial index
            trials <- trials + 1
        }

        samples[i] <- theta
        accept_rate[i] <- 1/trials

        if (i %% 50 == 0){
            cat(sprintf("[%s / %s]\n", i, N))
        }
    }
    
    return(list("samples" = samples, "trials" = accept_rate))
}
