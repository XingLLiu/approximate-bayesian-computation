source("src/blowfly/simulate_c.R")

init_blowfly_model <- function(n)
#' @description Initialize rprior, simulate and y for the blowfly model.
{
  simulate <- function(theta){
    P = theta[1]
    delta = theta[2]
    N0 = theta[3]
    sigma.d = theta[4]
    sigma.p = theta[5]
    tau = round(theta[6])
    burnin <- 50
    if (tau == 0){
      tau <- 1
    }
    lag <- tau

    outputs = matrix(NA, nrow = lag + burnin + n)
    outputs[1:lag] <- rep(180, lag)

    eps_t <- rgamma(burnin + n, shape = 1 / (sigma.d^2), rate = sigma.d^2)
    e_t <- rgamma(burnin + n, shape = 1 / (sigma.p^2), rate = sigma.p^2)

    outputs <- simulate_c(n, burnin, lag, eps_t, e_t, P, delta, N0,
                          sigma.d, sigma.p, tau, outputs)
      # return(simulate_c(n, burnin, lag, theta$P, theta$delta, theta$N0,
    #                   theta$sigma.d, theta$sigma.p, theta$tau)
    #       )
    # for (i in 1:(n + burnin)){
    #   t <- i + lag
    #   tau_t <- t - lag
    #   outputs[t] <- P * outputs[tau_t] * exp(- outputs[tau_t] / N0) * e_t[i] + 
    #                   outputs[t - 1] * exp(- delta * eps_t[i])
    # }
    return(matrix(outputs[-(1: (lag + burnin))], ncol = 1))
  }

  # rprior
  rprior <- function(n, hyperparams){
    sample <- exp(cbind(
                    2 + 2 * rnorm(n),
                    -1 + 0.4 * rnorm(n),
                    5 + 0.5 * rnorm(n),
                    -0.5 + rnorm(n),
                    -0.5 + rnorm(n),
                    2 + rnorm(n)
                  ))
    return(sample)
  }

  # log dprior
  dprior <- function(theta, hyperparams){
    eval <- dnorm(theta[, 1], 2, 2, log = TRUE) - theta[, 1] +
              dnorm(theta[, 2], -1, 0.4, log = TRUE) - theta[, 2] +
              dnorm(theta[, 3], 5, 0.5, log = TRUE) - theta[, 3] +
              dnorm(theta[, 4], -0.5, log = TRUE) - theta[, 4] +
              dnorm(theta[, 5], -0.5, log = TRUE) - theta[, 5] +
              dnorm(theta[, 6], 2, log = TRUE) - theta[, 6]
    for (i in 1:nrow(theta)){
      if (any(theta[i, ] < 0)){
        eval[i] <- -Inf   
      }
    }
    return(eval)
  }

  return(list(simulate = simulate, rprior = rprior, dprior = dprior))
}

