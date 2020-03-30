# Functions for the M\G\1 queueing model

init_queue_model_smc <- function(n, true_theta)
#' @description Initialize rprior, simulate and y for the M/G/1 queue model.
{
  simulate <- function(theta){
    theta1 <- theta[1]
    theta2 <- theta[2] + theta[1]
    theta3 <- theta[3]
    u <- runif(n, theta1, theta2)
    v <- y <- x <- rep(0, n)
    v[1] <- rexp(n = 1, rate = theta3)
    y[1] <- u[1] + v[1]
    x[1] <- y[1]
    for (t in 2:n){
      v[t] <- v[t-1] + rexp(n = 1, rate = theta3)
      y[t] <- u[t] + max(0, v[t] - x[t-1])
      x[t] <- x[t-1] + y[t]
    }
    return(matrix(y, ncol = 1))
  }

  y <- simulate(true_theta)

  # prior sampler
  upperbd <- min(y, 10)
  # upperbd <- 10
  rprior <- function(N, hyperparams){
    theta1 <- runif(n = N, min = 0, max = upperbd)
    theta2minus1 <- runif(n = N, min = 0, max = 10)
    theta3 <- runif(n = N, min = 0, max = 1/3)
    return(cbind(theta1, theta2minus1, theta3))
  }
  
  # prior density
  dprior <- function(theta, hyperparams){
    evals <- dunif(theta[, 1], min = 0, max = upperbd, log = TRUE) + 
              dunif(theta[, 2], min = 0, max = 10, log = TRUE) +
              dunif(theta[, 3], min = 0, max = 1/3, log = TRUE)
    return(evals)
  }


  return(list(simulate = simulate, rprior = rprior, dprior = dprior, y = y))
}


rearrange_queue_params <- function(samples)
#' @description Rearrange the parameters for the M/G/1 queue model from
#' (\theta_1, \theta_2 - \theta_1, \theta3) to (\theta_1, \theta_2, \theta3).
{
  samples[, 2] <- samples[, 2] + samples[, 1]
  return(samples)
}


rearrange_queue_params <- function(samples)
#' @description Rearrange the parameters for the M/G/1 queue model from
#' (\theta_1, \theta_2 - \theta_1, \theta3) to (\theta_1, \theta_2, \theta3).
{
  samples[, 2] <- samples[, 2] + samples[, 1]
  return(samples)
}


init_queue_model <- function(nobservation, true_theta)
#' @description Initialize rprior, simulate and y for the M/G/1 queue model.
{
  simulate <- function(n, theta){
    theta1 <- theta[1]
    theta2 <- theta[2] + theta[1]
    theta3 <- theta[3]
    u <- runif(n, theta1, theta2)
    v <- y <- x <- rep(0, n)
    v[1] <- rexp(n = 1, rate = theta3)
    y[1] <- u[1] + v[1]
    x[1] <- y[1]
    for (t in 2:n){
      v[t] <- v[t-1] + rexp(n = 1, rate = theta3)
      y[t] <- u[t] + max(0, v[t] - x[t-1])
      x[t] <- x[t-1] + y[t]
    }
    return(matrix(y, ncol = 1))
  }

  y <- simulate(nobservation, true_theta)
  # Gamma prior
  upperbd <- min(y, 10)
  # upperbd <- 10
  rprior <- function(n){
    theta1 <- runif(n = n, min = 0, max = upperbd)
    theta2minus1 <- runif(n = n, min = 0, max = 10)
    theta3 <- runif(n = n, min = 0, max = 1/3)
    return(cbind(theta1, theta2minus1, theta3))
  }
  return(list(simulate = simulate, rprior = rprior, y = y))
}

