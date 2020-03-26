library(rje)
source("src/metrics.R")

soft_abc <- function(N, epsilon, y, prior, simulate, kernel="gaussian", sumstat="original"){

  # initialize sample vector
  samples <- matrix(0, nrow = N, ncol = length(prior(1)))
  # initialize weights vector
  weights <- rep(0, N)
  # 
  # compute summary statistic of y
  y_summary <- sumstat_fun(y, sumstat)

  print(sprintf("Begin soft-ABC. Kernel: %s; epsilon: %f", kernel, epsilon))
  if (kernel == "uniform"){
    warning("Soft ABC is inefficient with a uniform kernel! Change to other kernels instead.")
  }
  for (i in 1:N) {
    # initialize criterion
    rho <- epsilon + 1
    accept_prob <- rho <= epsilon

    # sample candidate parameter
    theta <- prior(1)
    # sample synthetic data
    z <- simulate(n = length(y), theta = theta)
    # compute acceptance probability
    z_summary <- sumstat_fun(z, sumstat)
    weights[i] <- discrep_kernel(z_summary, y_summary, epsilon, kernel)
    samples[i, ] <- theta

    printPercentage(i, N)
  }

  # normalize weights
  weights <- weights / sum(weights)

  return(list("samples" = samples, "weights" = weights))
}