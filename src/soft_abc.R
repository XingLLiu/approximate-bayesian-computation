library(rje)
source("src/metrics.R")

soft_abc <- function(args)
#'Soft ABC
#' 
#' @description Soft ABC
#' 
#' @param y
#' @param nthetas
#' @param kernel 
#' @param epsilon FUN; lengthscale used for the kernel (decision threshold).
#' @param rprior
#' @param simulate
#' @param discrepancy FUN; discrepancy metric between s(y) and s(z)
#' @param sumstat FUN; default set to be the identity function
#'
#' @return samples and acceptance probability
{
  # initialize number of thetas to sample
  N <- args$nthetas
  # initialize sample vector
  samples <- matrix(0, nrow = N, ncol = length(rprior(1)))
  # initialize weights vector
  weights <- rep(0, N)
  # initialize kernel
  kernel <- args$kernel
  # initialize threshold
  epsilon <- args$epsilon

  # summary statistic function
  if (is.null(args$sumstat)) {
    sumstat_fun <- function(z){ return(z) }
  } else {
    sumstat_fun <- args$sumstat
  }
  # compute summary statistic of y
  y_summary <- sumstat_fun(args$y)

  # prior
  rprior <- args$rprior
  # data generating process
  simulate <- args$simulate
  # discrepancy measure
  discrepancy <- args$discrepancy

  print(sprintf("Begin soft-ABC. Kernel: %s; epsilon: %f", kernel, epsilon))
  if (kernel == "uniform"){
    warning("Soft ABC is inefficient with a uniform kernel! Change to other kernels instead.")
  }
  for (i in 1:N) {
    # sample candidate parameter
    theta <- rprior(1)
    # sample synthetic data
    z <- simulate(n = nrow(y), theta = theta)
    # compute acceptance probability
    z_summary <- sumstat_fun(z)
    # compute discrepancy measure
    rho <- discrepancy(y_summary, z_summary)
    weights[i] <- discrep_kernel(rho, epsilon, kernel)
    samples[i, ] <- theta

    printPercentage(i, N)
  }

  # normalize weights
  weights <- weights / sum(weights)

  return(list("samples" = samples, "weights" = weights))
}