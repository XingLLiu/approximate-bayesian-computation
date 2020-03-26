library(rje)
source("src/metrics.R")

# rej_abc <- function(N, epsilon, y, rprior, simulate, kernel="uniform", sumstat="original"){

#     # initialize sample vector
#     samples <- matrix(0, nrow = N, ncol = length(rprior(1)))
#     # initialize acceptance rate vector
#     accept_rate <- rep(0, N)
#     # compute summary statistic of y
#     y_summary <- sumstat_fun(y, sumstat)

#     print(sprintf("Begin soft-ABC. Kernel: %s; epsilon: %f", kernel, epsilon))
#     for (i in 1:N) {
#         # initialize criterion
#         rho <- epsilon + 1
#         accept_prob <- rho <= epsilon

#         # initialize trial index
#         trials <- 0
#         while (runif(n = 1) > accept_prob){
#             # sample candidate parameter
#             theta <- rprior(1)
#             # sample synthetic data
#             z <- simulate(n = length(y), theta = theta)
#             # compute acceptance probability
#             z_summary <- sumstat_fun(z, sumstat)
#             accept_prob <- discrep_kernel(z_summary, y_summary, epsilon, kernel)
#             # update trial index
#             trials <- trials + 1
#         }

#         samples[i, ] <- theta
#         accept_rate[i] <- 1/trials

#         printPercentage(i, N)
#     }
    
#     return(list("samples" = samples, "trials" = accept_rate))
# }


rej_abc <- function(args)
#'Rejection ABC
#' 
#' @description Rejection ABC
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
  # initialize acceptance rate vector
  accept_rate <- rep(0, N)

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

  print(sprintf("Begin rej-ABC. Kernel: %s; epsilon: %f", kernel, epsilon))
  if (kernel != "uniform"){
    warning("Rejection ABC is inefficient with a general kernel! Change to a unform kernel instead.")
  }
  for (i in 1:N) {
    # initialize criterion
    accept_prob <- 0

    # initialize trial index
    trials <- 0
    while (runif(1) > accept_prob){
        # sample candidate parameter
        theta <- rprior(1)
        # sample synthetic data
        z <- simulate(n = length(y), theta = theta)
        # compute acceptance probability
        z_summary <- sumstat_fun(z)
        # compute discrepancy measure
        rho <- discrepancy(y_summary, z_summary)
        accept_prob <- discrep_kernel(rho, epsilon, kernel)
        # update trial index
        trials <- trials + 1
    }

    # compute acceptance rate
    samples[i, ] <- theta
    accept_rate[i] <- 1/trials

    printPercentage(i, N)
  }

  return(list("samples" = samples, "trials" = accept_rate))
}