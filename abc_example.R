# Exponential with Gamma prior
library(kernlab)
library(winference)

set.seed(1)
lambda0 <- 2
alpha0 <- 1
beta0 <- 1
nobservation <- 500
y <- rexp(nobservation, lambda0)

# Rejection with identity kernel
source("src/abc_rej.R")
samples_id <- abc_sample(N = 1000, epsilon = 0.15, y = y, alpha = alpha0, beta = beta0)
# Rejection with Gaussian kernel
samples_rbf <- abc_sample(N = 1000, epsilon = 0.05, y = y, alpha = alpha0, beta = beta0, kernel = "gaussian")
# Rejection with Epanechnikov
samples_epan <- abc_sample(N = 1000, epsilon = 0.15, y = y, alpha = alpha0, beta = beta0, kernel = "epan")
# Rejection with identity kernel and sufficient summary statistic
samples_id_sufficient <-  abc_sample(N = 1000, epsilon = 0.01, y = y, alpha = alpha0, beta = beta0, sumstat = "mean")
# ABC MCMC with Gaussian proposal
source("src/abc_mcmc.R")
samples_mcmc <- abc_mcmc_sampler(N = 1000, epsilon = 0.1, y = y, sumstat = "mean", burnin = 30000, thinning = 100,
                                alpha = alpha0, beta = beta0)

# Plot results
samples_density_id <- density(samples_id$samples)
samples_density_rbf <- density(samples_rbf$samples)
samples_density_epan <- density(samples_epan$samples)
samples_density_id_sufficient <- density(samples_id_sufficient$samples)
samples_density_mcmc <- density(samples_mcmc)

plot(samples_density_id$x, dgamma(samples_density_id$x, shape = alpha0 + nobservation, rate = beta0 + sum(y)), type = "l")
lines(samples_density_id$x, samples_density_id$y, col = "blue")
lines(samples_density_rbf$x, samples_density_rbf$y, col = "red")
lines(samples_density_epan$x, samples_density_epan$y, col = "green")
lines(samples_density_id_sufficient$x, samples_density_id_sufficient$y, col = "orange")
lines(samples_density_mcmc$x, samples_density_mcmc$y, col = "cyan")

# convergence of MCMC
plot(cumsum(samples_mcmc)/1:length(samples_mcmc))
acf(samples_mcmcs)


# WABC
# function to generate from prior distribution
rprior <- function(N, parameters){
  particles <- matrix(nrow = N, ncol = 1)
  particles[,1] <- rgamma(N, shape = parameters$alpha0, rate = parameters$beta0)
  return(particles)
}

# function to evaluate prior log-density
dprior <- function(thetas, parameters){
  logdensities <- dgamma(thetas[,1], shape = parameters$alpha0, rate = parameters$beta0, log = TRUE)
  return(logdensities)
}

# generative model, given a parameter theta
simulate <- function(theta){
  observations <- rexp(nobservation, theta)
  return(observations)
}

# collect everything in a list
target <- list(simulate = simulate, rprior = rprior, dprior = dprior,
               parameter_names = c("lambda"),
               parameters = list(alpha0 = alpha0, beta0 = beta0),
               thetadim = 1, ydim = 1)

# function for calculating the distance
y_sorted <- sort(y)
# function to compute 1-Wasserstein distance between observed data and fake data given as argument
compute_distance <- function(y_fake){
  y_fake <- sort(y_fake)
  return(mean(abs(y_sorted - y_fake)))
} 

# specify hyper-parameters for SMC sampler
param_algo <- list(nthetas = 1024, nmoves = 1, proposal = mixture_rmixmod(),
                   minimum_diversity = 0.5, R = 2, maxtrials = 1e5)

# run wsmc
wsmcresults <- wsmc(compute_distance, target, param_algo, maxtime = 30)
# plot result
plot_marginal(wsmcresults, i = 1)

# continue running wsmc
wsmcresults_cont <- wsmc_continue(wsmcresults, maxtime = 30)
plot_marginal(wsmcresults_cont, i = 1)


# K2-ABC
source("src/k2abc/k2abc.R")
rbf <- rbfdot(0.5)
k2_epsilon <- 1e-4 # (4/(3 * nobservation))^(1/5)
k2_samples <- k2abc(list("nsample_out" = 500, "y" = y, "alpha" = alpha0, "beta" = beta0, "epsilon" = k2_epsilon, "kernel" = rbf))

# samples_density_k2abc <- density(sample(size = 1000, x = k2_samples$prior_samples, prob = k2_samples$weights, replace = TRUE))
samples_density_k2abc <- density(draw_k2_sample(1000, k2_samples))

plot(samples_density_id$x, dgamma(samples_density_id$x, shape = alpha0 + nobservation, rate = beta0 + sum(y)), type = "l")
lines(samples_density_k2abc, col = "red")


