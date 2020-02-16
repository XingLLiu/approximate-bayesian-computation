# Exponential with Gamma prior
set.seed(1)
lambda0 <- 2
alpha0 <- 1
beta0 <- 1
n <- 10
y <- rexp(n, lambda0)

# Rejection with identity kernel
source("src/metrics.R")
source("src/abc_rej.R")
samples_id <- abc_sample(N = 1000, epsilon = 0.15, y = y, alpha = alpha0, beta = beta0)
# Rejection with Gaussian kernel
samples_rbf <- abc_sample(N = 1000, epsilon = 0.05, y = y, alpha = alpha0, beta = beta0, kernel = "rbf")
# Rejection with Epanechnikov
samples_epan <- abc_sample(N = 1000, epsilon = 0.15, y = y, alpha = alpha0, beta = beta0, kernel = "epan")
# Rejection with identity kernel and sufficient summary statistic
samples_id_sufficient <-  abc_sample(N = 1000, epsilon = 0.01, y = y, alpha = alpha0, beta = beta0, summary_stat = "mean")
# ABC MCMC with Gaussian proposal
source("src/abc_mcmc.R")
samples_mcmc <- abc_mcmc_sampler(N = 1000, epsilon = 0.1, y = y, summary_stat = "mean", burnin = 30000, thinning = 100,
                                alpha = alpha0, beta = beta0)

# Plot results
samples_density_id <- density(samples_id$samples)
samples_density_rbf <- density(samples_rbf$samples)
samples_density_epan <- density(samples_epan$samples)
samples_density_id_sufficient <- density(samples_id_sufficient$samples)
samples_density_mcmc <- density(samples_mcmc)

plot(samples_density_id$x, dgamma(samples_density_id$x, shape = alpha0 + n, rate = beta0 + sum(y)), type = "l")
lines(samples_density_id$x, samples_density_id$y, col = "blue")
lines(samples_density_rbf$x, samples_density_rbf$y, col = "red")
lines(samples_density_epan$x, samples_density_epan$y, col = "green")
lines(samples_density_id_sufficient$x, samples_density_id_sufficient$y, col = "orange")
lines(samples_density_mcmc$x, samples_density_mcmc$y, col = "cyan")

# convergence of MCMC
plot(cumsum(samples_mcmc)/1:length(samples_mcmc))
acf(samples_mcmcs)

libs <- c("Rcpp", "RcppEigen", "doParallel", "doRNG", "foreach", "ggplot2", "ggthemes", "dplyr", "reshape2", "BH", "transport")
