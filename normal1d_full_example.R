# Normal(theta_\ast, \sigma_0^2) with Normal(m_0, \tau_0^2) prior
# Only theta_\ast is unknown
library(ggplot2)
set.seed(2020)
theta_star <- list(theta = 0, sigma = 2)
hyperparams <- list(m = 1, tau = 2)
epsilon <- c(0.05)   # c(0.05, 0.1, 0.5, 1)
nobservation <- 500
nthetas <- 1024
y <- matrix(rnorm(nobservation, mean = theta_star$theta, sd = theta_star$sigma), ncol = 1)


# Gamma prior
rprior <- function(n){ return(rnorm(n, mean = hyperparams$m, sd = hyperparams$tau)) }
# generating process
simulate <- function(n, theta){
  return(matrix(rnorm(n = n, mean = theta, sd = theta_star$sigma), ncol = 1))
}

# initialize dataframe to store the true density
thetavals <- seq(-2, 2, length.out = nthetas)
posterior_df <- data.frame(thetavals = thetavals,
                           true_posterior = NA,
                           abc_posterior = NA
                          )
# true posterior
post_sd <- 1 / sqrt( hyperparams$tau^(-2) + (theta_star$sigma^2 / nobservation)^(-1) )
post_mean <- (
              hyperparams$m * hyperparams$tau^(-2) + 
              mean(y) * (theta_star$sigma^2 / nobservation)^(-1)
             ) * post_sd^2
posterior_df$true_posterior <- dnorm(thetavals, mean = post_mean, sd = post_sd)

# exact abc posterior
abc_post_sd <- 1 / sqrt( hyperparams$tau^(-2) + (theta_star$sigma^2 / nobservation + epsilon[1]^2)^(-1) )
abc_post_mean <- (
              hyperparams$m * hyperparams$tau^(-2) + 
              mean(y) * (theta_star$sigma^2 / nobservation + epsilon[1]^2)^(-1)
             ) * post_sd^2

abcposterior_func <- function(theta){
  return(dnorm(theta, mean = abc_post_mean, sd = abc_post_sd))
}

posterior_df$abc_posterior <- abcposterior_func(thetavals)


# initialize arguments
source("src/rej_abc.R")
args_rej <- list(nthetas = nthetas, y = y,
                  rpiror = rprior,
                  simulate = simulate,
                  kernel = "uniform",
                  discrepancy = l2norm,
                  sumstat = mean,
                  epsilon = epsilon[1]
                )

# initialize dataframes to store samples
method_names <- c(paste("epsilon =", epsilon, " uniform"), "MMD", "Wasserstein", "KL divergence")
abc_df <- data.frame(method = rep(method_names, each = nthetas),
                     samples = NA
                    )

# function to set up the index
index <- function(i){
  return((1 + (i - 1) * nthetas): (i * nthetas))
}

# rejection abc
samples_df <- rej_abc(args_rej)
abc_df$samples[index(1)] <- samples_df$samples


# K2 ABC
source("src/mmd/mmdsq_c.R")
bandwidth <- median(apply(y, 1, l1norm))
mmdsq <- function(y, z){ 
  return(mmdsq_c(y, z, bandwidth)) 
}

args_mmd <- list(nthetas = nthetas, y = y,
                  rpiror = rprior,
                  simulate = simulate,
                  kernel = "uniform",
                  discrepancy = mmdsq,
                  epsilon = 0.01
                )

k2abc_out <- rej_abc(args_mmd)
abc_df$samples[index(2)] <- k2abc_out$samples


# KL ABC
# function to compute 1-Wasserstein distance between observed data and fake data given as argument
wdistance <- function(y_sorted, y_fake){
  y_fake <- sort(y_fake)
  return(mean(abs(y_sorted - y_fake)))
} 

args_wabc <- list(nthetas = nthetas, y = sort(y),
                  rpiror = rprior,
                  simulate = simulate,
                  kernel = "uniform",
                  discrepancy = wdistance,
                  epsilon = 0.1
                 ) 

wabc_out <- rej_abc(args_wabc)
abc_df$samples[index(3)] <- wabc_out$samples


# KL ABC
# function to compute 1-Wasserstein distance between observed data and fake data given as argument
kldist <- function(y, z){
  return(FNN::KLx.divergence(y, z, k = 1))
} 

args_kl <- list(nthetas = nthetas, y = sort(y),
                rpiror = rprior,
                simulate = simulate,
                kernel = "uniform",
                discrepancy = kldist,
                epsilon = 0.01
               ) 

klabc_out <- rej_abc(args_kl)
abc_df$samples[index(4)] <- klabc_out$samples


# plot results
plt <- ggplot(abc_df) +
        geom_density(aes(x = samples, colour = method)) +
        geom_line(data = posterior_df, aes(x = thetavals, y = true_posterior)) +
        geom_vline(xintercept = theta_star$theta, linetype = "dashed") +
        labs(x = "theta", y = "density") +
        theme(
              legend.position = c(.95, .95),
              legend.justification = c("right", "top"),
              legend.title=element_blank()
             ) +
        guides(color = guide_legend(override.aes = list(linetype = "solid")))

ggsave(plt, file = "plots/soft_abc/normal1d_full_eg.pdf", height = 5)

