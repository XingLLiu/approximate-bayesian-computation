# Normal(theta_\ast, \sigma_0^2) with Normal(m_0, \tau_0^2) prior
# Only theta_\ast is unknown
library(ggplot2)
library(winference)
source("src/rej_abc.R")
source("src/sabc.R")
set.seed(2020)
theta_star <- list(theta = 0, sigma = 2)
theta_names <- c("theta")
theta_star_list <- list(theta = theta_star$theta[1])
hyperparams <- list(m = 1, tau = 2)
nobservation <- 500
nthetas <- 1024
resultsprefix <- "results/soft_abc/"
plotsprefix <- "plots/soft_abc/"


# rprior
rprior <- function(n, params){ 
  return(matrix(rnorm(n, mean = hyperparams$m, sd = hyperparams$tau), ncol = 1)) 
}
# log dprior
dprior <- function(theta, params){
  eval <- dnorm(dnorm(theta[, 1], mean = hyperparams$m, sd = hyperparams$tau))
}
# generating process
simulate <- function(theta){
  return(matrix(rnorm(n = nobservation, mean = theta, sd = theta_star$sigma), ncol = 1))
}

y <- simulate(theta_star$theta)

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
abcposterior_func <- function(theta, epsilonval){
  abc_post_sd <- 1 / sqrt( hyperparams$tau^(-2) + (theta_star$sigma^2 / nobservation + epsilonval^2)^(-1) )
  abc_post_mean <- (
                hyperparams$m * hyperparams$tau^(-2) + 
                mean(y) * (theta_star$sigma^2 / nobservation + epsilonval^2)^(-1)
              ) * post_sd^2
  return(dnorm(theta, mean = abc_post_mean, sd = abc_post_sd))
}

posterior_df$abc_posterior <- abcposterior_func(thetavals, 0.1)


# euclidean discrepancy
y_summary <- mean(y)
eucdiscrep <- function(z){
  z_summary <- mean(z)
  return(abs(z_summary - y_summary))
}


# initialize dataframes to store samples
method_names <- c("Euc. Summary", "MMD", "Wasserstein", "KL Divergence")
abc_df <- data.frame(methods = rep(method_names, each = nthetas),
                     samples.theta = NA
                    )

# rej ABC
args_rej <- list(nthetas = nthetas,
                  rprior = rprior,
                  dprior = dprior,
                  simulate = simulate,
                  discrepancy = eucdiscrep,
                  parameter_names = theta_names,
                  thetadim = length(theta_names),
                  ydim = ncol(y)
                )

rej_out <- sabc(args_rej, maxsimulation = maxsimulation)
abc_df[index(1, nthetas), 2:ncol(abc_df)] <- sabc_get_last_samples(rej_out)[, theta_names]


# K2 ABC
source("src/mmd/mmdsq_c.R")
bandwidth <- median(apply(y, 1, l1norm))
mmdsq <- function(z){ 
  return(mmdsq_c(y, z, bandwidth)) 
}

args_mmd <- list(nthetas = nthetas,
                  rprior = rprior,
                  dprior = dprior,
                  simulate = simulate,
                  discrepancy = mmdsq,
                  parameter_names = theta_names,
                  thetadim = length(theta_names),
                  ydim = ncol(y)
                )

mmd_out <- sabc(args_mmd, maxsimulation = maxsimulation)
abc_df[index(2, nthetas), 2:ncol(abc_df)] <- sabc_get_last_samples(mmd_out)[, theta_names]


# WABC
# function to compute 1-Wasserstein distance between observed data and fake data given as argument
y_sorted <- sort(y)
wdistance <- function(y_fake){
  y_fake <- sort(y_fake)
  return(mean(abs(y_sorted - y_fake)))
} 

args_wabc <- list(nthetas = nthetas,
                  rprior = rprior,
                  dprior = dprior,
                  simulate = simulate,
                  discrepancy = wdistance,
                  parameter_names = theta_names,
                  thetadim = length(theta_names),
                  ydim = ncol(y)
                )

wabc_out <- sabc(args_wabc, maxsimulation= maxsimulation)
abc_df[index(3, nthetas), 2:ncol(abc_df)] <- sabc_get_last_samples(wabc_out)[, theta_names]


# KL ABC
kldist <- function(z){
  return(FNN::KLx.divergence(y, z, k = 1))
} 

args_kl <- list(nthetas = nthetas,
                  rprior = rprior,
                  dprior = dprior,
                  simulate = simulate,
                  discrepancy = kldist,
                  parameter_names = theta_names,
                  thetadim = length(theta_names),
                  ydim = ncol(y)
                )

klabc_out <- sabc(args_kl, maxsimulation= maxsimulation)
abc_df[index(4, nthetas), 2:ncol(abc_df)] <- sabc_get_last_samples(klabc_out)[, theta_names]

# save results
write.csv(abc_df, paste0(resultsprefix, "abc_df.csv"), row.names = FALSE)
abc_df <- read.csv(paste0(resultsprefix, "/abc_df.csv"))


# plot results
my_colours <- init_colours()
plt <- ggplot(abc_df) +
        geom_density(aes(x = samples.theta, fill = methods, colour = methods), alpha = 0.5) +
        geom_line(data = posterior_df, aes(x = thetavals, y = true_posterior)) +
        geom_ribbon(data = posterior_df, aes(x = thetavals, ymax = true_posterior), ymin = 0, alpha = 0.5) +
        geom_vline(xintercept = theta_star$theta, linetype = "dashed") +
        labs(x = "theta", y = "density") +
        xlim(-1, 1) +
        scale_color_manual(name = "", values = my_colours) +
        scale_fill_manual(name = "", values = my_colours) +
        change_sizes(16, 20) +
        geom_vline(xintercept = theta_star$theta[1], linetype = 2) +
        add_legend(0.95, 0.95)

ggsave(plt, file = paste0(plotsprefix, "normal1d_full_eg.pdf"), height = 5)

