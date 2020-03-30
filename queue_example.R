# M/G/1 queuing model
# Only theta_\ast = (\theta_1, \theta2 - \theta1, \theta3)
library(winference)
library(ggplot2)
source("src/rej_abc.R")
source("src/soft_abc.R")
source("src/sabc.R")
source("src/mgqueue/mgqueue.R")
set.seed(2020)
theta_star <- list(theta = c(1, 4, 0.2))
theta_names <- c("theta1", "theta2minus1", "theta2")
theta_star_list <- list(theta1 = theta_star$theta[1],
                          theta2minus1 = theta_star$theta[2],
                          theta3 = theta_star$theta[3]
                        )
nobservation <- 50
nthetas <- 1024
maxsimulation <- 10^6
resultsprefix <- "results/queue/"
plotprefix <- "plots/queue/"


# abc posterior
queue_model <- init_queue_model_smc(nobservation, theta_star$theta)

# rprior 
rprior <- queue_model$rprior

# log dprior
dprior <- queue_model$dprior

# generating process
simulate <- queue_model$simulate

# compute or load observed y
# y <- queue_model$y
# save(y, file = "results/queue/y.RData")
load(file = "results/queue/y.RData")

# summary statistic
sumstat <- function(z){
  return(matrix(quantile(z, probs = seq(0, 1, length.out = 20)), ncol = 1))
}

# euclidean discrepancy
y_summary <- sumstat(y)
eucdiscrep <- function(z){
  z_summary <- sumstat(z)
  return(l2norm(z_summary, y_summary))
}


# laod PMMH for the true posteriors 
# uncomment the following two lines upon first run
# source("src/mgqueue/queue_wsmc_marginal_intermediate.R")
# source("src/mgqueue/queue_pmmh_intermediate.R")
filename = paste0(resultsprefix, "queue_intermediate_pmmh.RData")
load(file = filename)
results_pmmh = res
rm(res)
mcmc.df = mhchainlist_to_dataframe(results_pmmh$chains)
mcmc.df = mcmc.df %>% filter(iteration %% 10 == 1)
names(mcmc.df) = c("ichain", "iteration", "theta1", "theta2minus1", "theta3")
mcmc.df$theta2 = mcmc.df$theta2minus1 + mcmc.df$theta1


# initialize dataframes to store samples
method_names <- c(paste("Euc. Summary"), "MMD", "Wasserstein", "KL Divergence")
abc_df <- data.frame(methods = rep(method_names, each = nthetas),
                     samples.theta1 = NA,
                     samples.theta2minus1 = NA,
                     samples.theta3 = NA
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

mmd_out <- sabc(args_mmd, maxsimulation= maxsimulation)
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
abc_df$samples.theta2 <- abc_df$samples.theta2minus1 + abc_df$samples.theta1
write.csv(abc_df, paste0(resultsprefix, "abc_df.csv"), row.names = FALSE)
abc_df <- read.csv(paste0(resultsprefix, "abc_df.csv"))


# plot results
my_colours <- init_colours()
pdf(paste0(plotprefix, "posterior_densities.pdf"), width = 14)
g1 <- ggplot(data = mcmc.df, aes(x = theta1, fill = "Posterior", colour = "Posterior"), alpha = 0.5) +
        geom_density() +
        geom_density(data = abc_df, aes(x = samples.theta1, fill = methods, colour = methods), alpha = 0.5) +
        scale_color_manual(name = "", values = my_colours) +
        scale_fill_manual(name = "", values = my_colours) +
        labs(x = "theta1") +
        change_sizes(16, 20) +
        geom_vline(xintercept = theta_star$theta[1], linetype = 2) +
        add_legend(0.65, 0.95)
g2 <- ggplot(data = mcmc.df, aes(x = theta2minus1, fill = "Posterior", colour = "Posterior"), alpha = 0.5) +
        geom_density() +
        geom_density(data = abc_df, aes(x = samples.theta2minus1, fill = methods, colour = methods), alpha = 0.5) +
        scale_color_manual(name = "", values = my_colours) +
        scale_fill_manual(name = "", values = my_colours) +
        labs(x = "theta2 - theta1") +
        change_sizes(16, 20) +
        geom_vline(xintercept = theta_star$theta[2], linetype = 2) +
        theme(legend.position = "none") 
g3 <- ggplot(data = mcmc.df, aes(x = theta3, fill = "Posterior", colour = "Posterior"), alpha = 0.5) +
        geom_density() +
        geom_density(data = abc_df, aes(x = samples.theta3, fill = methods, colour = methods), alpha = 0.5) +
        scale_color_manual(name = "", values = my_colours) +
        scale_fill_manual(name = "", values = my_colours) +
        labs(x = "theta3") +
        change_sizes(16, 20) +
        geom_vline(xintercept = theta_star$theta[3], linetype = 2) +
        theme(legend.position = "none")
gridExtra::grid.arrange(g1, g2, g3, ncol = 3)
dev.off()



