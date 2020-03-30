# blowfly data
library(winference)
library(ggplot2)
source("src/rej_abc.R")
source("src/soft_abc.R")
source("src/sabc.R")
source("src/blowfly/blowfly.R")
set.seed(2020)
# theta_star <- list(theta = c(29, 0.2, 260, 0.6, 0.3, 7))
theta_names <- c("P", "delta", "N0", "sigma.d", "sigma.p", "tau")
# theta_star_list <- list(P = theta_star$theta[1],
#                         delta = theta_star$theta[2],
#                         N0 = theta_star$theta[3],
#                         sigma.d = theta_star$theta[4],
#                         sigma.p = theta_star$theta[5],
#                         tau = theta_star$theta[6]
#                         )
# hyperparams <- list(m = 1, tau = 2)
nobservation <- 180
nthetas <- 1024
maxsimulation <- 5 * 10^6
resultsprefix <- "results/blowfly/"
plotprefix <- "plots/blowfly/"


# abc posterior
blowfly_model <- init_blowfly_model(nobservation)

# rprior 
rprior <- blowfly_model$rprior

# log dprior
dprior <- blowfly_model$dprior

# generating process
simulate <- blowfly_model$simulate

# load observed y
load(file = paste0(resultsprefix, "blowfly.rda"))
y <- matrix(blowfly[, 1], ncol = 1)

# summary statistic
sumstat <- function(z){
  s <- rep(NA, 8)
  s[1:4] <- quantile(z, probs = (1:4)/4) / 1000
  s[5:8] <- quantile(diff(z), probs = (1:4) / 4) / 1000
  return(matrix(s, ncol = 1))
}
# euclidean discrepancy
y_summary <- sumstat(y)
eucdiscrep <- function(z){
  z_summary <- sumstat(z)
  return(l2norm(z_summary, y_summary))
}


# initialize dataframes to store samples
method_names <- c("Euc. Summary", "MMD", "Wasserstein", "KL Divergence")
abc_df <- data.frame(methods = rep(method_names, each = nthetas),
                     samples.P = NA,
                     samples.delta = NA,
                     samples.N0 = NA,
                     samples.sigma.d = NA,
                     samples.sigma.p = NA,
                     samples.tau = NA
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
# function to compute 1-Wasserstein distance for 2d data
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

wabc_out <- sabc(args_wabc, maxsimulation = maxsimulation)
abc_df[index(3, nthetas), 2:ncol(abc_df)] <- sabc_get_last_samples(wabc_out)[, theta_names]


# # KL ABC
# kldist <- function(z){
#   return(FNN::KLx.divergence(y, z, k = 1))
# } 

# args_kl <- list(nthetas = nthetas,
#                   rprior = rprior,
#                   dprior = dprior,
#                   simulate = simulate,
#                   discrepancy = kldist,
#                   parameter_names = theta_names,
#                   thetadim = length(theta_names),
#                   ydim = ncol(y)
#                 )

# klabc_out <- sabc(args_kl, maxsimulation= maxsimulation)
# abc_df[index(4, nthetas), 2:ncol(abc_df)] <- sabc_get_last_samples(klabc_out)[, theta_names]


# save results
write.csv(abc_df, paste0(resultsprefix, "abc_df.csv"), row.names = FALSE)
abc_df <- read.csv(paste0(resultsprefix, "/abc_df.csv"))


# plot results
my_colours <- init_colours()
pdf(paste0(plotprefix, "posterior_densities.pdf"), width = 14)
g1 <- ggplot(data = abc_df, aes(x = samples.P, colour = methods, fill = methods)) +
        geom_density(alpha = 0.5) +
        labs(x = "P") +
        change_sizes(16, 20) +
        scale_color_manual(name = "", values = my_colours) +
        scale_fill_manual(name = "", values = my_colours) +
        # geom_vline(xintercept = theta_star_list$P, linetype = 2) +
        theme(legend.position = "none") 
g2 <- ggplot(data = abc_df, aes(x = samples.N0, colour = methods, fill = methods)) +
        geom_density(alpha = 0.5) +
        labs(x = "N0") +
        change_sizes(16, 20) +
        scale_color_manual(name = "", values = my_colours) +
        scale_fill_manual(name = "", values = my_colours) +
        # geom_vline(xintercept = theta_star_list$N0, linetype = 2) +
        theme(legend.position = "none") 
g3 <- ggplot(data = abc_df, aes(x = samples.sigma.d, colour = methods, fill = methods)) +
        geom_density(alpha = 0.5) +
        labs(x = "sigma_d") +
        change_sizes(16, 20) +
        scale_color_manual(name = "", values = my_colours) +
        scale_fill_manual(name = "", values = my_colours) +
        # geom_vline(xintercept = theta_star_list$sigma.d, linetype = 2) +
        theme(legend.position = "none") 
g4 <- ggplot(data = abc_df, aes(x = samples.sigma.p, colour = methods, fill = methods)) +
        geom_density(alpha = 0.5) +
        labs(x = "sigma_p") +
        change_sizes(16, 20) +
        scale_color_manual(name = "", values = my_colours) +
        scale_fill_manual(name = "", values = my_colours) +
        # geom_vline(xintercept = theta_star_list$sigma.p, linetype = 2) +
        theme(legend.position = "none") 
g5 <- ggplot(data = abc_df, aes(x = samples.tau, colour = methods, fill = methods)) +
        geom_density(alpha = 0.5) +
        labs(x = "tau") +
        change_sizes(16, 20) +
        scale_color_manual(name = "", values = my_colours) +
        scale_fill_manual(name = "", values = my_colours) +
        # geom_vline(xintercept = theta_star_list$tau, linetype = 2) +
        theme(legend.position = "none") 
g6 <- ggplot(data = abc_df, aes(x = samples.delta, colour = methods, fill = methods)) +
        geom_density(alpha = 0.5) +
        labs(x = "delta") +
        change_sizes(16, 20) +
        scale_color_manual(name = "", values = my_colours) +
        scale_fill_manual(name = "", values = my_colours) +
        # geom_vline(xintercept = theta_star_list$delta, linetype = 2) +
        theme(legend.position = "none") 
gridExtra::grid.arrange(g1, g2, g3, g4, g5, g6, ncol = 3)
dev.off()

# pot realizations
obs_df <- data.frame(time = blowfly[, 2], y = blowfly[, 1], lab = "real")
posterior_means <- matrix(NA, nrow = length(method_names), ncol = length(theta_names))
realizations_colours <- c(real = "grey50", "1" = "chocolate", "2" = "purple", "3" = "darkblue")
for (i in 1:(length(method_names) - 1)){
  posterior_means[i, ] <- colMeans(filter(abc_df, methods == method_names[i])[, 2:7])
  obsfake_df <- rbind(obs_df,
                      data.frame(time = blowfly[, 2], y = simulate(posterior_means[i, ]), lab = "1")
                     )
  plt <- ggplot(obsfake_df, aes(x = time, y = y, colour = lab)) +
          geom_line(size = 1.2) +
          scale_color_manual(name = "", values = realizations_colours) +
          labs(x = "time", y = "number of blowflies") +
          change_sizes(16, 20) +
          theme(legend.position = "none") 
  ggsave(plt, file = paste0(plotprefix, "realizations_", method_names[i], ".pdf"), width = 14)
}




