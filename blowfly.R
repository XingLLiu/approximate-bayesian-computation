# blowfly data
library(winference)
library(ggplot2)
source("src/rej_abc.R")
source("src/soft_abc.R")
source("src/sabc.R")
source("src/blowfly/blowfly.R")
set.seed(2020)
theta_names <- c("P", "delta", "N0", "sigma.d", "sigma.p", "tau")
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
save(rej_out, file = paste0(resultsprefix, "rej_out.RData"))


# K2 ABC
source("src/mmd/mmdsq_c.R")
bandwidth <- median(dist(y, method = "manhattan"))
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
save(mmd_out, file = paste0(resultsprefix, "mmd_out.RData"))


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
save(wabc_out, file = paste0(resultsprefix, "wabc_out.RData"))


# # KL ABC
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

klabc_out <- sabc(args_kl, maxsimulation = maxsimulation)
abc_df[index(4, nthetas), 2:ncol(abc_df)] <- sabc_get_last_samples(klabc_out)[, theta_names]
save(klabc_out, file = paste0(resultsprefix, "klabc_out.RData"))


# save results
write.csv(abc_df, paste0(resultsprefix, "abc_df.csv"), row.names = FALSE)
abc_df <- read.csv(paste0(resultsprefix, "/abc_df.csv"))


# plot results
my_colours <- init_colours()
pdf(paste0(plotprefix, "posterior_densities.pdf"), width = 18)
g1 <- ggplot(data = abc_df, aes(x = samples.P, colour = methods, fill = methods)) +
        geom_density(alpha = 0.5) +
        labs(x = "P") +
        change_sizes(16, 20) +
        scale_color_manual(name = "", values = my_colours) +
        scale_fill_manual(name = "", values = my_colours) +
        theme(legend.position = "none") +
        add_legend(0.95, 0.95)
g2 <- ggplot(data = abc_df, aes(x = samples.N0, colour = methods, fill = methods)) +
        geom_density(alpha = 0.5) +
        labs(x = "N0") +
        change_sizes(16, 20) +
        scale_color_manual(name = "", values = my_colours) +
        scale_fill_manual(name = "", values = my_colours) +
        theme(legend.position = "none") 
g3 <- ggplot(data = abc_df, aes(x = samples.sigma.d, colour = methods, fill = methods)) +
        geom_density(alpha = 0.5) +
        labs(x = "sigma_d") +
        change_sizes(16, 20) +
        scale_color_manual(name = "", values = my_colours) +
        scale_fill_manual(name = "", values = my_colours) +
        theme(legend.position = "none") 
g4 <- ggplot(data = abc_df, aes(x = samples.sigma.p, colour = methods, fill = methods)) +
        geom_density(alpha = 0.5) +
        labs(x = "sigma_p") +
        change_sizes(16, 20) +
        scale_color_manual(name = "", values = my_colours) +
        scale_fill_manual(name = "", values = my_colours) +
        theme(legend.position = "none") 
g5 <- ggplot(data = abc_df, aes(x = samples.tau, colour = methods, fill = methods)) +
        geom_density(alpha = 0.5) +
        labs(x = "tau") +
        change_sizes(16, 20) +
        scale_color_manual(name = "", values = my_colours) +
        scale_fill_manual(name = "", values = my_colours) +
        theme(legend.position = "none") 
g6 <- ggplot(data = abc_df, aes(x = samples.delta, colour = methods, fill = methods)) +
        geom_density(alpha = 0.5) +
        labs(x = "delta") +
        change_sizes(16, 20) +
        scale_color_manual(name = "", values = my_colours) +
        scale_fill_manual(name = "", values = my_colours) +
        theme(legend.position = "none") 
gridExtra::grid.arrange(g1, g2, g3, g4, g5, g6, ncol = 3)
dev.off()

# plot realizations
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
          change_sizes(28, 36) +
          theme(legend.position = "none") 
  ggsave(plt, file = paste0(plotprefix, "realizations_", method_names[i], ".pdf"), width = 14)
}


for (i in 2:(length(method_names) - 1)){
  # posterior_means[i, ] <- colMeans(filter(abc_df, methods == method_names[i])[, 2:7]) c(8, 0.2, 14, 0.8, 0.3, 7)
  obsfake_df <- rbind(obs_df,
                      data.frame(time = blowfly[, 2], y = simulate(posterior_means[1, ]), lab = "1")
                     )
  plt <- ggplot(obsfake_df, aes(x = time, y = y, colour = lab)) +
          geom_line(size = 1.2) +
          scale_color_manual(name = "", values = realizations_colours) +
          labs(x = "time", y = "number of blowflies") +
          change_sizes(16, 20) +
          theme(legend.position = "none") 
  ggsave(plt, file = paste0(plotprefix, "haha_realizations_", method_names[i], ".pdf"), width = 14)
}


# thresholds
threshold_history <- rbind(
                            cbind(method_names[1], cumsum(rej_out$ncomputed), rej_out$threshold_history),
                            cbind(method_names[2], cumsum(mmd_out$ncomputed), mmd_out$threshold_history),
                            cbind(method_names[3], cumsum(wabc_out$ncomputed), wabc_out$threshold_history)
                          )
threshold_history <- data.frame(
                                methods = threshold_history[, 1],
                                nsimulations = as.numeric(threshold_history[, 2]),
                                thresholds = as.numeric(threshold_history[, 3])
                               )

draw_thresholds <- function(method){
  g1 <- ggplot(data = threshold_history %>% filter(methods == method), 
               aes(y = thresholds, x = nsimulations)
              ) +
        geom_line(color = my_colours[method]) +
        geom_point(color = my_colours[method]) +
        labs(x = "number of model simulations", y = "threshold") +
        xlim(0, maxsimulation * 1.2) +
        change_sizes(16, 20) +
        add_legend(0.95, 0.95)
  return(g1)
}

pdf(paste0(plotprefix, "thresholds.pdf"), height = 5, width = 24)
g1 <- draw_thresholds(method_names[1])
g2 <- draw_thresholds(method_names[2])
g3 <- draw_thresholds(method_names[3])
gridExtra::grid.arrange(g1, g2, g3, nrow = 1)
dev.off()


# computational times
ztemp <- simulate(colMeans(filter(abc_df, methods == "Wasserstein")[, 2:7]))
print("Computational times of one evaluation:")
microbenchmark::microbenchmark(
                                eucdiscrep(ztemp),
                                mmdsq(ztemp),
                                wdistance(ztemp),
                                times = 1000
                              )


