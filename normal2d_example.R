# Normal(theta_\ast, \sigma_0^2) with Normal(m_0, \tau_0^2) prior
# Only theta_\ast is unknown
library(ggplot2)
library(winference)
source("src/rej_abc.R")
source("src/soft_abc.R")
source("src/sabc.R")
set.seed(2020)
theta_star <- list(theta = c(0.3, 0.7, 0.7, -0.7, -0.7),
                   cov1 = matrix(c(0.5, -0.3, -0.3, 0.5), ncol = 2),
                   cov2 = matrix(c(0.25, 0, 0, 0.25), ncol = 2)
                  )
theta_names <- c("p", "mu01", "mu02", "mu11", "mu12")
theta_star_list <- list(p = theta_star$theta[1],
                          mu01 = theta_star$theta[2],
                          mu02 = theta_star$theta[3],
                          mu11 = theta_star$theta[4],
                          mu12 = theta_star$theta[5]
                        )
# hyperparams <- list(m = 1, tau = 2)
epsilon <- c(1)   # c(0.05, 0.1, 0.5, 1)
nobservation <- 500
nthetas <- 1024
maxsimulation <- 10^5
resultsprefix <- "results/normal2d/"
plotprefix <- "plots/normal2d/"


# rprior
rprior <- function(n, hyperparams){
  return(matrix(c(runif(n), runif(4 * n, -1, 1)), nrow = n))
}
# log dprior
dprior <- function(theta, hyperparams){
  eval <- dunif(theta[, 1], log = TRUE) +
            dunif(theta[, 2], min = -1, max = 1, log = TRUE) +
            dunif(theta[, 3], min = -1, max = 1, log = TRUE) +
            dunif(theta[, 4], min = -1, max = 1, log = TRUE) +
            dunif(theta[, 5], min = -1, max = 1, log = TRUE)
  return(eval)
}
# generating process
simulate <- function(theta){
  u <- runif(nobservation)
  ind <- u < theta[1]
  z <- matrix(NA, ncol = 2, nrow = nobservation)
  z[!ind] <- fast_rmvnorm(sum(!ind), mean = theta[2:3], covariance = theta_star$cov1)
  z[ind] <- fast_rmvnorm(sum(ind), mean = theta[4:5], covariance = theta_star$cov2)
  return(z)
}

y <- simulate(theta_star$theta)

# summary statistic
sumstat <- function(z){
  s <- rep(NA, 5)
  s[1:2] <- colMeans(z)
  s[3:5] <- as.vector( t(y) %*% y / nrow(y) )[c(1, 4, 2)]
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
                     samples.p = NA,
                     samples.mu01 = NA,
                     samples.mu02 = NA,
                     samples.mu11 = NA,
                     samples.mu12 = NA
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
wdistance <- function(y_fake){
  sink("/dev/null")
  dist <- exact_transport_distance(t(y), t(y_fake))
  sink()
  return(dist)
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
abc_df <- read.csv(paste0(resultsprefix, "abc_df.csv"))


# plot results
my_colours <- init_colours()
pdf(paste0(plotprefix, "posterior_densities.pdf"), width = 14)
g1 <- ggplot(data = abc_df, aes(x = samples.p, colour = methods, fill = methods)) +
        geom_density(alpha = 0.5) +
        labs(x = "p") +
        change_sizes(16, 20) +
        scale_color_manual(name = "", values = my_colours) +
        scale_fill_manual(name = "", values = my_colours) +
        geom_vline(xintercept = theta_star$theta[1], linetype = 2) +
        theme(legend.position = "none") 
g2 <- ggplot(data = abc_df, aes(x = samples.mu01, colour = methods, fill = methods)) +
        geom_density(alpha = 0.5) +
        labs(x = "mu01") +
        change_sizes(16, 20) +
        scale_color_manual(name = "", values = my_colours) +
        scale_fill_manual(name = "", values = my_colours) +
        geom_vline(xintercept = theta_star$theta[2], linetype = 2) +
        theme(legend.position = "none") 
g3 <- ggplot(data = abc_df, aes(x = samples.mu02, colour = methods, fill = methods)) +
        geom_density(alpha = 0.5) +
        labs(x = "mu02") +
        change_sizes(16, 20) +
        scale_color_manual(name = "", values = my_colours) +
        scale_fill_manual(name = "", values = my_colours) +
        geom_vline(xintercept = theta_star$theta[3], linetype = 2) +
        theme(legend.position = "none") 
g4 <- ggplot(data = abc_df, aes(x = samples.mu11, colour = methods, fill = methods)) +
        geom_density(alpha = 0.5) +
        labs(x = "mu11") +
        change_sizes(16, 20) +
        scale_color_manual(name = "", values = my_colours) +
        scale_fill_manual(name = "", values = my_colours) +
        geom_vline(xintercept = theta_star$theta[4], linetype = 2) +
        theme(legend.position = "none") 
g5 <- ggplot(data = abc_df, aes(x = samples.mu12, colour = methods, fill = methods)) +
        geom_density(alpha = 0.5) +
        labs(x = "mu12") +
        change_sizes(16, 20) +
        scale_color_manual(name = "", values = my_colours) +
        scale_fill_manual(name = "", values = my_colours) +
        geom_vline(xintercept = theta_star$theta[5], linetype = 2) +
        theme(legend.position = "none") 
g6 <- ggplot(data = abc_df, aes(x = samples.mu01, colour = methods, fill = methods)) +
        geom_density(alpha = 0.5) +
        scale_color_manual(name = "", values = my_colours) +
        scale_fill_manual(name = "", values = my_colours) +
        theme(
          legend.position = c(.5, .6),
          legend.justification = c("center"),
          legend.title = element_blank(),
          legend.key.size = unit(2,"line"),
          legend.text = element_text(size = 20)
        )
gridExtra::grid.arrange(g1, g2, g3, g4, g5, g_legend(g6), ncol = 3)
dev.off()

# plot contours
plot_and_save_contour <- function(method){
  pdf(paste0(plotprefix, "contour_", method, ".pdf"), width = 14)
  g1 <- ggplot(filter(abc_df, methods == method),
              aes(x = samples.mu01, y = samples.mu02)
              ) +
          geom_density_2d(aes(color = ..level..), size = 1.5) +
          scale_color_viridis_c() +
          labs(x = "mu01") +
          labs(y = "mu02") +
          xlim(-1, 1) +
          ylim(-1, 1) +
          change_sizes(16, 20) +
          geom_vline(xintercept = theta_star$theta[2], linetype = 2) +
          geom_hline(yintercept = theta_star$theta[3], linetype = 2) +
          theme(legend.position = "none") 
  g2 <- ggplot(filter(abc_df, methods == method),
              aes(x = samples.mu11, y = samples.mu12)
              ) +
          geom_density_2d(aes(color = ..level..), size = 1.5) +
          scale_color_viridis_c() +
          labs(x = "mu11") +
          labs(y = "mu12") +
          xlim(-1, 1) +
          ylim(-1, 1) +
          change_sizes(16, 20) +
          geom_vline(xintercept = theta_star$theta[4], linetype = 2) +
          geom_hline(yintercept = theta_star$theta[5], linetype = 2) +
          theme(legend.position = "none") 
  gridExtra::grid.arrange(g1, g2, ncol = 2)
  dev.off()
}

for (method in method_names){
  plot_and_save_contour(method)
}

# plot data
pdf(paste0(plotprefix, "data.pdf"), width = 14)
g1 <- ggplot(data.frame(y1 = y[, 1], y2 = y[, 2]),
              aes(x = y1, y = y2)
            ) +
      geom_point(size = 1.5) +
      scale_color_manual(values = "grey") +
      labs(x = "y1") +
      labs(y = "y2") +
      xlim(-1, 1) +
      ylim(-1, 1) +
      change_sizes(16, 20) +
      geom_vline(xintercept = theta_star$theta[2], linetype = 2) +
      geom_hline(yintercept = theta_star$theta[3], linetype = 2) +
      geom_vline(xintercept = theta_star$theta[4], linetype = 2) +
      geom_hline(yintercept = theta_star$theta[5], linetype = 2) +
      theme(legend.position = "none") 
dev.off()


# computational times
ztemp <- simulate(theta_star$theta)
print("Computational times of one evaluation:")
microbenchmark::microbenchmark(
                                eucdiscrep(ztemp),
                                mmdsq(ztemp),
                                wdistance(ztemp),
                                kldist(ztemp)
                              )
