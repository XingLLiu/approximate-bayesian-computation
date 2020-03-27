# Normal(theta_\ast, \sigma_0^2) with Normal(m_0, \tau_0^2) prior
# Only theta_\ast is unknown
library(ggplot2)
library(winference)
source("src/rej_abc.R")
source("src/soft_abc.R")
set.seed(2020)
theta_star <- list(theta = c(0.3, 0.7, 0.7, -0.7, -0.7),
                   cov1 = matrix(c(0.5, -0.3, -0.3, 0.5), ncol = 2),
                   cov2 = matrix(c(0.25, 0, 0, 0.25), ncol = 2)
                  )
# hyperparams <- list(m = 1, tau = 2)
epsilon <- c(1)   # c(0.05, 0.1, 0.5, 1)
nobservation <- 500
nthetas <- 1024


# Gamma prior
rprior <- function(n){
  return(c(runif(1), runif(4, -1, 1)))
}
# generating process
simulate <- function(n, theta){
  u <- runif(n)
  ind <- u < theta[1]
  z <- matrix(NA, ncol = 2, nrow = n)
  z[!ind] <- fast_rmvnorm(sum(ind), mean = theta[2:3], covariance = theta_star$cov1)
  z[ind] <- fast_rmvnorm(sum(!ind), mean = theta[4:5], covariance = theta_star$cov2)
  return(z)
}

y <- simulate(nobservation, theta_star$theta)

# summary statistic
sumstat <- function(z){
  s <- rep(NA, 5)
  s[1:2] <- colMeans(z)
  # s[3:5] <- as.vector(cov(z))[c(1, 4, 2)]
  s[3:5] <- as.vector( t(y) %*% y / nrow(y) )[c(1, 4, 2)]
  return(s)
}


# initialize dataframes to store samples
method_names <- c(paste("Rejection"), "MMD", "Wasserstein", "KL divergence")
abc_df <- data.frame(methods = rep(method_names, each = nthetas),
                     samples.p = NA,
                     samples.mu01 = NA,
                     samples.mu02 = NA,
                     samples.mu11 = NA,
                     samples.mu12 = NA
                    )

# rej ABC
args_rej <- list(nthetas = nthetas, y = y,
                  rpiror = rprior,
                  simulate = simulate,
                  kernel = "uniform",
                  discrepancy = l2norm,
                  sumstat = sumstat,
                  epsilon = 0.1
                )

samples_df <- rej_abc(args_rej)
abc_df[index(1, nthetas), 2:ncol(abc_df)] <- samples_df$samples


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
                  epsilon = 0.25
                )

k2abc_out <- rej_abc(args_mmd)
abc_df[index(2, nthetas), 2:ncol(abc_df)] <- k2abc_out$samples


# WABC
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
abc_df[index(3, nthetas), 2:ncol(abc_df)] <- wabc_out$samples


# KL ABC
kldist <- function(y, z){
  return(FNN::KLx.divergence(y, z, k = 1))
  # return(KLx.divergence(y, z, k = 1))
} 

args_kl <- list(nthetas = nthetas, y = y,
                rpiror = rprior,
                simulate = simulate,
                kernel = "uniform",
                discrepancy = kldist,
                epsilon = 0.5
               ) 

klabc_out <- rej_abc(args_kl)
abc_df[index(4, nthetas), 2:ncol(abc_df)] <- klabc_out$samples


# save results
write.csv(abc_df, "results/normal2d/abc_df.csv", row.names = FALSE)
abc_df <- read.csv("results/normal2d/abc_df.csv")


# plot results
pdf("plots/normal2d/posterior_densities.pdf", width = 14)
g1 <- ggplot(data = abc_df, aes(x = samples.p, colour = methods)) +
        geom_density() +
        labs(x = "p") +
        geom_vline(xintercept = theta_star$theta[1], linetype = 2) +
        theme(legend.position = "none") 
g2 <- ggplot(data = abc_df, aes(x = samples.mu01, colour = methods)) +
        geom_density() +
        labs(x = "mu01") +
        geom_vline(xintercept = theta_star$theta[2], linetype = 2) +
        theme(legend.position = "none") 
g3 <- ggplot(data = abc_df, aes(x = samples.mu02, colour = methods)) +
        geom_density() +
        labs(x = "mu02") +
        geom_vline(xintercept = theta_star$theta[3], linetype = 2) +
        theme(legend.position = "none") 
g4 <- ggplot(data = abc_df, aes(x = samples.mu11, colour = methods)) +
        geom_density() +
        labs(x = "mu11") +
        geom_vline(xintercept = theta_star$theta[4], linetype = 2) +
        theme(legend.position = "none") 
g5 <- ggplot(data = abc_df, aes(x = samples.mu12, colour = methods)) +
        geom_density() +
        labs(x = "mu12") +
        geom_vline(xintercept = theta_star$theta[5], linetype = 2) +
        theme(legend.position = "none") 
g6 <- ggplot(data = abc_df, aes(x = samples.mu01, colour = methods)) +
        geom_density() +
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
  pdf(paste0("plots/normal2d/contour_", method, ".pdf"), width = 14)
  g1 <- ggplot(filter(abc_df, methods == method),
              aes(x = samples.mu01, y = samples.mu02)
              ) +
          geom_density_2d(aes(color = ..level..), size = 1.5) +
          scale_color_viridis_c() +
          labs(x = "mu01") +
          labs(y = "mu02") +
          xlim(-1, 1) +
          ylim(-1, 1) +
          ggtitle(method) +
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
          geom_vline(xintercept = theta_star$theta[4], linetype = 2) +
          geom_hline(yintercept = theta_star$theta[5], linetype = 2) +
          theme(legend.position = "none") 
  gridExtra::grid.arrange(g1, g2, ncol = 2)
  dev.off()
}

for (method in method_names){
  plot_and_save_contour(method)
}




