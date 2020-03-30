# g-and-k distribution with uniform prior on [0, 10]^4
source("src/mcmc/mh.R")
source("src/gandk/gandk.R")
source("src/rej_abc.R")
source("src/sabc.R")
set.seed(2020)
theta_star <- list(a = 3, b = 1, g = 2, k = 0.5)
theta_star_vec <- c(theta_star$a, theta_star$b, theta_star$g, theta_star$k)
theta_names <- c("a", "b", "g", "k")
theta_star_list <- list(a = theta_star$theta[1],
                          b = theta_star$theta[2],
                          g = theta_star$theta[3],
                          k = theta_star$theta[4]
                        )
hyperparams <- list(burnin = 50000, thinning = 10)
nobservation <- 1000
nthetas <- 1024
maxsimulation <- 10^5
resultsprefix <- "results/gandk1d/"
plotprefix <- "plots/gandk1d/"

# prior is uniform [0,10] on each parameter
rprior <- function(N, hyperparams){
  return(matrix(runif(N * 4, min = 0, max = 10), ncol = 4))
}
# log dprior
dprior <- function(thetaparticles, hyperparams){
  densities <- rep(0, nrow(thetaparticles))
  for (i in 1:nrow(thetaparticles)){
    if (any(thetaparticles[i, ] > 10) || any(thetaparticles[i, ] < 0)){
      densities[i] <- -Inf
    }
  }
  return(densities)
}

# function to generate a dataset given a parameter
# simulate <- function(theta){
#   observations <- qgandk(runif(nobservation), theta)
#   return(observations)
# }
simulate <- function(theta){
  return(matrix(qgandk(runif(nobservation), theta), ncol = 1))
}

# loglikelihood of ys, given each theta 
loglikelihood <- function(theta, ys, ...){
  return(logdgandk(ys, theta, ...))
}

y <- simulate(theta_star_vec)

target <- list(rprior = rprior,
               dprior = dprior,
               simulate = simulate,
               loglikelihood = loglikelihood,
               parameter_names = c("a", "b", "g", "k"),
               parameters = list(),
               thetadim = length(theta_names),
               ydim = ncol(y)
              )


# tuning params for MCMC
source("src/mcmc/mh.R")
tuning_parameters <- list(niterations = (nthetas - 1) * hyperparams$thinning + hyperparams$burnin + 1, 
                          nchains = 1,
                          cov_proposal = diag(1e-2, 4),
                          adaptation = 2000,
                          init_chains = matrix(theta_star_vec, nrow = 1)
                          )
# mhout <- mh(y, target, tuning_parameters)
# mhout.df <- mhchainlist_to_dataframe(mhout$chains)
# mhout.df <- mhout.df %>% dplyr::filter(iteration > hyperparams$burnin,
#                                        iteration %% hyperparams$thinning == 1)

# write.csv(mhout.df, paste0(resultsprefix, "gandk_mcmc.csv"), row.names = FALSE)
mhout.df <- read.csv(paste0(resultsprefix, "gandk_mcmc.csv"), header = TRUE)

# summary statistic
sumstat <- function(z){
  s <- rep(0, 4)
  octiles <- quantile(z, probs = (1:7)/8) # quantile(z, probs = (1:7)/10)
  s[1] <- octiles[4]
  s[2] <- octiles[6] - octiles[2]
  s[3] <- (octiles[6] + octiles[2] - 2 * octiles[4]) / s[2]
  s[4] <- (octiles[7] - octiles[5] + octiles[3] - octiles[1]) / s[2]
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
                     samples.a = NA,
                     samples.b = NA,
                     samples.g = NA,
                     samples.k = NA
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
abc_df <- read.csv(paste0(resultsprefix, "abc_df.csv"))

# plot results
my_colours <- init_colours()
pdf(paste0(plotprefix, "gandk1d_eg.pdf"), width = 14)
g1 <- ggplot(mhout.df, aes(x = X.1, fill = "Posterior", colour = "Posterior"), alpha = 0.5) + 
        geom_density() + 
        geom_vline(xintercept = theta_star_vec[1], linetype = 2)
g2 <- ggplot(mhout.df, aes(x = X.2, fill = "Posterior", colour = "Posterior"), alpha = 0.5) + 
        geom_density() + 
        geom_vline(xintercept = theta_star_vec[2], linetype = 2)
g3 <- ggplot(mhout.df, aes(x = X.3, fill = "Posterior", colour = "Posterior"), alpha = 0.5) + 
        geom_density() + 
        geom_vline(xintercept = theta_star_vec[3], linetype = 2)
g4 <- ggplot(mhout.df, aes(x = X.4, fill = "Posterior", colour = "Posterior"), alpha = 0.5) + 
        geom_density() + 
        geom_vline(xintercept = theta_star_vec[4], linetype = 2)

g1 <- g1 + geom_density(data = abc_df, aes(x = samples.a, fill = methods, colour = methods), alpha = 0.5) +
          scale_color_manual(name = "", values = my_colours) +
          scale_fill_manual(name = "", values = my_colours) +
          labs(x = "a") +
          xlim(2.5, 3.5) +
          change_sizes(16, 20) +
          theme(legend.position = "none")
g2 <- g2 + geom_density(data = abc_df, aes(x = samples.b, fill = methods, colour = methods), alpha = 0.5) +
          scale_color_manual(name = "", values = my_colours) +
          scale_fill_manual(name = "", values = my_colours) +
          labs(x = "b") +
          xlim(0, 2.5) +
          change_sizes(16, 20) +
          add_legend(0.95, 0.95)
g3 <- g3 + geom_density(data = abc_df, aes(x = samples.g, fill = methods, colour = methods), alpha = 0.5) +
          scale_color_manual(name = "", values = my_colours) +
          scale_fill_manual(name = "", values = my_colours) +
          labs(x = "g") +
          xlim(0, 5) +
          change_sizes(16, 20) +
          theme(legend.position = "none") 
g4 <- g4 + geom_density(data = abc_df, aes(x = samples.k, fill = methods, colour = methods), alpha = 0.5) +
          scale_color_manual(name = "", values = my_colours) +
          scale_fill_manual(name = "", values = my_colours) +
          labs(x = "k") +
          xlim(0, 1.5) +
          change_sizes(16, 20) +
          theme(legend.position = "none") 

gridExtra::grid.arrange(g1, g2, g3, g4, ncol = 2)
dev.off()

