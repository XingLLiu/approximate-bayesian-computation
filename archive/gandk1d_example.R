# g-and-k distribution with uniform prior on [0, 10]^4
source("src/mcmc/mh.R")
source("src/gandk/gandk.R")
source("src/rej_abc.R")
set.seed(2020)
theta_star <- list(a = 3, b = 1, g = 2, k = 0.5)
theta_star_vec <- c(theta_star$a, theta_star$b, theta_star$g, theta_star$k)
hyperparams <- list(burnin = 50000, thinning = 10)
epsilon <- c(0.1)
nobservation <- 1000
nthetas <- 1024
resultsprefix <- "results/gandk1d/"
plotprefix <- "plots/gandk1d/"

# prior is uniform [0,10] on each parameter
rprior <- function(N, parameters){
  return(matrix(runif(N * 4, min = 0, max = 10), ncol = 4))
}
# evaluate the log-density of the prior, for each particle
dprior <- function(thetaparticles, parameters){
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

target <- list(rprior = rprior,
               dprior = dprior,
               simulate = simulate,
               loglikelihood = loglikelihood,
               parameter_names = c("a", "b", "g", "k"),
               parameters = list(),
               thetadim = 4,
               ydim = 1)

y <- target$simulate(theta_star_vec)

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


# data-generating process
simulate <- function(n, theta){
  return(matrix(qgandk(runif(n), theta), ncol = 1))
}

# summary statistic
sumstat <- function(z){
  s <- rep(0, 4)
  octiles <- quantile(z, probs = (1:7)/10)
  s[1] <- octiles[4]
  s[2] <- octiles[6] - octiles[2]
  s[3] <- (octiles[6] + octiles[2] - 2 * octiles[4]) / s[2]
  s[4] <- (octiles[7] - octiles[5] + octiles[3] - octiles[1]) / s[2]
  return(s)
}

# initialize dataframes to store samples
method_names <- c(paste("epsilon =", epsilon, " uniform"), "MMD", "Wasserstein", "KL Divergence")
abc_df <- data.frame(methods = rep(method_names, each = nthetas),
                     samples.a = NA,
                     samples.b = NA,
                     samples.g = NA,
                     samples.k = NA
                    )

# rej ABC
args_rej <- list(nthetas = nthetas, y = y,
                  rpiror = rprior,
                  simulate = simulate,
                  kernel = "uniform",
                  discrepancy = l2norm,
                  sumstat = sumstat,
                  epsilon = 0.5 #epsilon[1]
                )

samples_df <- rej_abc(args_rej)
abc_df[index(1, nthetas), 2:5] <- samples_df$samples


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
                  epsilon = 0.05
                )

k2abc_out <- rej_abc(args_mmd)
abc_df[index(2, nthetas), 2:5] <- k2abc_out$samples


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
                  epsilon = 0.5
                 ) 

wabc_out <- rej_abc(args_wabc)
abc_df[index(3, nthetas), 2:5] <- wabc_out$samples


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
                epsilon = 0.3
               ) 

klabc_out <- rej_abc(args_kl)
abc_df[index(4, nthetas), 2:5] <- klabc_out$samples


# save results
write.csv(abc_df, paste0(resultsprefix, "abc_df.csv"), row.names = FALSE)
abc_df <- read.csv(paste0(resultsprefix, "abc_df.csv"))

# plot results
pdf(paste0(plotprefix, "gandk1d_eg.pdf"))
g1 <- ggplot(mhout.df, aes(x = X.1)) + 
        geom_density() + 
        geom_vline(xintercept = theta_star_vec[1], linetype = 2)
g2 <- ggplot(mhout.df, aes(x = X.2)) + 
        geom_density() + 
        geom_vline(xintercept = theta_star_vec[2], linetype = 2)
g3 <- ggplot(mhout.df, aes(x = X.3)) + 
        geom_density() + 
        geom_vline(xintercept = theta_star_vec[3], linetype = 2)
g4 <- ggplot(mhout.df, aes(x = X.4)) + 
        geom_density() + 
        geom_vline(xintercept = theta_star_vec[4], linetype = 2)

g1 <- g1 + geom_density(data = abc_df, aes(x = samples.a, colour = methods)) +
          labs(x = "a") +
          xlim(1.5, 5) +
          theme(legend.position = "none")
g2 <- g2 + geom_density(data = abc_df, aes(x = samples.b, colour = methods)) +
          labs(x = "b") +
          xlim(0, 5) +
          theme(
            legend.position = c(.95, .95),
            legend.justification = c("right", "top"),
            legend.title = element_blank()
          )
g3 <- g3 + geom_density(data = abc_df, aes(x = samples.g, colour = methods)) +
          labs(x = "g") +
          xlim(0, 5) +
          theme(legend.position = "none") 
g4 <- g4 + geom_density(data = abc_df, aes(x = samples.k, colour = methods)) +
          labs(x = "k") +
          xlim(0, 4) +
          theme(legend.position = "none") 

gridExtra::grid.arrange(g1, g2, g3, g4, ncol = 2)
dev.off()

