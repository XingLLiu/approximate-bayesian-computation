# g-and-k distribution with uniform prior on [0, 10]^4
source("src/mcmc/mh.R")
source("src/gandk/gandk.R")
set.seed(2020)
theta_star <- list(a = 3, b = 1, g = 2, k = 0.5)
hyperparams <- list(burnin = 50000, thinning = 10)
epsilon <- c(1)
nobservation <- 250
nthetas <- 1024

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
simulate <- function(theta){
  observations <- qgandk(runif(nobservation), )
  return(observations)
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

theta_star_vec <- c(theta_star$a, theta_star$b, theta_star$g, theta_star$k)
y <- target$simulate(theta_star_vec)

# tuning params for MCMC
tuning_parameters <- list(niterations = (nthetas - 1) * hyperparams$thinning + hyperparams$burnin + 1, 
                          nchains = 1,
                          cov_proposal = diag(1e-2, 4),
                          adaptation = 2000,
                          init_chains = matrix(theta_star_vec, nrow = 1)
                          )
mhout <- mh(y, target, tuning_parameters)
mhout.df <- mhchainlist_to_dataframe(mhout$chains)
mhout.df <- mhout.df %>% filter(iteration > hyperparams$burnin,
                                iteration %% hyperparams$thinning == 1)

write.csv(mhout.df, "results/gandk/mcmcsamples.R")

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
gridExtra::grid.arrange(g1, g2, g3, g4, ncol = 2)



# rejection abc
# Gamma prior
rprior <- function(n){ return(rnorm(4 * n, mean = hyperparams$m, sd = hyperparams$tau)) }

# generating process
simulate <- function(n, theta){
  observations <- qgandk(runif(n), )
  return(observations)
}

# initialize dataframes to store samples
method_names <- paste("epsilon =", epsilon)
abc_df <- data.frame(method = rep(method_names, each = nthetas),
                     samples = NA
                    )

# rejection abc
source("src/abc_rej.R")
for (i in 1:length(epsilon)){
  samples_df <- soft_abc(N = nthetas, epsilon = epsilon[i], y = y, prior = rprior,
                          simulate = simulate, sumstat = "mean")
  abc_df$samples[(1 + (i - 1) * nthetas): (i * nthetas)] <- samples_df$samples
}



