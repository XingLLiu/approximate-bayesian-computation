# g-and-k distribution with uniform prior on [0, 10]^4
source("src/mcmc/mh.R")
source("src/gandk/gandk.R")
set.seed(2020)
theta_star <- list(a = 3, b = 1, g = 2, k = 0.5)
theta_star_vec <- c(theta_star$a, theta_star$b, theta_star$g, theta_star$k)
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
  observations <- qgandk(runif(nobservation), theta)
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
# generating process
simulate_abc <- function(n, theta){
  observations <- qgandk(runif(n), theta)
  return(observations)
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
method_names <- c(paste("epsilon =", epsilon), "gaussian")
abc_df <- data.frame(method = rep(method_names, each = nthetas),
                     samples.a = NA,
                     samples.b = NA,
                     samples.g = NA,
                     samples.k = NA
                    )

# rejection abc
source("src/rej_abc.R")
for (i in 1:length(epsilon)){
  samples_df <- rej_abc(N = nthetas, epsilon = epsilon[i], y = y, prior = rprior,
                        simulate = simulate_abc, sumstat = sumstat)
  abc_df[(1 + (i - 1) * nthetas): (i * nthetas), -1] <- samples_df$samples
}

# soft abc with gaussian kernel
source("src/soft_abc.R")
samples_df <- soft_abc(N = nthetas, epsilon = epsilon[i], y = y, prior = rprior,
                        simulate = simulate_abc, sumstat = sumstat)
ind <- (1 + length(epsilon) * nthetas): ((length(epsilon) + 1) * nthetas)
abc_df[ind, -1] <- apply(samples_df$samples, 2, sample,
                          size = nthetas, prob = samples_df$weights, replace = TRUE)


# plot results
# plt_color <- scales::seq_gradient_pal(rgb(1, 0.5, 0.5), "darkblue")(seq(0, 1, length.out = length(epsilon) + 1))
g1 + geom_density(data = abc_df, aes(x = samples.a, colour = method))

g2 + geom_density(data = abc_df, aes(x = samples.b, colour = method))

g3 + geom_density(data = abc_df, aes(x = samples.g, colour = method))

g4 + geom_density(data = abc_df, aes(x = samples.k, colour = method))



plt <- gridExtra::grid.arrange(g1, g2, g3, g4, ncol = 2)

ggsave(plt, file = "plots/soft_abc/gandk1d_eg.pdf", height = 5)

