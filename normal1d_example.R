# Normal(theta_\ast, \sigma_0^2) with Normal(m_0, \tau_0^2) prior
# Only theta_\ast is unknown
library(ggplot2)
set.seed(2020)
theta_star <- list(theta = 2, sigma = 2)
hyperparams <- list(m = 1, tau = 2)
epsilon <- c(0.1, 0.5, 1, 2.5, 5)   # c(0.05, 0.1, 0.5, 1)
nobservation <- 500
nthetas <- 1024
y <- matrix(rnorm(nobservation, mean = theta_star$theta, sd = theta_star$sigma), ncol = 1)


# Gamma prior
rprior <- function(n){ return(rnorm(n, mean = hyperparams$m, sd = hyperparams$tau)) }
# generating process
simulate <- function(n, theta){
  return(matrix(rnorm(n = n, mean = theta, sd = theta_star$sigma), ncol = 1))
}

# initialize dataframe to store the true density
thetavals <- seq(-4, 4, length.out = nthetas)
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
abc_post_sd <- 1 / sqrt( hyperparams$tau^(-2) + (theta_star$sigma^2 / nobservation + epsilon[1]^2)^(-1) )
abc_post_mean <- (
              hyperparams$m * hyperparams$tau^(-2) + 
              mean(y) * (theta_star$sigma^2 / nobservation + epsilon[1]^2)^(-1)
             ) * post_sd^2

abcposterior_func <- function(theta){
  return(dnorm(theta, mean = abc_post_mean, sd = abc_post_sd))
}

posterior_df$abc_posterior <- abcposterior_func(thetavals)


# initialize arguments
source("src/rej_abc.R")
source("src/soft_abc.R")
args <- list(nthetas = nthetas, y = y,
              rpiror = rprior,
              simulate = simulate,
              kernel = "gaussian",
              discrepancy = l2norm,
              sumstat = mean
              )

# initialize dataframes to store samples
method_names <- c(paste("epsilon =", epsilon), "gaussian")
abc_df <- data.frame(methods = rep(method_names, each = nthetas),
                     samples = NA
                    )

# rejection abc
for (i in 1:length(method_names)){
  if (i <= length(epsilon)){
    # samples_df <- rej_abc(N = nthetas, epsilon = epsilon[i], y = y, rprior = rprior, simulate = simulate, sumstat = "mean")
    args_rej <- args
    args_rej$epsilon <- epsilon[i]
    args_rej$kernel <- "uniform"
    samples_df <- rej_abc(args_rej)
  } else {
    args$epsilon <- epsilon[i]
    samples_df <- soft_abc(args)
  }
  abc_df$samples[(1 + (i - 1) * nthetas): (i * nthetas)] <- samples_df$samples
}


# plot results
plt_color <- scales::seq_gradient_pal(rgb(1, 0.5, 0.5), "darkblue")(seq(0, 1, length.out = length(method_names)))
plt <- ggplot(abc_df) +
        geom_density(aes(x = samples, colour = methods)) +
        scale_color_manual(values = plt_color) +
        geom_line(data = posterior_df, aes(x = thetavals, y = true_posterior)) +
        geom_vline(xintercept = theta_star$theta, linetype = "dashed") +
        labs(x = "theta", y = "density") +
        theme(
              legend.position = c(.95, .95),
              legend.justification = c("right", "top"),
              legend.title=element_blank()
             ) +
        guides(color = guide_legend(override.aes = list(linetype = "solid")))

ggsave(plt, file = "plots/soft_abc/normal1d_eg.pdf", height = 5)
