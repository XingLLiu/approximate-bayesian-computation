# Exponential with Gamma prior
library(ggplot2)
set.seed(2020)
theta_star <- list(lambda = 2)
hyperparams <- list(alpha = 2, beta = 2)
epsilon <- c(0.01, 0.05, 0.1, 0.5, 1, 2)
nobservation <- 1
nthetas <- 1024
y <- matrix(rexp(nobservation, theta_star$lambda), ncol = 1)


# Gamma prior
rprior <- function(n){ return(rgamma(n, shape = hyperparams$alpha, rate = hyperparams$beta)) }
# generating process
simulate <- function(n, theta){ return(matrix(rexp(n, rate = theta), ncol = 1)) }

# initialize dataframe to store the true density
thetavals <- seq(0, 7.5, length.out = nthetas)
posterior_df <- data.frame(thetavals = thetavals,
                           true_posterior = NA,
                           abc_posterior = NA
                          )
# true posterior
posterior_df$true_posterior <- dgamma(thetavals,
                                      shape = hyperparams$alpha + nobservation,
                                      rate = hyperparams$beta + sum(y))

# exact abc posterior (for n = 1)
abcposterior_func <- function(theta, epsilonval){
  const <- gamma(hyperparams$alpha) * (1 / (hyperparams$beta + c(y) - epsilonval)^hyperparams$alpha - 
                                       1 / (hyperparams$beta + c(y) + epsilonval)^hyperparams$alpha
                                      )
  return(
          theta^(hyperparams$alpha - 1) * exp(- (hyperparams$beta + c(y)) * theta) * 
                ( exp(epsilonval * theta) - exp(- epsilonval * theta) ) * (theta > 0) / const 
        )
}

posterior_df$abc_posterior <- abcposterior_func(thetavals, epsilon[1])



# initialize arguments
source("src/rej_abc.R")
args <- list(nthetas = nthetas, y = y,
              rpiror = rprior,
              simulate = simulate,
              kernel = "uniform",
              discrepancy = l2norm
              )

# initialize dataframes to store samples
method_names <- paste("epsilon =", epsilon)
abc_df <- data.frame(methods = rep(method_names, each = nthetas),
                     samples = NA
                    )
posterior_df2 <- data.frame(methods = rep(method_names, each = nthetas),
                     thetas = rep(thetavals, length(method_names)),
                     densities = NA
                    )

# rejection abc
for (i in 1:length(epsilon)){
  # samples_df <- rej_abc(N = nthetas, epsilon = epsilon[i], y = y, rprior = rprior, simulate = simulate)
  args$epsilon <- epsilon[i]
  samples_df <- rej_abc(args)
  abc_df$samples[(1 + (i - 1) * nthetas): (i * nthetas)] <- samples_df$samples
  # exact abc posterior densities
  posterior_df2$densities[(1 + (i - 1) * nthetas): (i * nthetas)] <- abcposterior_func(thetavals, epsilon[i])
}


# plot results
plt_color <- init_discrete_grad_colours(length(epsilon))
plt <- ggplot(posterior_df) +
        geom_line(aes(x = thetavals, y = true_posterior)) +
        geom_ribbon(aes(x = thetavals, ymax = true_posterior), ymin = 0, alpha = 0.5) +
        geom_density(data = abc_df, aes(x = samples, colour = methods)) +
        scale_color_manual(values = plt_color) +
        geom_vline(xintercept = theta_star$lambda, linetype = "dashed") +
        labs(x = "theta", y = "density") +
        xlim(0, 7.5) +
        change_sizes(16, 20) +
        add_legend(0.95, 0.95)

ggsave(plt, file = "plots/soft_abc/exponential_eg.pdf", height = 5)

# exact abc posterior densities
plt <- ggplot(posterior_df) +
        geom_line(aes(x = thetavals, y = true_posterior)) +
        geom_ribbon(aes(x = thetavals, ymax = true_posterior), ymin = 0, alpha = 0.5) +
        geom_line(data = posterior_df2, aes(x = thetas, y = densities, colour = methods)) +
        scale_color_manual(values = plt_color) +
        geom_vline(xintercept = theta_star$lambda, linetype = "dashed") +
        labs(x = "theta", y = "density") +
        xlim(0, 7.5) +
        change_sizes(16, 20) +
        add_legend(0.95, 0.95)

ggsave(plt, file = "plots/soft_abc/exponential_eg_exact.pdf", height = 5)
