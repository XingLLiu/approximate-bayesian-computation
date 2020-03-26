# Exponential with Gamma prior
library(ggplot2)
set.seed(2020)
theta_star <- list(lambda = 2)
hyperparams <- list(alpha = 1, beta = 1)
epsilon <- c(0.01, 0.05, 0.1, 0.5, 1)
nobservation <- 1
nthetas <- 1024
y <- matrix(rexp(nobservation, theta_star$lambda), nocl = 1)


# Gamma prior
rprior <- function(n){ return(rgamma(n, shape = hyperparams$alpha, rate = hyperparams$beta)) }
# generating process
simulate <- function(n, theta){ return(matrix(rexp(n, rate = theta), ncol = 1)) }

# initialize dataframe to store the true density
thetavals <- seq(0, 10, length.out = nthetas)
posterior_df <- data.frame(thetavals = thetavals,
                           true_posterior = NA,
                           abc_posterior = NA
                          )
# true posterior
posterior_df$true_posterior <- dgamma(thetavals,
                                      shape = hyperparams$alpha + nobservation,
                                      rate = hyperparams$beta + sum(y))

# exact abc posterior (for n = 1)
abcposterior_func <- function(theta){
  const <- gamma(hyperparams$alpha) * (1 / (hyperparams$beta + y - epsilon[1])^hyperparams$alpha - 
                                       1 / (hyperparams$beta + y + epsilon[1])^hyperparams$alpha
                                      )
  return(
          theta^(hyperparams$alpha - 1) * exp(- (hyperparams$beta + y) * theta) * 
                ( exp(epsilon[1] * theta) - exp(- epsilon[1] * theta) ) * (theta > 0) / const 
        )
}

posterior_df$abc_posterior <- abcposterior_func(thetavals)

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
abc_df <- data.frame(method = rep(method_names, each = nthetas),
                     samples = NA
                    )

# rejection abc
for (i in 1:length(epsilon)){
  # samples_df <- rej_abc(N = nthetas, epsilon = epsilon[i], y = y, rprior = rprior, simulate = simulate)
  args$epsilon <- epsilon[i]
  samples_df <- rej_abc(args)
  abc_df$samples[(1 + (i - 1) * nthetas): (i * nthetas)] <- samples_df$samples
}

plt_color <- scales::seq_gradient_pal(rgb(1, 0.5, 0.5), "darkblue")(seq(0, 1, length.out = length(epsilon)))
plt <- ggplot(abc_df) +
        # geom_line(data = posterior_df, aes(x = thetavals, y = abc_posterior), colour = "grey") +
        geom_density(aes(x = samples, colour = method)) +
        scale_color_manual(values = plt_color) +
        geom_line(data = posterior_df, aes(x = thetavals, y = true_posterior), colour = "black") +
        geom_vline(xintercept = theta_star$lambda, linetype = "dashed") +
        labs(x = "theta", y = "density") +
        theme(
              legend.position = c(.95, .95),
              legend.justification = c("right", "top"),
              legend.title=element_blank()
             ) +
        guides(color = guide_legend(override.aes = list(linetype = "solid")))

ggsave(plt, file = "plots/soft_abc/exponential_eg.pdf", height = 5)
