# Exponential with Gamma prior
library(ggplot2)
set.seed(2020)
theta_star <- list(lambda = 2)
hyperparams <- list(alpha = 1, beta = 1)
epsilon <- c(0.05, 0.1, 0.5, 1)
nobservation <- 1
nthetas <- 1024
y <- rexp(nobservation, theta_star$lambda)


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

# # rejection with identity kernel
# source("src/abc_rej.R")
# samples_unif <- soft_abc(N = 1000, epsilon = epsilon, y = y, alpha = hyperparams$alpha, beta = hyperparams$beta)
# # rejection with Gaussian kernel
# samples_gaus <- soft_abc(N = 1000, epsilon = epsilon, y = y, alpha = hyperparams$alpha, beta = hyperparams$beta,
#                           kernel = "gaussian")


# initialize dataframes to store samples
method_names <- paste("epsilon =", epsilon)
abc_df <- data.frame(method = rep(method_names, each = nthetas),
                     samples = NA
                    )

# rejection abc
source("src/abc_rej.R")
for (i in 1:length(epsilon)){
  samples_df <- soft_abc(N = nthetas, epsilon = epsilon[i], y = y, alpha = hyperparams$alpha, beta = hyperparams$beta)
  abc_df$samples[(1 + (i - 1) * nthetas): (i * nthetas)] <- samples_df$samples
}

plt <- ggplot(abc_df) +
        # geom_line(data = posterior_df, aes(x = thetavals, y = abc_posterior), colour = "grey") +
        geom_density(aes(x = samples, colour = method)) +
        scale_color_manual(values = scales::seq_gradient_pal(rgb(1, 0.5, 0.5), "darkblue")(seq(0, 1, length.out = length(epsilon)))) +
        geom_line(data = posterior_df, aes(x = thetavals, y = true_posterior), colour = "black") +
        labs(x = "theta", y = "density") +
        theme(
              legend.position = c(.95, .95),
              legend.justification = c("right", "top"),
              legend.title=element_blank()
             ) +
        guides(color = guide_legend(override.aes = list(linetype = "solid")))

ggsave(plt, file = "plots/soft_abc/exponential_eg.pdf", height = 5)



# # integrate data into dataframes
# abc_df <- data.frame(thetavals = thetavals,
#                      true_posterior = true_posterior,
#                      abcposterior = abcposterior,
#                      unif = samples_unif$samples,
#                      gaus = samples_gaus$samples
#                     )
# samples_df <- data.frame(unif = samples_unif$samples, gaus = samples_gaus$samples)

# # plot results
# plt  <- ggplot(abc_df) +
#         geom_line(aes(x = thetavals, y = true_posterior, colour = "true posterior"), size = 1.1) +
#         geom_line(aes(x = thetavals, y = abcposterior, colour = "exact abc posterior (uniform)"), size = 1.1) + 
#         geom_density(data = abc_df, aes(x = unif, colour = "estimated abc posterior (uniform)"), size = 1.1) +
#         geom_density(data = abc_df, aes(x = gaus, colour = "estimated abc posterior (Gaussian)"), size = 1.1) +
#         labs(x = "theta", y = "density", colour = "Legend") +
#         scale_color_manual(values = c("red", "blue", "grey", "black")) +
#         theme(
#           legend.position = c(.95, .95),
#           legend.justification = c("right", "top"),
#           legend.text = element_text(size = 15),
#           legend.title = element_text(size = 15),
#           axis.text.x = element_text(size = 15),
#           axis.text.y = element_text(size = 15),
#           axis.title.x = element_text(size = 18),
#           axis.title.y = element_text(size = 18)
#         )


# ggsave(plt, file = "plots/soft_abc/exponential_eg.pdf", height = 5)

