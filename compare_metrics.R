library(winference)
library(rje)
source("src/metrics.R")
source("src/mmd/mmdsq_c.R")
# register parallel cores
registerDoParallel(cores = detectCores())
# apply preferences for ggplotting
require("gridExtra")
resultsprefix <- "results/compare_distances/"
plotprefix <- "plots/compare_distances/"
theme_set(theme_bw())

nobservation <- 1000
ydim <- c(1, 2, 5, 10)
theta0 <- 4
theta_vec <- seq(0.1, 10, length.out = 500)


compare_distances_fun <- function(i){
  set.seed(2020)
  y <- fast_rmvnorm(nobservation, mean = rep(0, ydim[i]), covariance = theta0 * diag(1, ydim[i]))

  dist <- list()
  dist$MMD <- dist$Wasserstein <- dist$kl.divergence <- rep(NA, length(theta_vec))

  # generate synthetic data
  obs_fake <- array(NA, c(nobservation, ydim[i], length(theta_vec)))
  for (j in 1:length(theta_vec)){
    obs_fake[, , j] <- fast_rmvnorm(nobservation, 
                                    mean = rep(0, ydim[i]),
                                    covariance = theta_vec[j] * diag(1, ydim[i]))
  }

  w1 <- rep(1/nobservation, nobservation)
  w2 <- rep(1/nobservation, nobservation)
  bandwidth <- median(apply(y, 1, l1norm))
  for (j in 1:length(theta_vec)){
    C <- cost_matrix_L2(t(y), t(matrix(obs_fake[, , j], ncol = ydim[i])))
    dist$MMD[j] <- mmdsq_c(y, matrix(obs_fake[, , j], ncol = ydim[i]), bandwidth)
    sink("/dev/null")
    dist$Wasserstein[j] <- as.numeric(exact_transport_given_C(w1, w2, C, p = 1))
    sink()
    dist$kl.divergence[j] <- FNN::KLx.divergence(y, matrix(obs_fake[, , j], ncol = ydim[i]), k = 1)[1]
    printPercentage(j, length(theta_vec))
  }

  # save results
  write.csv(dist, paste0(resultsprefix, "dim", ydim[i],"distances.csv"), row.names = FALSE)
  dist <- read.csv(paste0(resultsprefix, "dim", ydim[i],"distances.csv"))


  # integrate results into one dataframe
  method_names <- c("MMD", "Wasserstein", "KL Divergence")
  dist.df <- data.frame(methods = rep(method_names, each = length(theta_vec)),
                        thetas = rep(theta_vec, length(method_names)),
                        distances = c(dist$MMD, dist$Wasserstein, dist$kl.divergence)
                      )


  # plot results
  # my_colours <- init_colours()
  g1 <-  ggplot(data = dist.df, aes(x = thetas, y = distances, colour = methods)) +
            geom_point() +
            labs(x = "theta") +
            change_sizes(16, 20) +
            # scale_color_manual(name = "", values = my_colours) +
            geom_vline(xintercept = theta0, linetype = 2) +
            add_legend(0.95, 0.95)
  ggsave(g1, file = paste0(plotprefix, "dim", ydim[i], ".pdf"))
  # plot mmd alone
  g2 <- g1 + ylim() 
}


# run function
for (i in 1:length(ydim)){
  compare_distances_fun(i)
}



# check assumption 5
check_assumption_plots <- function(i){
  dist2 <- read.csv(paste0(resultsprefix, "dim", ydim[i],"distances.csv"))
  dist2$rho_theta <- sapply(theta_vec, FUN = function(a){ l2norm(a, theta0) })
  pdf(paste0(plotprefix, "check_assumption_dim", ydim[i], ".pdf"), width = 14)
  g1 <- ggplot(dist2, aes(x = Wasserstein, y = rho_theta)) +
          geom_point(colour = "darkblue") +
          labs(x = "Wasserstein distances", y = "|theta - theta_star|") +
          change_sizes(16, 20) +
          theme(legend.position = "none")
  g2 <- ggplot(dist2, aes(x = MMD, y = rho_theta)) +
          geom_point(colour = "orange") +
          labs(x = "MMD^2", y = "|theta - theta_star|") +
          change_sizes(16, 20) +
          theme(legend.position = "none")
  gridExtra::grid.arrange(g1, g2, ncol = 2)
  dev.off()
  # log scale
  pdf(paste0(plotprefix, "check_assumption_dim", ydim[i], "_log.pdf"), width = 14)
  g1 <- ggplot(dist2, aes(x = log(Wasserstein), y = log(rho_theta))) +
          geom_point(colour = "darkblue") +
          labs(x = "log Wasserstein distances", y = "log |theta - theta_star|") +
          change_sizes(16, 20) +
          theme(legend.position = "none")
  g2 <- ggplot(dist2, aes(x = log(MMD), y = log(rho_theta))) +
          geom_point(colour = "orange") +
          labs(x = "log MMD^2", y = "log |theta - theta_star|") +
          change_sizes(16, 20) +
          theme(legend.position = "none")
  gridExtra::grid.arrange(g1, g2, ncol = 2)
  dev.off()
}

# run function2
for (i in (1:length(ydim))){
  check_assumption_plots(i)
}



