library("mvtnorm")
library("FNN")
library("winference")
# register parallel cores
registerDoParallel(cores = detectCores())
# apply preferences for ggplotting
require("gridExtra")
theme_set(theme_bw())

set.seed(1)

nsample <- 1000
dim <- 2
theta0 <- 4
theta_vec <- seq(0.1, 9, length.out = 100)
obs <- rmvnorm(nsample, mean = rep(0, dim), sigma = theta0 * diag(1, dim))

dist <- list()
dist$wasserstein <- dist$sinkhorn <- dist$kl <- rep(NA, length(theta_vec))

obs_fake <- array(NA, c(nsample, dim, length(theta_vec)))

for (j in 1:length(theta_vec)){
  obs_fake[, , j] <- rmvnorm(nsample, mean = rep(0, dim), sigma = theta_vec[j] * diag(1, dim))
}

w1 <- rep(1/nsample, nsample)
w2 <- rep(1/nsample, nsample)
for (j in 1:length(theta_vec)){
  C <- cost_matrix_L2(t(obs), t(obs_fake[, , j]))

  dist$wasserstein[j] <- as.numeric(exact_transport_given_C(w1, w2, C, p = 1))
  dist$sinkhorn[j] <- sinkhorn_given_C(w1, w2, C, p = 1, eps = 0.05, niterations = 100)$corrected  
  dist$kl[j] <- KLx.divergence(obs, obs_fake[, , j])[1]
  if (j %% 10 == 0){
    print(paste("iter:", j))
  }
}

plot(theta_vec, dist$wasserstein, ylim = c(0, range(dist$wasserstein)[2]), col = "red")
points(theta_vec, dist$sinkhorn_distance, col = "green")
points(theta_vec, dist$kl * (dist$kl > 0), col = "blue")
legend("topright", legend = c("wasserstein", "sinkhorn", "KL"), lty = c(1, 1, 1),
        col = c("red", "green", "blue"))
abline(v = theta0)

# qplot(x = theta_vec, y = dist$wasserstein, geom = "blank") +
#      geom_point(aes(colour = "wasserstein")) +
#      geom_vline(xintercept = theta0) +
#      geom_point(aes(x = theta_vec, y = dist$sinkhorn_distance, colour = "sinkhorn")) +
#      geom_point(aes(x = theta_vec, y = dist$kl, colour = "KL")) +
#      xlab("theta") + 
#      ylab("distance") + 
#      ylim(c(0, range(dist$wasserstein)[2])) +
#      scale_colour_manual(name = "", values = c("red", "blue", "green"))


