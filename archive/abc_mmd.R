library(kernlab)
source("src/mmd/mmd.R")

gamma <- 2
rbf <- rbfdot(1/gamma)

nvec <- c(10, 100, 200, 300, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000)
mmd_sq_vec <- rep(NA, length(nvec))
mmd_sq_rf_vec <- rep(NA, length(nvec))
sigma1 <- 0.5
sigma2 <- 0.25
rbf <- rbfdot(0.5)
set.seed(2018)
# Assuming we use a rbf kernel of sigma = 0.5 following the notation of rbfdot, and
# Y ~ N(0, \sigma_1^2), Z ~ N(0, \sigma_2^2). Then we have the following MMD:
mmd_sq_true <- 1/sqrt(1 + 2 * sigma1^2) + 1/sqrt(1 + 2 * sigma2^2) - 2/sqrt(1 + sigma1^2 + sigma2^2)
for (i in 1:length(nvec)){
  y <- rnorm(nvec[i], sd = sigma1)
  z <- rnorm(nvec[i], sd = sigma2)
  mmd_sq_vec[i] <- mmd_sq(rbf, y, z)
  # mmd_sq_rf_vec[i] <- mmd_sq(rbf, y, z, method = "rf", size = 1e5)
  print(nvec[i])
  # print(c(mmd_sq_true, mmd_sq_vec[i], mmd_sq_rf_vec[i]))
}

plot(nvec, mmd_sq_vec)
points(nvec, mmd_sq_rf_vec, col = "blue")
abline(h = mmd_sq_true)



# plot \rho(\theta, \theta_*) against MMD^2(\mu_\theta, \mu_\ast)
set.seed(2018)
gamma <- 2
rbf <- rbfdot(1/gamma)
nsample <- 5000
sigma1 <- 0.5
sigma2 <- seq(0.1, 0.9, by = 0.05)
mmd_sq_vec <- rep(NA, length(sigma2))
rho_vec <- rep(NA, length(sigma2))
rbf <- rbfdot(0.5)

for (i in 1:length(sigma2)){
  y <- rnorm(nsample, sd = sigma1)
  z <- rnorm(nsample, sd = sigma2[i])
  mmd_sq_vec[i] <- mmd_sq(rbf, y, z)
  rho_vec[i] <- (sigma1 - sigma2[i])^2
  print(sprintf("[%i / %i]", i, length(sigma2)))
}

plot(mmd_sq_vec, rho_vec)
