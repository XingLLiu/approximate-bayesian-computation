randFourierRBF <- function(size, y, z, rbf)
# Compute the kernel matrix of a RBF kernel with random Fourer feature.
# For RBF \exp( - \| x - y \| / (2\sigma^2) ), W ~ normal(0, \sigma^2).
# Input :
#         size: Length of random Fourier map to approximate the kernel.
#         y, z: data of dimension 1.
#         rbf: A rbf object instantiated by 'rbfdot'.
# Output:
#         A (size, length(y)) matrix, with (l, i)th entry being \sqrt(2/size) * cos(w_l^T y^(i) + u_l).
{
  # if (is.null(dim(y))){
  #   y <- matrix(y, ncol = 1)
  #   z <- matrix(z, ncol = 1)
  # }
  # dim <- dim(y)[2]
  n <- length(y)

  w <- matrix(rnorm(size * n, sd = sqrt(2 * rbf@kpar$sigma)), ncol = n)
  u <- matrix(runif(size * n, -pi, pi), ncol = n)

  return( list(
                # phi.y = sqrt(2 / size) * cos(w %*% t(y) + matrix(rep(u, dim(y)[1]), ncol = dim(y)[1])),
                # phi.z = sqrt(2 / size) * cos(w %*% t(z) + matrix(rep(u, dim(z)[1]), ncol = dim(z)[1]))
                # phi.y = sqrt(2 / size) * cos(w %*% t(y) + u),
                # phi.z = sqrt(2 / size) * cos(w %*% t(z) + u)

                phi.y = sqrt(2 / size) * cos(w %*% diag(y) + u),
                phi.z = sqrt(2 / size) * cos(w %*% diag(z) + u)
              )
        )
}

rf <- randFourierRBF(1000, y, z, rbf)
c(t(rf$phi.y[, 1]) %*% rf$phi.z[, 1], rbf(y[1], z[1]))
c(t(rf$phi.y[, 2]) %*% rf$phi.z[, 2], rbf(y[2], z[2]))

