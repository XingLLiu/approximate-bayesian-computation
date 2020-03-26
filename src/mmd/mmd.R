library(kernlab)
source("src/mmd/off_diag_sum.R")
source("src/mmd/rand_fourier_maps.R")

mmd_sq <- function(kernel, y, z, method="ordinary", size=1000)
# kernel: A S4 object of class 'kernel' as in 'kernlab'.
# size: Length of the random Fourier feature vector. Must be specified if method == rf.
{

  ny <- length(y)
  nz <- length(z)

  if (method == "ordinary"){
    # compute kernels
    kernelmat_y <- kernelMatrix(kernel, y)
    yy <- off_diag_sum(kernelmat_y) / (ny * (ny - 1))
    rm(kernelmat_y)

    kernelmat_z <- kernelMatrix(kernel, z)
    zz <- off_diag_sum(kernelmat_z) / (nz * (nz - 1))
    rm(kernelmat_z)

    kernelmat_yz <- kernelMatrix(kernel, y, z)
    yz <- sum(kernelmat_yz) / (ny * nz)
    rm(kernelmat_yz)

    mmd_sq <- yy + zz - 2 * yz

  } else if (method == "rf") {
    if (class(rbf) != "rbfkernel") {
      stop("kernel has to be of class 'rbfkernel' to apply the RBF random Fourier feature!")
    }
    rf <- randFourierRBF(size, y, z, kernel)
    mmd_sq <- sum(rowMeans(rf$phi.y - rf$phi.z)^2)
    rm(rf)
  }

  return(mmd_sq)
}
