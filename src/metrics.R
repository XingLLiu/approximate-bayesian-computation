source("src/kernels.R")

sumstat <- function(z, family="original"){
    if (family == "original"){
        return(z)
    } else if (family == "mean"){
        return(mean(z))
    }
}

discrepancy <- function(z, y){
    return( sqrt(sum( abs( z - y )^2 )) )
}

discrep_kernel <- function(z, y, epsilon, kernel="uniform"){
    # compute discrepancy
    rho <- discrepancy(z, y)
    if (kernel == "uniform"){
        output <- uniform_kernel(rho / epsilon) / epsilon
    } else if (kernel == "gaussian"){
        output <- gaussian_kernel(rho / epsilon) / epsilon
    } else if (kernel == "epan"){
        output <- epanechnikov_kernel(rho / epsilon) / epsilon
    }
    return(output)
}
