source("src/kernels.R")

sumstat_fun <- function(z, FUN="original"){
    if (is.function(FUN)){
        return(FUN(z))
    } else if (FUN == "original"){
        return(z)
    } else if (FUN == "mean"){
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
