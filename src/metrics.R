source("src/kernels.R")

summary_stat <- function(z, family="original"){
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
        # output <- rho <= epsilon
        output <- uniform_kernel(rho / epsilon) / epsilon
    } else if (kernel == "gaussian"){
        # output <- exp(- rho^2 / (2 * epsilon^2)) / (sqrt(2 * pi) * epsilon)
        output <- gaussian_kernel(rho / epsilon) / epsilon
    } else if (kernel == "epan"){
        # output <- 0.75 * (1 - rho^2 / epsilon^2) * (rho <= epsilon)
        output <- epanechnikov_kernel(rho / epsilon) / epsilon
    }
    return(output)
}
