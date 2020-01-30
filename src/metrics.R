summary_stat <- function(z, family="original"){
    if (family == "original"){
        return(z)
    } else if (family == "mean"){
        return(mean(z))
    }
}

discrepancy <- function(z, y){
    return( mean( abs( z - y ) ) )
}

discrep_kernel <- function(z, y, epsilon, kernel="identity"){
    # compute discrepancy
    rho <- discrepancy(z, y)
    if (kernel == "identity"){
        output <- rho <= epsilon
    } else if (kernel == "rbf"){
        output <- exp(- rho^2 / (2 * epsilon^2)) / (sqrt(2 * pi) * epsilon)
    } else if (kernel == "epan"){
        output <- 0.75 * (1 - rho^2 / epsilon^2) * (rho <= epsilon)
    }
    return(output)
}
