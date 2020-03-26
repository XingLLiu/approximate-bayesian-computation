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

l2norm <- function(z, y){
    return( sqrt(sum( abs( z - y )^2 )) )
}

l1norm <- function(y){
  return(sum(abs(y)))
}

discrepancy <- function(z, y){
    return( sqrt(sum( abs( z - y )^2 )) )
}

discrep_kernel <- function(rho, epsilon, kernel="uniform"){
    if (kernel == "uniform"){
        output <- uniform_kernel(rho / epsilon)
    } else if (kernel == "gaussian"){
        output <- gaussian_kernel(rho / epsilon) / epsilon
    } else if (kernel == "epan"){
        output <- epanechnikov_kernel(rho / epsilon) / epsilon
    }
    return(output)
}


# function to set up the index to store results in a data frame
index <- function(i, nthetas){
  return((1 + (i - 1) * nthetas): (i * nthetas))
}