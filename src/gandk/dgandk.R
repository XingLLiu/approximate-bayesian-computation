source("src/gandk/pgandk.R")

logdgandk <- function(ys, thetas, ...){
  n <- length(ys)
  evals <- rep(0, nrow(thetas))
  for (itheta in 1:nrow(thetas)){
    ll <- function(ys, h = 1e-5, tolerance = 1e-10){
      all_ys <- c(ys-h, ys+h) # for finite difference differentiation
      o <- order(all_ys)
      x <- rep(0, length(all_ys))
      x[o[1]] <- pgandk(y = all_ys[o[1]], theta = thetas[itheta,], tolerance = tolerance)
      for (i in 2:length(all_ys)){
        x[o[i]] <- pgandk(y = all_ys[o[i]], theta = thetas[itheta,], tolerance = tolerance, lower = x[o[i-1]])
      }
      return(sum(log((x[(n+1):(2*n)] - x[1:n])/(2*h))))
    }
    evals[itheta] <- ll(ys)
  }
  return(evals)
}

dgandk <- function(ys, thetas, ...){
  return(exp(logdgandk(ys, thetas)))
}