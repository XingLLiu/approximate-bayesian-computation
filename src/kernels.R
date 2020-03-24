uniform_kernel <- function(u){
  return(0.5 * (abs(u) <= 1))
}


triangular_kernel <- function(u){
  return( (1 - abs(u)) * (abs(u) <= 1) )
}


epanechnikov_kernel <- function(u){
  return( 0.75 * (1 - u^2) * (abs(u) <= 1) )
}


gaussian_kernel <- function(u){
  return( 1 / sqrt(2 * pi) * exp( - 0.5 * u^2) )
}


laplacian_kernel <- function(u){
  return( 0.4 * exp(- 0.5 * abs(u)) )
}

