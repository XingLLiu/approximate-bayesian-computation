off_diag_sum <- function(X)
# compute the off-diagonal sum of matrix X
{
  if (!is.matrix(X)){
    stop(sprintf("Class must be matrix but got %s instead!", class(X)))
  }  

  return(sum(X) - sum(diag(X)))
}
