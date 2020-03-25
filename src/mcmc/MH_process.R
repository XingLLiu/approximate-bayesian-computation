mh_process <- function(samples, nsamples, burnin, thinning){
  # Apply burnin and thinning to a chain of MH samples
  if (class(samples) == "matrix"){
    return( samples[-(1:burnin), ][seq(from = 1, length.out = nsamples, by = thinning), ] )  
  }
  return( samples[-(1:burnin)][seq(from = 1, length.out = nsamples, by = thinning)] )
}

