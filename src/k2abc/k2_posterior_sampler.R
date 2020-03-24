draw_k2_sample <- function(size, k2abc)
# draw posterior samples from the K2ABC weighted posterior particles.
{
  return( sample(size = size, x = k2abc$prior_samples, prob = k2abc$weights, replace = TRUE) )
}
