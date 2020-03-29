library(winference)

sabc <- function(args, ...)
#'@description  Implement the SMC-ABC with a given data discrepancy metric.
{
  param_algo <- list(nthetas = args$nthetas, nmoves = 1, proposal = mixture_rmixmod(),
                   minimum_diversity = 0.5, R = 2, maxtrials = 1e5)
  smcresults <- wsmc(args$discrepancy, args, param_algo, ...)  
}


sabc_get_last_samples <- function(sabc.results)
#'@description  Extract the samples from the last step in 'sabc.result'.
{
  samples.df <- wsmc_to_dataframe(sabc.results)
  samples.df <- filter(samples.df, step == samples.df$step[nrow(samples.df)])
}

