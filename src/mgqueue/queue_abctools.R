library(abctools)
library(winference)
rm(list = ls())
registerDoParallel(cores = detectCores())
setmytheme()
set.seed(2018)
target <- get_queue()

prefix = "results/queue/"

nobservations <- 50
load(paste0(prefix, "50.intermediateobs.neal.RData"))
obs = matrix(obs, nrow = 1)


target$simulate <- function(theta){
  return(matrix(target$robservation(nobservations, theta, target$parameters), nrow = 1))
}


load(paste0(prefix, "queue_intermediate_wsmc_marginal.RData"))
results_was = results
#ndists = cumsum(results$ncomputed)
#nsim = ndists[34]
# nsim_approx = 10^6
# nsim = ndists[which.min(abs(nsim_approx - ndists))]
nsim = 10^7
m = 20 #determines number of order stats in initial summary stat
l = 1 #determines how many powers of the order stats are included in the initial summary stat

obs_subset = sort(obs)[c(1,(1:nobservations)[(nobservations/m)*(1:(m-2))+0.5],nobservations)]
