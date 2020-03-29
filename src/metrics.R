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
    return( sqrt(sum( ( z - y )^2 )) )
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

# index <- function(method, dataframe){
#   return(dataframe %>% dplyr::filter(methods = method))
# }


# extract legend as a separate plot
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}


# initialize colours
init_colours <- function()
{
  return(c(Wasserstein = "darkblue", MMD = "#56B4E9", `KL divergence` = "darkgreen", Posterior = "grey50",
            Euclidean = "orange", `Euc. Summary` = "chocolate", Sinkhorn = "purple"))
}


# initialize discrete gradient colours
init_discrete_grad_colours <- function(n.out)
{
  return(scales::seq_gradient_pal(rgb(1, 0.5, 0.5), "darkblue")(seq(0, 1, length.out = n.out)))
}


# change curve and filling colours 
change_colour_and_fill <- function(colours)
{
  scale_color_manual(name = "", values = colours) +
  scale_fill_manual(name = "", values = colours)
}


# customize legend
add_legend <- function(xpos, ypos)
{
  theme(
        legend.position = c(xpos, ypos),
        legend.justification = c("right", "top"),
        legend.title=element_blank(),
        legend.key.size = unit(1.5, "line"),
        legend.text = element_text(size = 18)
        )
}


# customize label sizes
change_sizes <- function(test.size, title.size)
{
  theme(axis.text = element_text(size = test.size),
        axis.title = element_text(size = title.size, face = "bold"))
}

