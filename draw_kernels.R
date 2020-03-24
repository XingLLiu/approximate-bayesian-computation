# draw plots of 5 commonly used kernels

library("ggplot2")
source("src/kernels.R")

kernel_funcs <- c(uniform_kernel, triangular_kernel, epanechnikov_kernel, gaussian_kernel, laplacian_kernel)
xvals <- seq(-2.5, 2.5, by = 0.01)
yvals <- data.frame(x = rep(NA, length(kernel_funcs) * length(xvals)),
                    y = rep(NA, length(kernel_funcs) * length(xvals)),
                    kernel = rep(NA, length(kernel_funcs) * length(xvals)))
yvals$x <- rep(xvals, length(kernel_funcs))
yvals$kernel <- rep(c("Uniform", "Triangular", "Epanechnikov", "Gaussian", "Laplacian"),
                    each = length(xvals))


for (i in 1:length(kernel_funcs)){
  kernel <- kernel_funcs[[i]]
  yvals[(1 + (i - 1) * length(xvals)):(i * length(xvals)), "y"] <- kernel(xvals)
}

# create and save plot
# x = yvals$x, y = yvals$y, geom = "line", color = yvals$kernel, size = 1.2) +
plt <- ggplot(yvals, aes(x = x)) +
       geom_line(aes(y = y, color = kernel), size = 1.2) +
       xlab("u") +
       ylab("k(u)") +
       labs(color = "kernels") + 
       theme(
         legend.position = c(.95, .95),
         legend.justification = c("right", "top"),
         legend.text = element_text(size = 15),
         legend.title = element_text(size = 15),
         axis.text.x = element_text(size = 15),
         axis.text.y = element_text(size = 15),
         axis.title.x = element_text(size = 18),
         axis.title.y = element_text(size = 18),
       )

ggsave(plt, file = "plots/soft_abc/compare_kernels.pdf", height = 5)


