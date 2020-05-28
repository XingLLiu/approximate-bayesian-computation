# Approximate Bayesian Computation
Part III essay (Lent 2020). Supervised by Dr. Sergio Bacallado.

## Dependencies
```{r}
install.packages("Rcpp", "RcppEigen", "doParallel", "doRNG", "foreach", "ggplot2", "ggthemes", "dplyr", "reshape2", "BH", "transport", "Rmixmod", "gridExtra", "mclust")

# wabc
devtools::install_github("pierrejacob/winference")
```

## Folder Structure
    .
    ├── ...
    ├── archive             # archived files
    ├── plots               # generated plots
    ├── results             # saved results
    ├── src                 # source code for each algorithm
    └── README.md

## Abstract
The complexity of many data generating processes in real life either defies the access to the likelihood function or renders it too expansive to be evaluated. In this case, standard Bayesian inference techniques, such as Markov chain Monte Carlo, can no longer be used. A popular roundabout is Approximate Bayesian Computation (ABC). ABC only assumes one has a generative model from which data can be drawn. It relies on a user-specified discrepancy metric that compares some summaries of the observation and the generated data. However, an improperly selected metric or summary may bias the discrimination between models. Optimal transport (OT) metrics have recently been proposed to remedy this issue. OT metrics are flexible, admit decent convergence properties and are often able to capture all differences between distributions. In this essay, we review and compare two OT metrics and one information-based measure that arose in the ABC literature, namely the Wasserstein distances, the maximum mean discrepancy (MMD) and the Kullback-Leibler (KL) divergence. We summarize the theoretical studies of their posterior concentration in the present literature, and discuss how these metrics can be adapted to large-scale data sets. We also compare these methods through four benchmark experiments, including a real-life study on ecological dynamic systems.