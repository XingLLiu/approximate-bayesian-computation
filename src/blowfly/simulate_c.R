# Function to simulate from the blowfly model
Rcpp::cppFunction("
  NumericMatrix simulate_c(int n, int burnin, int lag, NumericVector eps_t, NumericVector e_t, 
                            double P, double delta,
                            double N0, double sigmad, double sigmap, int tau, NumericMatrix dat){
    //NumericVector eps_t = rgamma(burnin + n, 1 / std::pow(sigmad, 2), std::pow(sigmad, 2));
    //NumericVector e_t = rgamma(burnin + n, 1 / std::pow(sigmap, 2), std::pow(sigmap, 2));
    int t;
    int tau_t;

    for (int i = 0; i < (burnin + n); i ++){
      t = i + lag;
      tau_t = t - lag;

      dat(t, 0) = P * dat(tau_t, 0) * exp(- dat(tau_t, 0)/N0) * e_t(i) + dat(t-1, 0) * exp(- delta*eps_t(i));
    }

    return dat;
  } 
")
