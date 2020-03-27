# Function to compute the MMD^2 by the unbiased estimator
Rcpp::cppFunction("
  double mmdsq_c(NumericMatrix y, NumericMatrix z, double epsilon){
    int nobs = y.rows();
    int dimension = y.cols();
    double result = 0;
    double cost_yy = 0;
    double cost_zz = 0;
    double cost_yz = 0;
    double zz = 0;
    double yy = 0;
    double yz = 0;

    for (int i1 = 0; i1 < nobs; i1 ++){
      for (int i2 = 0; i2 < nobs; i2 ++){
        cost_yy = 0;
        cost_zz = 0;
        cost_yz = 0;

        for (int j = 0; j < dimension; j ++){
          if (i1 != i2){
            cost_yy += std::pow(y(i1, j) - y(i2, j), 2);
            cost_zz += std::pow(z(i1, j) - z(i2, j), 2);
          }
          
          yy += exp(- cost_yy / (2 * epsilon * epsilon));
          zz += exp(- cost_zz / (2 * epsilon * epsilon));
          cost_yz += std::pow(y(i1, j) - z(i2, j), 2);
          yz += exp(- cost_yz / (2 * epsilon * epsilon));
        }
      }
    }
    result = yy / (nobs * (nobs - 1)) + zz / (nobs * (nobs - 1)) - 2 * yz / (nobs * nobs);
    return result;
  }
")


