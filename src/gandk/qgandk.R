# include <Rcpp.h>
Rcpp::cppFunction("
  NumericVector qgandk(NumericVector x, NumericVector theta){
    NumericVector z = qnorm(x);
    double A = theta(0);
    double B = theta(1);
    double c = 0.8;
    double g = theta(2);
    double k = theta(3);
    return A + B * (1 + c * (1 - exp(- g * z)) / (1 + exp(- g * z))) * pow((1 + pow(z, 2.0)), k) * z;
  }
")