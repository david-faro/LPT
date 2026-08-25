#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame cpp_update_position_inflow(double x, double y, double ux, double uy, double up, double K, double dt, std::string drift_destination) {
  NumericVector L = Rcpp::rnorm(2, 0, 1);
  
  double upx = 0;
  double upy = 0;
  
  double xn, yn;
  
  if (drift_destination == "in flow") {
    xn = x + ux * dt + upx * dt + L[0] * std::sqrt(2 * K * dt);
    yn = y + uy * dt + upy * dt + L[1] * std::sqrt(2 * K * dt);
  } else if (drift_destination == "outside boundary") {
    xn = x + upx * dt + L[0] * std::sqrt(2 * K * dt);
    yn = y + upy * dt + L[1] * std::sqrt(2 * K * dt);
  }
  
  return DataFrame::create(Named("x") = xn,
                           Named("y") = yn);
}
