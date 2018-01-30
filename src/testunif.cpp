#include <Rcpp.h>

#include "mexp_unif.hpp"

using namespace Rcpp;
using namespace marlib;

// [[Rcpp::export]]

NumericMatrix test_mexp_unif(char trans, NumericMatrix P, double qv, double t, NumericMatrix x) {
  int right = poisson<NumericVector>::rightbound(qv*t, 1.0e-8);
  NumericVector prob(right+1);
  poisson<NumericVector> pois(prob.size(), prob);
  pois.set(qv*t, 0, right);
  mexp_unif_notrans(P, qv, pois, x, x, 0);
  return x;
}
