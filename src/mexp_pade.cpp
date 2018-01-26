#include <Rcpp.h>

#include "gaussinte.hpp"

using namespace Rcpp;
using namespace mylib;

// ff1 <- function(t0, t1, lam, a, mu, b) {
// #	cat("ff1", t0, t1, lam, a, mu, b, "\n")
// 	f <- function(x, tk, lam, a, mu, b) {
// 		dgamma(x=x, shape=a, rate=lam) * pgamma(q=tk-x, shape=b, rate=mu)
// 	}
// 	res0 <- integrate(f=f, lower=t0, upper=t1, tk=t1, a=a, b=b, lam=lam, mu=mu, abs.tol=0)$value
// #	res <- deformula.moneone(f=f, lower=t0, upper=t1, tk=t1, a=a, b=b, lam=lam, mu=mu)$value
// #	cat(res0, "\n")
// 	res0
// }

double ff1integrand(double x, double tk, double lam, int a, double mu, int b) {
  return R::dgamma(x, a, 1/lam, false) * R::pgamma(tk-x, b, 1/mu, true, false);
}

double ff1(double t0, double t1, double lam, int a, double mu, int b, const gaussinte<NumericVector>& g, NumericVector fx) {
  double c = g.comp_fx(t0, t1, fx);
  for (int i=0; i<g.get_n(); i++) {
    fx[i] = ff1integrand(fx[i], t1, lam, a, mu, b);
  }
  return g.comp_value(c, fx);
}

// [[Rcpp::export]]

double test(int n, double t0, double t1, double lam, int a, double mu, int b) {
  NumericVector x(n);
  NumericVector w(n);
  NumericVector fx(n);
  gaussinte<NumericVector> g(n, x, w, 1.0e-8);

  return ff1(t0, t1, lam, a, mu, b, g, fx);
}
