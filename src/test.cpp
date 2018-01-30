#include <Rcpp.h>

#include "cppblas.hpp"

using namespace Rcpp;
using namespace marlib;

// [[Rcpp::export]]

NumericVector test_daxpy(double a, NumericVector x, NumericVector y) {
  return marlib::daxpy(a, x, y);
}

// [[Rcpp::export]]

NumericVector test_dscal(double a, NumericVector x) {
  return marlib::dscal(a, x);
}

// [[Rcpp::export]]

double test_dasum(NumericMatrix x) {
  return marlib::dasum(x);
}

// [[Rcpp::export]]

NumericVector test_dgemv(char trans, double a, NumericMatrix A, NumericVector x, double b, NumericVector y) {
  switch (trans) {
    case 'N':
    case 'n':
      return marlib::dgemvN(a, A, x, b, y);
      break;
    case 'T':
    case 't':
      return marlib::dgemvT(a, A, x, b, y);
      break;
    default:
      return NULL;
  }
}

// [[Rcpp::export]]

NumericMatrix test_dgemm(char transA, char transB, double a, NumericMatrix A, NumericMatrix B, double b, NumericMatrix C) {
  switch (transA) {
    case 'N':
    case 'n':
      switch (transB) {
        case 'N':
        case 'n':
          return marlib::dgemmNN(a, A, B, b, C);
          break;
        case 'T':
        case 't':
          return marlib::dgemmNT(a, A, B, b, C);
          break;
      }
      break;
    case 'T':
    case 't':
      switch (transB) {
        case 'N':
        case 'n':
          return marlib::dgemmTN(a, A, B, b, C);
          break;
        case 'T':
        case 't':
          return marlib::dgemmTT(a, A, B, b, C);
          break;
      }
      break;
  }
}
