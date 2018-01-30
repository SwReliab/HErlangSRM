/*
  blas.hpp
*/

#include "dblas.h"

namespace marlib {

  // dfill

  template <typename VectorT>
  inline
  VectorT& dfill(VectorT& x, double alpha) {
    for (auto xptr = x.begin(); xptr != x.end(); xptr++) {
      *xptr = alpha;
    }
    return x;
  }

  // dcopy

  template <typename VectorT>
  inline
  VectorT& dcopy(const VectorT& x, VectorT& y) {
    dblas::dcopy(x.size(), &x[0], 1, &y[0], 1);
    return y;
  }

  // BLAS level 1

  // daxpy

  template <typename VectorT>
  inline
  VectorT& daxpy(double alpha, const VectorT& x, VectorT& y) {
    dblas::daxpy(x.size(), alpha, &x[0], 1, &y[0], 1);
    return y;
  }

  // dscal

  template <typename VectorT>
  inline
  VectorT& dscal(double alpha, VectorT& x) {
    dblas::dscal(x.size(), alpha, &x[0], 1);
    return x;
  }

  // dasum

  template <typename VectorT>
  inline
  double dasum(const VectorT& x) {
    return dblas::dasum(x.size(), &x[0], 1);
  }

  // BLAS level 2

  // dgemv

  template <typename VectorT, typename MatrixT>
  inline
  VectorT& dgemvN(double alpha, const MatrixT& A, const VectorT& x, double beta, VectorT& y) {
    int m = y.size();
    int n = x.size();
    dblas::dgemv('N', m, n, alpha, &A[0], m, &x[0], 1, beta, &y[0], 1);
    return y;
  }

  template <typename VectorT, typename MatrixT>
  inline
  VectorT& dgemvT(double alpha, const MatrixT& A, const VectorT& x, double beta, VectorT& y) {
    int m = x.size();
    int n = y.size();
    dblas::dgemv('T', m, n, alpha, &A[0], m, &x[0], 1, beta, &y[0], 1);
    return y;
  }

//   template <typename ValueT>
//   inline
//   vector<ValueT>& dgemv(const trans_t& trans,
//     const ValueT alpha, const csr_matrix<ValueT>& A,
//     const vector<ValueT>& x, const ValueT beta, vector<ValueT>& y) {
//
//     switch (trans) {
//       case NoTrans:
//       y *= beta;
//       for (size_type i=0; i<A.nrow(); i++) {
//         for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
//           size_type j = A.colind(z);
//           y(i + y.begin()) += alpha * A.value(z) * x(j + x.begin());
//         }
//       }
//       break;
//       case Trans:
//       y *= beta;
//       for (size_type i=0; i<A.nrow(); i++) {
//         for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
//           size_type j = A.colind(z);
//           y(j + y.begin()) += alpha * A.value(z) * x(i + x.begin());
//         }
//       }
//       break;
//     }
//     return y;
//   }

//   // dger
//
//   template <typename ValueT>
//   inline
//   dense_matrix<ValueT>& dger(const trans_t& trans,
//     const ValueT alpha, const vector<ValueT>& x, const vector<ValueT>& y,
//     dense_matrix<ValueT>& A) {
//
//     switch (trans) {
//       case NoTrans:
//       assert(A.nrow() == x.size());
//       assert(A.ncol() == y.size());
//       for (index_type ja=A.cbegin(), j=y.begin(); j<=y.end(); ja++, j++) {
//         for (index_type ia=A.rbegin(), i=x.begin(); i<=x.end(); ia++, i++) {
//           A(ia,ja) += alpha * x(i) * y(j);
//         }
//       }
//       break;
//       case Trans:
//       assert(A.ncol() == x.size());
//       assert(A.nrow() == y.size());
//       for (index_type ja=A.cbegin(), j=x.begin(); j<=x.end(); ja++, j++) {
//         for (index_type ia=A.rbegin(), i=y.begin(); i<=y.end(); ia++, i++) {
//           A(ia,ja) += alpha * x(j) * y(i);
//         }
//       }
//       break;
//     }
//     return A;
//   }
//
//   template <typename ValueT>
//   inline
//   csr_matrix<ValueT>& dger(const trans_t& trans,
//     const ValueT alpha, const vector<ValueT>& x, const vector<ValueT>& y,
//     csr_matrix<ValueT>& A) {
//
//     switch (trans) {
//       case NoTrans:
//       for (size_type i=0; i<A.nrow(); i++) {
//         for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
//           size_type j = A.colind(z);
//           ValueT tmp = alpha;
//           tmp *= x(i + x.begin());
//           tmp *= y(j + y.begin());
//           A.value(z) += tmp;
//         }
//       }
//       break;
//       case Trans:
//       for (size_type i=0; i<A.nrow(); i++) {
//         for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
//           size_type j = A.colind(z);
//           ValueT tmp = alpha;
//           tmp *= x(j + x.begin());
//           tmp *= y(i + y.begin());
//           A.value(z) += tmp;
//         }
//       }
//       break;
//     }
//     return A;
//   }
//
// #ifdef F77BLAS
//   template <>
//   inline
//   dense_matrix<double>& dger(
//     const trans_t& trans,
//     const double alpha, const vector<double>& x, const vector<double>& y,
//     dense_matrix<double>& A) {
//     switch (trans) {
//       case NoTrans:
//         dblas::dger(A.nrow(), A.ncol(), alpha, x.ptr(), x.inc(), y.ptr(), y.inc(), A.ptr(), A.ld());
//         break;
//       case Trans:
//         dblas::dger(A.nrow(), A.ncol(), alpha, y.ptr(), y.inc(), x.ptr(), x.inc(), A.ptr(), A.ld());
//         break;
//     }
//     return A;
//   }
// #endif
//

  // BLAS level 3

  // dgemm

  template <typename MatrixT, typename MatrixT2>
  inline
  MatrixT& dgemmNN(double alpha, const MatrixT2& A, const MatrixT& B, double beta, MatrixT& C) {
    dblas::dgemm('N', 'N', C.nrow(), C.ncol(), A.ncol(), alpha, &A[0], A.nrow(), &B[0], B.nrow(),
      beta, &C[0], C.nrow());
    return C;
  }

  template <typename MatrixT, typename MatrixT2>
  inline
  MatrixT& dgemmNT(double alpha, const MatrixT2& A, const MatrixT& B, double beta, MatrixT& C) {
    dblas::dgemm('N', 'T', C.nrow(), C.ncol(), A.ncol(), alpha, &A[0], A.nrow(), &B[0], B.nrow(),
      beta, &C[0], C.nrow());
    return C;
  }

  template <typename MatrixT, typename MatrixT2>
  inline
  MatrixT& dgemmTN(double alpha, const MatrixT2& A, const MatrixT& B, double beta, MatrixT& C) {
    dblas::dgemm('T', 'N', C.nrow(), C.ncol(), A.nrow(), alpha, &A[0], A.nrow(), &B[0], B.nrow(),
      beta, &C[0], C.nrow());
    return C;
  }

  template <typename MatrixT, typename MatrixT2>
  inline
  MatrixT& dgemmTT(double alpha, const MatrixT2& A, const MatrixT& B, double beta, MatrixT& C) {
    dblas::dgemm('T', 'T', C.nrow(), C.ncol(), A.nrow(), alpha, &A[0], A.nrow(), &B[0], B.nrow(),
      beta, &C[0], C.nrow());
    return C;
  }

  ///

  template <typename MatrixT>
  inline
  Rcpp::NumericVector& dgemmNN(double alpha, const MatrixT& A, const Rcpp::NumericVector& B,
    double beta, Rcpp::NumericVector& C) {
    return dgemvN(alpha, A, B, beta, C);
  }

  template <typename MatrixT>
  inline
  Rcpp::NumericVector& dgemmNT(double alpha, const MatrixT& A, const Rcpp::NumericVector& B,
    double beta, Rcpp::NumericVector& C) {
    return dgemvN(alpha, A, B, beta, C);
  }

  template <typename MatrixT>
  inline
  Rcpp::NumericVector& dgemmTN(double alpha, const MatrixT& A, const Rcpp::NumericVector& B,
    double beta, Rcpp::NumericVector& C) {
    return dgemvT(alpha, A, B, beta, C);
  }

  template <typename MatrixT>
  inline
  Rcpp::NumericVector& dgemmTT(double alpha, const MatrixT& A, const Rcpp::NumericVector& B,
    double beta, Rcpp::NumericVector& C) {
    return dgemvT(alpha, A, B, beta, C);
  }

//   template <typename ValueT>
//   inline
//   dense_matrix<ValueT>& dgemm(
//     const trans_t& transA,
//     const trans_t& transB,
//     const ValueT alpha,
//     const csr_matrix<ValueT>& A,
//     const dense_matrix<ValueT>& B,
//     const ValueT beta,
//     dense_matrix<ValueT>& C) {
//
//     switch (transB) {
//       case NoTrans:
//       switch (transA) {
//         case NoTrans:
//         C *= beta;
//         for (size_type i=0; i<A.nrow(); i++) {
//           for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
//             size_type j = A.colind(z);
//             for (index_type jx=B.cbegin(), jy=C.cbegin(); jx<=B.cend(); jx++, jy++) {
//               C(i + C.rbegin(), jy) += alpha * A.value(z) * B(j + B.rbegin(), jx);
//             }
//           }
//         }
//         break;
//         case Trans:
//         C *= beta;
//         for (size_type i=0; i<A.nrow(); i++) {
//           for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
//             size_type j = A.colind(z);
//             for (index_type jx=B.cbegin(), jy=C.cbegin(); jx<=B.cend(); jx++, jy++) {
//               C(j + C.rbegin(), jy) += alpha * A.value(z) * B(i + B.rbegin(), jx);
//             }
//           }
//         }
//         break;
//       }
//       break;
//       case Trans:
//       switch (transA) {
//         case NoTrans:
//         C *= beta;
//         for (size_type i=0; i<A.nrow(); i++) {
//           for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
//             size_type j = A.colind(z);
//             for (index_type jx=B.rbegin(), jy=C.rbegin(); jx<=B.rend(); jx++, jy++) {
//               C(jy, i + C.cbegin()) += alpha * A.value(z) * B(jx, j + B.cbegin());
//             }
//           }
//         }
//         break;
//         case Trans:
//         C *= beta;
//         for (size_type i=0; i<A.nrow(); i++) {
//           for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
//             size_type j = A.colind(z);
//             for (index_type jx=B.rbegin(), jy=C.rbegin(); jx<=B.rend(); jx++, jy++) {
//               C(jy, j + C.cbegin()) += alpha * A.value(z) * C(jx, i + C.cbegin());
//             }
//           }
//         }
//         break;
//       }
//       break;
//     }
//     return C;
//   }
//
// #ifdef F77BLAS
//   template <>
//   inline
//   dense_matrix<double>& dgemm(
//     const trans_t& transA,
//     const trans_t& transB,
//     const double alpha,
//     const dense_matrix<double>& A,
//     const dense_matrix<double>& B,
//     const double beta,
//     dense_matrix<double>& C) {
//       switch (transA) {
//         case NoTrans:
//         switch (transB) {
//           case NoTrans:
//           assert(C.nrow() == A.nrow() && C.ncol() == B.ncol() && A.ncol() == B.nrow());
//           dblas::dgemm(transA, transB, C.nrow(), C.ncol(), A.ncol(), alpha, A.ptr(), A.ld(), B.ptr(), B.ld(), beta, C.ptr(), C.ld());
//           break;
//           case Trans:
//           assert(C.nrow() == A.nrow() && C.ncol() == B.nrow() && A.ncol() == B.ncol());
//           dblas::dgemm(transA, transB, C.nrow(), C.ncol(), A.ncol(), alpha, A.ptr(), A.ld(), B.ptr(), B.ld(), beta, C.ptr(), C.ld());
//           break;
//         }
//         break;
//         case Trans:
//         switch (transB) {
//           case NoTrans:
//           assert(C.nrow() == A.ncol() && C.ncol() == B.ncol() && A.nrow() == B.nrow());
//           dblas::dgemm(transA, transB, C.nrow(), C.ncol(), A.nrow(), alpha, A.ptr(), A.ld(), B.ptr(), B.ld(), beta, C.ptr(), C.ld());
//           break;
//           case Trans:
//           assert(C.nrow() == A.ncol() && C.ncol() == B.nrow() && A.nrow() == B.ncol());
//           dblas::dgemm(transA, transB, C.nrow(), C.ncol(), A.nrow(), alpha, A.ptr(), A.ld(), B.ptr(), B.ld(), beta, C.ptr(), C.ld());
//           break;
//         }
//         break;
//       }
//     return C;
//   }
// #endif
//
// // dgemm for vector
//
// template <typename ValueT>
// inline
// vector<ValueT>& dgemm(
//   const trans_t& transA,
//   const trans_t& transB,
//   const ValueT alpha,
//   const dense_matrix<ValueT>& A,
//   const vector<ValueT>& B,
//   const ValueT beta,
//   vector<ValueT>& C) {
//   return dgemv<ValueT>(transA, alpha, A, B, beta, C);
// }
//
// template <typename ValueT>
// inline
// vector<ValueT>& dgemm(
//   const trans_t& transA,
//   const trans_t& transB,
//   const ValueT alpha,
//   const csr_matrix<ValueT>& A,
//   const vector<ValueT>& B,
//   const ValueT beta,
//   vector<ValueT>& C) {
//   return dgemv<ValueT>(transA, alpha, A, B, beta, C);
// }
//
//   // lapack
//
//   // dgesv: solve A X = B
//   // The solution is assgined to B
//
//   template <typename ValueT>
//   dense_matrix<ValueT>& dgesv(
//     const dense_matrix<ValueT>& A, dense_matrix<ValueT>& B);
//
//   template <typename ValueT>
//   vector<ValueT>& dgesv(
//     const dense_matrix<ValueT>& A, vector<ValueT>& B);
//
// #ifdef F77BLAS
//   template <>
//   inline
//   dense_matrix<double>& dgesv(
//     const dense_matrix<double>& A, dense_matrix<double>& B) {
//     assert(A.nrow() == A.ncol());
//     dense_matrix<double> MA = A.clone();
//     array<int> ipiv(A.nrow());
//     int info = dblas::dgesv(A.nrow(), B.ncol(), MA.ptr(), MA.ld(), &ipiv[0], B.ptr(), B.ld());
//     assert(info == 0);
//     return B;
//   }
//
//   template <>
//   inline
//   vector<double>& dgesv(
//     const dense_matrix<double>& A, vector<double>& B) {
//     assert(A.nrow() == A.ncol());
//     dense_matrix<double> MA = A.clone();
//     array<int> ipiv(A.nrow());
//     int info = dblas::dgesv(A.nrow(), 1, MA.ptr(), MA.ld(), &ipiv[0], B.ptr(), B.size());
//     assert(info == 0);
//     return B;
//   }
// #endif

}
