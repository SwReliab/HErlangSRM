#include <Rcpp.h>

#include "gaussinte.hpp"

using namespace Rcpp;
using namespace marlib;

namespace herlang {

  // ff1 <- function(t0, t1, lam, a, mu, b) {
  // 	f <- function(x, tk, lam, a, mu, b) {
  // 		dgamma(x=x, shape=a, rate=lam) * pgamma(q=tk-x, shape=b, rate=mu)
  // 	}
  // 	integrate(f=f, lower=t0, upper=t1, tk=t1, a=a, b=b, lam=lam, mu=mu, abs.tol=0)$value
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

  // ff2 <- function(t0, t1, tl0, tl1, lam, a, mu, b) {
  // 	f <- function(x, tl0, tl1, lam, a, mu, b) {
  // 		dgamma(x=x, shape=a, rate=lam) * (pgamma(q=tl1-x, shape=b, rate=mu) - pgamma(q=tl0-x, shape=b, rate=mu))
  // 	}
  // 	integrate(f=f, lower=t0, upper=t1, tl0=tl0, tl1=tl1, a=a, b=b, lam=lam, mu=mu)$value
  // }

  double ff2integrand(double x, double tl0, double tl1, double lam, int a, double mu, int b) {
    return R::dgamma(x, a, 1/lam, false) * (R::pgamma(tl1-x, b, 1/mu, true, false) - R::pgamma(tl0-x, b, 1/mu, true, false));
  }

  double ff2(double t0, double t1, double tl0, double tl1, double lam, int a, double mu, int b, const gaussinte<NumericVector>& g, NumericVector fx) {
    double c = g.comp_fx(t0, t1, fx);
    for (int i=0; i<g.get_n(); i++) {
      fx[i] = ff2integrand(fx[i], tl0, tl1, lam, a, mu, b);
    }
    return g.comp_value(c, fx);
  }

  // ff3 <- function(t0, t1, tl, lam, a, mu, b) {
  // 	f <- function(x, tl, lam, a, mu, b) {
  // 		dgamma(x=x, shape=a, rate=lam) * pgamma(q=tl-x, shape=b, rate=mu, lower.tail=FALSE)
  // 	}
  // 	integrate(f=f, lower=t0, upper=t1, tl=tl, a=a, b=b, lam=lam, mu=mu)$value
  // }

  double ff3integrand(double x, double tl, double lam, int a, double mu, int b) {
    return R::dgamma(x, a, 1/lam, false) * R::pgamma(tl-x, b, 1/mu, false, false);
  }

  double ff3(double t0, double t1, double tl, double lam, int a, double mu, int b, const gaussinte<NumericVector>& g, NumericVector fx) {
    double c = g.comp_fx(t0, t1, fx);
    for (int i=0; i<g.get_n(); i++) {
      fx[i] = ff3integrand(fx[i], tl, lam, a, mu, b);
    }
    return g.comp_value(c, fx);
  }

  // ff4 <- function(t, lam, a, mu, b) {
  // 	pgamma(q=t, shape=a, rate=lam, lower.tail=FALSE)
  // }

  double ff4(double t, double lam, int a) {
    return R::pgamma(t, a, 1/lam, false, false);
  }

  ////////

  // ff1x <- function(t0, t1, lam, a, mu, b) {
  // 	f <- function(x, tk, lam, a, mu, b) {
  // 		x * dgamma(x=x, shape=a, rate=lam) * pgamma(q=tk-x, shape=b, rate=mu)
  // 	}
  // 	integrate(f=f, lower=t0, upper=t1, tk=t1, a=a, b=b, lam=lam, mu=mu)$value
  // }

  double ff1xintegrand(double x, double tk, double lam, int a, double mu, int b) {
    return x * R::dgamma(x, a, 1/lam, false) * R::pgamma(tk-x, b, 1/mu, true, false);
  }

  double ff1x(double t0, double t1, double lam, int a, double mu, int b, const gaussinte<NumericVector>& g, NumericVector fx) {
    double c = g.comp_fx(t0, t1, fx);
    for (int i=0; i<g.get_n(); i++) {
      fx[i] = ff1xintegrand(fx[i], t1, lam, a, mu, b);
    }
    return g.comp_value(c, fx);
  }

  // ff2x <- function(t0, t1, tl0, tl1, lam, a, mu, b) {
  // 	f <- function(x, tl0, tl1, lam, a, mu, b) {
  // 		x * dgamma(x=x, shape=a, rate=lam) * (pgamma(q=tl1-x, shape=b, rate=mu) - pgamma(q=tl0-x, shape=b, rate=mu))
  // 	}
  // 	integrate(f=f, lower=t0, upper=t1, tl0=tl0, tl1=tl1, a=a, b=b, lam=lam, mu=mu)$value
  // }

  double ff2xintegrand(double x, double tl0, double tl1, double lam, int a, double mu, int b) {
    return x * R::dgamma(x, a, 1/lam, false) * (R::pgamma(tl1-x, b, 1/mu, true, false) - R::pgamma(tl0-x, b, 1/mu, true, false));
  }

  double ff2x(double t0, double t1, double tl0, double tl1, double lam, int a, double mu, int b, const gaussinte<NumericVector>& g, NumericVector fx) {
    double c = g.comp_fx(t0, t1, fx);
    for (int i=0; i<g.get_n(); i++) {
      fx[i] = ff2xintegrand(fx[i], tl0, tl1, lam, a, mu, b);
    }
    return g.comp_value(c, fx);
  }

  // ff3x <- function(t0, t1, tl, lam, a, mu, b) {
  // 	f <- function(x, tl, lam, a, mu, b) {
  // 		x * dgamma(x=x, shape=a, rate=lam) * pgamma(q=tl-x, shape=b, rate=mu, lower.tail=FALSE)
  // 	}
  // 	integrate(f=f, lower=t0, upper=t1, tl=tl, a=a, b=b, lam=lam, mu=mu)$value
  // }

  double ff3xintegrand(double x, double tl, double lam, int a, double mu, int b) {
    return x * R::dgamma(x, a, 1/lam, false) * R::pgamma(tl-x, b, 1/mu, false, false);
  }

  double ff3x(double t0, double t1, double tl, double lam, int a, double mu, int b, const gaussinte<NumericVector>& g, NumericVector fx) {
    double c = g.comp_fx(t0, t1, fx);
    for (int i=0; i<g.get_n(); i++) {
      fx[i] = ff3xintegrand(fx[i], tl, lam, a, mu, b);
    }
    return g.comp_value(c, fx);
  }

  // ff4x <- function(t, lam, a, mu, b) {
  // 	a * pgamma(q=t, shape=a+1, rate=lam, lower.tail=FALSE) / lam
  // }

  double ff4x(double t, double lam, int a) {
    return a * R::pgamma(t, a+1, 1/lam, false, false) / lam;
  }

  // //////////////##

  // ff1s <- function(t0, t1, lam, a, mu, b) {
  // 	f <- function(x, tk, lam, a, mu, b) {
  // 		dgamma(x=x, shape=a, rate=lam) * b * pgamma(q=tk-x, shape=b+1, rate=mu) / mu
  // 	}
  // 	integrate(f=f, lower=t0, upper=t1, tk=t1, a=a, b=b, lam=lam, mu=mu)$value
  // }

  double ff1sintegrand(double x, double tk, double lam, int a, double mu, int b) {
    return R::dgamma(x, a, 1/lam, false) * b * R::pgamma(tk-x, b+1, 1/mu, true, false) / mu;
  }

  double ff1s(double t0, double t1, double lam, int a, double mu, int b, const gaussinte<NumericVector>& g, NumericVector fx) {
    double c = g.comp_fx(t0, t1, fx);
    for (int i=0; i<g.get_n(); i++) {
      fx[i] = ff1sintegrand(fx[i], t1, lam, a, mu, b);
    }
    return g.comp_value(c, fx);
  }

  // ff2s <- function(t0, t1, tl0, tl1, lam, a, mu, b) {
  // 	f <- function(x, tl0, tl1, lam, a, mu, b) {
  // 		dgamma(x=x, shape=a, rate=lam) * b * (pgamma(q=tl1-x, shape=b+1, rate=mu) - pgamma(q=tl0-x, shape=b+1, rate=mu)) / mu
  // 	}
  // 	integrate(f=f, lower=t0, upper=t1, tl0=tl0, tl1=tl1, a=a, b=b, lam=lam, mu=mu)$value
  // }

  double ff2sintegrand(double x, double tl0, double tl1, double lam, int a, double mu, int b) {
    return R::dgamma(x, a, 1/lam, false) * b * (R::pgamma(tl1-x, b+1, 1/mu, true, false) - R::pgamma(tl0-x, b+1, 1/mu, true, false)) / mu;
  }

  double ff2s(double t0, double t1, double tl0, double tl1, double lam, int a, double mu, int b, const gaussinte<NumericVector>& g, NumericVector fx) {
    double c = g.comp_fx(t0, t1, fx);
    for (int i=0; i<g.get_n(); i++) {
      fx[i] = ff2sintegrand(fx[i], tl0, tl1, lam, a, mu, b);
    }
    return g.comp_value(c, fx);
  }

  // ff3s <- function(t0, t1, tl, lam, a, mu, b) {
  // 	f <- function(x, tl, lam, a, mu, b) {
  // 		dgamma(x=x, shape=a, rate=lam) * b * pgamma(q=tl-x, shape=b+1, rate=mu, lower.tail=FALSE) / mu
  // 	}
  // 	integrate(f=f, lower=t0, upper=t1, tl=tl, a=a, b=b, lam=lam, mu=mu)$value
  // }

  double ff3sintegrand(double x, double tl, double lam, int a, double mu, int b) {
    return R::dgamma(x, a, 1/lam, false) * b * R::pgamma(tl-x, b+1, 1/mu, false, false) / mu;
  }

  double ff3s(double t0, double t1, double tl, double lam, int a, double mu, int b, const gaussinte<NumericVector>& g, NumericVector fx) {
    double c = g.comp_fx(t0, t1, fx);
    for (int i=0; i<g.get_n(); i++) {
      fx[i] = ff3sintegrand(fx[i], tl, lam, a, mu, b);
    }
    return g.comp_value(c, fx);
  }

  // ff4s <- function(t, lam, a, mu, b) {
  // 	b * pgamma(q=t, shape=a, rate=lam, lower.tail=FALSE) / mu
  // }

  double ff4s(double t, double lam, int a, double mu, int b) {
    return b * R::pgamma(t, a, 1/lam, false, false) / mu;
  }

  double Cestep(NumericVector t, IntegerVector i, IntegerVector j, IntegerVector n,
    double omega, NumericMatrix p, NumericVector lam, IntegerVector a, NumericVector mu, IntegerVector b,
    NumericMatrix em, NumericMatrix ex, NumericMatrix es,
    const gaussinte<NumericVector>& g, NumericVector& fx) {

    double total;
    double llf = 0;
    int K = t.length() - 1;

    NumericMatrix m(a.length(), b.length());
    NumericMatrix x(a.length(), b.length());
    NumericMatrix s(a.length(), b.length());

    for (int ii=0; ii<a.length(); ii++) {
      for (int jj=0; jj<b.length(); jj++) {
        em(ii,jj) = 0;
        ex(ii,jj) = 0;
        es(ii,jj) = 0;
      }
    }

    for (int h=0; h<n.length(); h++) {
      int nn = n[h];
      int k = i[h];
      int l = j[h];
      if (k == l) {
        total = 0;
        for (int ii=0; ii<a.length(); ii++) {
          for (int jj=0; jj<b.length(); jj++) {
            m(ii,jj) = p(ii,jj) * ff1(t[k], t[k+1], lam[ii], a[ii], mu[jj], b[jj], g, fx);
            total += m(ii,jj);
            x(ii,jj) = p(ii,jj) * ff1x(t[k], t[k+1], lam[ii], a[ii], mu[jj], b[jj], g, fx);
            s(ii,jj) = p(ii,jj) * ff1s(t[k], t[k+1], lam[ii], a[ii], mu[jj], b[jj], g, fx);
          }
        }
        for (int ii=0; ii<a.length(); ii++) {
          for (int jj=0; jj<b.length(); jj++) {
            em(ii,jj) += nn * m(ii,jj) / total;
            ex(ii,jj) += nn * x(ii,jj) / total;
            es(ii,jj) += nn * s(ii,jj) / total;
          }
        }
        llf += nn * log(omega * total) - lgamma(nn+1);
      }
      if (k != l && l != K) {
        total = 0;
        for (int ii=0; ii<a.length(); ii++) {
          for (int jj=0; jj<b.length(); jj++) {
            m(ii,jj) = p(ii,jj) * ff2(t[k], t[k+1], t[l], t[l+1], lam[ii], a[ii], mu[jj], b[jj], g, fx);
            total += m(ii,jj);
            x(ii,jj) = p(ii,jj) * ff2x(t[k], t[k+1], t[l], t[l+1], lam[ii], a[ii], mu[jj], b[jj], g, fx);
            s(ii,jj) = p(ii,jj) * ff2s(t[k], t[k+1], t[l], t[l+1], lam[ii], a[ii], mu[jj], b[jj], g, fx);
          }
        }
        for (int ii=0; ii<a.length(); ii++) {
          for (int jj=0; jj<b.length(); jj++) {
            em(ii,jj) += nn * m(ii,jj) / total;
            ex(ii,jj) += nn * x(ii,jj) / total;
            es(ii,jj) += nn * s(ii,jj) / total;
          }
        }
        llf += nn * log(omega * total) - lgamma(nn+1);
      }
      if (l == K) {
        total = 0;
        for (int ii=0; ii<a.length(); ii++) {
          for (int jj=0; jj<b.length(); jj++) {
            m(ii,jj) = p(ii,jj) * ff3(t[k], t[k+1], t[l], lam[ii], a[ii], mu[jj], b[jj], g, fx);
            total += m(ii,jj);
            x(ii,jj) = p(ii,jj) * ff3x(t[k], t[k+1], t[l], lam[ii], a[ii], mu[jj], b[jj], g, fx);
            s(ii,jj) = p(ii,jj) * ff3s(t[k], t[k+1], t[l], lam[ii], a[ii], mu[jj], b[jj], g, fx);
          }
        }
        for (int ii=0; ii<a.length(); ii++) {
          for (int jj=0; jj<b.length(); jj++) {
            em(ii,jj) += nn * m(ii,jj) / total;
            ex(ii,jj) += nn * x(ii,jj) / total;
            es(ii,jj) += nn * s(ii,jj) / total;
          }
        }
        llf += nn * log(omega * total) - lgamma(nn+1);
      }
  	}
    total = 0;
    for (int ii=0; ii<a.length(); ii++) {
      for (int jj=0; jj<b.length(); jj++) {
        m(ii,jj) = p(ii,jj) * ff4(t[K], lam[ii], a[ii]);
        total += m(ii,jj);
        x(ii,jj) = p(ii,jj) * ff4x(t[K], lam[ii], a[ii]);
        s(ii,jj) = p(ii,jj) * ff4s(t[K], lam[ii], a[ii], mu[jj], b[jj]);
      }
    }
    for (int ii=0; ii<a.length(); ii++) {
      for (int jj=0; jj<b.length(); jj++) {
        em(ii,jj) += omega * m(ii,jj);
        ex(ii,jj) += omega * x(ii,jj);
        es(ii,jj) += omega * s(ii,jj);
      }
    }
    llf -= omega * (1 - total);

    return llf;
  }

  void Cmstep(NumericMatrix em, NumericMatrix ex, NumericMatrix es,
    double& omega, NumericMatrix& p, NumericVector& lam, IntegerVector a, NumericVector& mu, IntegerVector b) {

    omega = 0;
    NumericVector lamx(a.length());
    NumericVector mus(b.length());

    for (int ii=0; ii<a.length(); ii++) {
      lam[ii] = 0;
    }
    for (int jj=0; jj<b.length(); jj++) {
      mu[jj] = 0;
    }

    for (int ii=0; ii<a.length(); ii++) {
      for (int jj=0; jj<b.length(); jj++) {
        omega += em(ii,jj);
        lam[ii] += em(ii,jj);
        mu[jj] += em(ii,jj);
        lamx[ii] += ex(ii,jj);
        mus[jj] += es(ii,jj);
      }
    }
    for (int ii=0; ii<a.length(); ii++) {
      for (int jj=0; jj<b.length(); jj++) {
        p(ii,jj) = em(ii,jj) / omega;
      }
    }
    for (int ii=0; ii<a.length(); ii++) {
      lam[ii] = a[ii] * lam[ii] / lamx[ii];
    }
    for (int jj=0; jj<b.length(); jj++) {
      mu[jj] = b[jj] * mu[jj] / mus[jj];
    }
  }
}

// [[Rcpp::export]]

List Cherlang_emest(NumericVector t, IntegerVector i, IntegerVector j, IntegerVector n,
  double omega, NumericMatrix p_, NumericVector lam_, IntegerVector a, NumericVector mu_, IntegerVector b,
  int gn = 15, double geps = 1.0e-8, double atol=1.0e-3, double rtol=1.0e-6, int maxiter=2000) {

  NumericVector gx(gn);
  NumericVector gw(gn);
  NumericVector fx(gn);
  gaussinte<NumericVector> g(gn, gx, gw, geps);

  NumericMatrix p = clone(p_);
  NumericVector lam = clone(lam_);
  NumericVector mu = clone(mu_);
  NumericMatrix em(a.length(), b.length());
  NumericMatrix ex(a.length(), b.length());
  NumericMatrix es(a.length(), b.length());

  double llf = herlang::Cestep(t, i, j, n, omega, p, lam, a, mu, b, em, ex, es, g, fx);
  herlang::Cmstep(em, ex, es, omega, p, lam, a, mu, b);
  int iter = 1;
  double aerror, rerror;
  bool conv = false;

  while (true) {
    double prev = llf;
    llf  = herlang::Cestep(t, i, j, n, omega, p, lam, a, mu, b, em, ex, es, g, fx);
    herlang::Cmstep(em, ex, es, omega, p, lam, a, mu, b);

    // if (!Rcpp::is_finite(llf)) {
    //   stop("LLF is NaN");
    // }

    if (llf - prev < 0) {
      warning("LLF decreases at iter=%d", iter);
    }

    aerror = std::abs(llf - prev);
    rerror = std::abs(aerror / llf);

    if (aerror < atol && rerror < rtol) {
      conv = true;
      break;
    }

    if (iter >= maxiter) {
      break;
    }

    iter++;
    checkUserInterrupt();
}

  return List::create(
    Named("omega") = omega,
    Named("p") = p,
    Named("lam") = lam,
    Named("a") = a,
    Named("mu") = mu,
    Named("b") = b,
    Named("llf") = llf,
    Named("aerror") = aerror,
    Named("rerror") = rerror,
    Named("iter") = iter
  );
}

// [[Rcpp::export]]

NumericVector Cherlang_mvfd(NumericVector t, List model) {
  const double omega = model["omega"];
  NumericMatrix p = model["p"];
  NumericVector lam = model["lam"];
  IntegerVector a = model["a"];
  NumericVector mu = model["mu"];
  IntegerVector b = model["b"];

  NumericVector result(t.length());
  NumericVector pa(a.length());
  for (int ii=0; ii<a.length(); ii++) {
    for (int jj=0; jj<b.length(); jj++) {
      pa[ii] += p(ii,jj);
    }
  }
  for (int k=0; k<t.length(); k++) {
    for (int ii=0; ii<a.length(); ii++) {
      result[k] += omega * pa[ii] * R::pgamma(t[k], a[ii], 1/lam[ii], true, false);
    }
  }
  return result;
}

// [[Rcpp::export]]

NumericVector Cherlang_mvfc(NumericVector t, List model, int gn = 15, double geps = 1.0e-8) {
  const double omega = model["omega"];
  NumericMatrix p = model["p"];
  NumericVector lam = model["lam"];
  IntegerVector a = model["a"];
  NumericVector mu = model["mu"];
  IntegerVector b = model["b"];

  NumericVector gx(gn);
  NumericVector gw(gn);
  NumericVector fx(gn);
  gaussinte<NumericVector> g(gn, gx, gw, geps);

  NumericVector result(t.length());
  double prevt = 0;
  double prevx = 0;
  for (int k=0; k<t.length(); k++) {
    for (int ii=0; ii<a.length(); ii++) {
      for (int jj=0; jj<b.length(); jj++) {
        result[k] += omega * p(ii,jj) * herlang::ff1(prevt, t[k], lam[ii], a[ii], mu[jj], b[jj], g, fx);
      }
    }
    for (int l=k+1; l<t.length(); l++) {
      for (int ii=0; ii<a.length(); ii++) {
        for (int jj=0; jj<b.length(); jj++) {
          result[l] += omega * p(ii,jj) * herlang::ff2(prevt, t[k], t[l-1], t[l], lam[ii], a[ii], mu[jj], b[jj], g, fx);
        }
      }
    }
    prevt = t[k];
    result[k] += prevx;
    prevx = result[k];
  }
  return result;
}

// reli <- function(t, s, res) {
//   p <- apply(res$p, 1, sum)
//   v <- numeric(length(t))
//   for (i in 1:length(res$a)) {
//     v <- v + p[i] * (pgamma(q=t, shape=res$a[i], rate=res$lam[i]) - pgamma(q=s, shape=res$a[i], rate=res$lam[i]))
//  }
//  exp(-res$omega * v)
// }

// [[Rcpp::export]]

NumericVector Cherlang_reli(NumericVector t, double s, List model) {
  const double omega = model["omega"];
  NumericMatrix p = model["p"];
  NumericVector lam = model["lam"];
  IntegerVector a = model["a"];
  NumericVector mu = model["mu"];
  IntegerVector b = model["b"];

  NumericVector result(t.length());
  NumericVector pa(a.length());
  for (int ii=0; ii<a.length(); ii++) {
    for (int jj=0; jj<b.length(); jj++) {
      pa[ii] += p(ii,jj);
    }
  }
  for (int k=0; k<t.length(); k++) {
    for (int ii=0; ii<a.length(); ii++) {
      result[k] += pa[ii] * (R::pgamma(t[k], a[ii], 1/lam[ii], true, false) - R::pgamma(s, a[ii], 1/lam[ii], true, false));
    }
    result[k] = exp(-omega * result[k]);
  }
  return result;
}

// [[Rcpp::export]]

NumericVector Cherlang_pc(NumericVector t, double s, List model, int gn = 15, double geps = 1.0e-8) {
  const double omega = model["omega"];
  NumericMatrix p = model["p"];
  NumericVector lam = model["lam"];
  IntegerVector a = model["a"];
  NumericVector mu = model["mu"];
  IntegerVector b = model["b"];

  NumericVector gx(gn);
  NumericVector gw(gn);
  NumericVector fx(gn);
  gaussinte<NumericVector> g(gn, gx, gw, geps);

  NumericVector result(t.length());
  double prevt = s;
  double prevx = 0;
  for (int k=0; k<t.length(); k++) {
    for (int ii=0; ii<a.length(); ii++) {
      for (int jj=0; jj<b.length(); jj++) {
        result[k] += p(ii,jj) * herlang::ff1(prevt, t[k], lam[ii], a[ii], mu[jj], b[jj], g, fx);
      }
    }
    for (int l=k+1; l<t.length(); l++) {
      for (int ii=0; ii<a.length(); ii++) {
        for (int jj=0; jj<b.length(); jj++) {
          result[l] += p(ii,jj) * herlang::ff2(prevt, t[k], t[l-1], t[l], lam[ii], a[ii], mu[jj], b[jj], g, fx);
        }
      }
    }
    prevt = t[k];
    result[k] += prevx;
    prevx = result[k];
    result[k] = exp(omega * result[k]);
  }
  return result;
}

// relic <- function(t, s, res) {
// 	reli(t, s, res) * ppp(t, s, res)
// }

// [[Rcpp::export]]

NumericVector Cherlang_relic(NumericVector t, double s, List model, int gn = 15, double geps = 1.0e-8) {
  NumericVector r = Cherlang_reli(t, s, model);
  NumericVector p = Cherlang_pc(t, s, model, gn, geps);
  for (int k=0; k<t.length(); k++) {
    r[k] *= p[k];
  }
  return r;
}
