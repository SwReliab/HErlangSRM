
tofrm <- function(dd) {
	nmax <- dim(dd$counts)[1]
	t <- c(0, cumsum(dd$days))
	mat <- as(dd$counts, "TsparseMatrix")
	list(ij=cbind(mat@i+1,mat@j+1), n=mat@x, nmax=nmax, t=t)
}

### Erlang

makeErlang <- function(lam, a, mu, b) {
	i <- 1:a
	j <- 1:a
	x <- rep(-lam,a)

	i <- c(i, a + 1:b)
	j <- c(j, a + 1:b)
	x <- c(x, rep(-mu,b))

	i <- c(i, 1:a)
	j <- c(j, 2:(a+1))
	x <- c(x, rep(lam,a))

	if (b >= 2) {
		i <- c(i, a + 1:(b-1))
		j <- c(j, a + 1 + 1:(b-1))
		x <- c(x, rep(mu, b-1))
	}

	T <- sparseMatrix(i=i, j=j, x=x)
	tau <- -apply(T, 1, sum)

	list(alpha=c(1, rep(0, a+b-1)), T=T, tau=tau)
}

#' @export

ff1 <- function(t0, t1, lam, a, mu, b) {
#	cat("ff1", t0, t1, lam, a, mu, b, "\n")
	f <- function(x, tk, lam, a, mu, b) {
		dgamma(x=x, shape=a, rate=lam) * pgamma(q=tk-x, shape=b, rate=mu)
	}
	res0 <- integrate(f=f, lower=t0, upper=t1, tk=t1, a=a, b=b, lam=lam, mu=mu, abs.tol=0)$value
#	res <- deformula.moneone(f=f, lower=t0, upper=t1, tk=t1, a=a, b=b, lam=lam, mu=mu)$value
#	cat(res0, "\n")
	res0
}

ff2 <- function(t0, t1, tl0, tl1, lam, a, mu, b) {
	f <- function(x, tl0, tl1, lam, a, mu, b) {
		dgamma(x=x, shape=a, rate=lam) * (pgamma(q=tl1-x, shape=b, rate=mu) - pgamma(q=tl0-x, shape=b, rate=mu))
	}
	integrate(f=f, lower=t0, upper=t1, tl0=tl0, tl1=tl1, a=a, b=b, lam=lam, mu=mu)$value
}

ff3 <- function(t0, t1, tl, lam, a, mu, b) {
	f <- function(x, tl, lam, a, mu, b) {
		dgamma(x=x, shape=a, rate=lam) * pgamma(q=tl-x, shape=b, rate=mu, lower.tail=FALSE)
	}
	integrate(f=f, lower=t0, upper=t1, tl=tl, a=a, b=b, lam=lam, mu=mu)$value
}

ff4 <- function(t, lam, a, mu, b) {
	pgamma(q=t, shape=a, rate=lam, lower.tail=FALSE)
}

##

ff1x <- function(t0, t1, lam, a, mu, b) {
	f <- function(x, tk, lam, a, mu, b) {
		x * dgamma(x=x, shape=a, rate=lam) * pgamma(q=tk-x, shape=b, rate=mu)
	}
	integrate(f=f, lower=t0, upper=t1, tk=t1, a=a, b=b, lam=lam, mu=mu)$value
}

ff2x <- function(t0, t1, tl0, tl1, lam, a, mu, b) {
	f <- function(x, tl0, tl1, lam, a, mu, b) {
		x * dgamma(x=x, shape=a, rate=lam) * (pgamma(q=tl1-x, shape=b, rate=mu) - pgamma(q=tl0-x, shape=b, rate=mu))
	}
	integrate(f=f, lower=t0, upper=t1, tl0=tl0, tl1=tl1, a=a, b=b, lam=lam, mu=mu)$value
}

ff3x <- function(t0, t1, tl, lam, a, mu, b) {
	f <- function(x, tl, lam, a, mu, b) {
		x * dgamma(x=x, shape=a, rate=lam) * pgamma(q=tl-x, shape=b, rate=mu, lower.tail=FALSE)
	}
	integrate(f=f, lower=t0, upper=t1, tl=tl, a=a, b=b, lam=lam, mu=mu)$value
}

ff4x <- function(t, lam, a, mu, b) {
	a * pgamma(q=t, shape=a+1, rate=lam, lower.tail=FALSE) / lam
}

##

ff1s <- function(t0, t1, lam, a, mu, b) {
	f <- function(x, tk, lam, a, mu, b) {
		dgamma(x=x, shape=a, rate=lam) * b * pgamma(q=tk-x, shape=b+1, rate=mu) / mu
	}
	integrate(f=f, lower=t0, upper=t1, tk=t1, a=a, b=b, lam=lam, mu=mu)$value
}

ff2s <- function(t0, t1, tl0, tl1, lam, a, mu, b) {
	f <- function(x, tl0, tl1, lam, a, mu, b) {
		dgamma(x=x, shape=a, rate=lam) * b * (pgamma(q=tl1-x, shape=b+1, rate=mu) - pgamma(q=tl0-x, shape=b+1, rate=mu)) / mu
	}
	integrate(f=f, lower=t0, upper=t1, tl0=tl0, tl1=tl1, a=a, b=b, lam=lam, mu=mu)$value
}

ff3s <- function(t0, t1, tl, lam, a, mu, b) {
	f <- function(x, tl, lam, a, mu, b) {
		dgamma(x=x, shape=a, rate=lam) * b * pgamma(q=tl-x, shape=b+1, rate=mu, lower.tail=FALSE) / mu
	}
	integrate(f=f, lower=t0, upper=t1, tl=tl, a=a, b=b, lam=lam, mu=mu)$value
}

ff4s <- function(t, lam, a, mu, b) {
	b * pgamma(q=t, shape=a, rate=lam, lower.tail=FALSE) / mu
}

##

estep <- function(frm, omega, p, lam, a, mu, b) {
	K <- frm$nmax
	t <- frm$t
	llf <- 0
	indexlist <- makelist(length(a), length(b))
	em <- array(0, length(indexlist))
	ex <- array(0, length(indexlist))
	es <- array(0, length(indexlist))
	for (h in 1:length(frm$n)) {
		K <- frm$nmax
		nn <- frm$n[h]
		k <- frm$ij[h,1]
		l <- frm$ij[h,2]
#cat(c(k,l), "\n")
		if (k == l) {
#cat("p", p, "\n")
			u <- sapply(indexlist, function(idx) p[idx[1],idx[2]] * ff1(t[k], t[k+1], lam[idx[1]], a[idx[1]], mu[idx[2]], b[idx[2]]))
			total <- sum(u)
			x <- sapply(indexlist, function(idx) p[idx[1],idx[2]] * ff1x(t[k], t[k+1], lam[idx[1]], a[idx[1]], mu[idx[2]], b[idx[2]]))
			s <- sapply(indexlist, function(idx) p[idx[1],idx[2]] * ff1s(t[k], t[k+1], lam[idx[1]], a[idx[1]], mu[idx[2]], b[idx[2]]))
#cat("k=", k, "u=", u, "x=", x, "s=", s, "\n")
			em <- em + nn * u / total
			ex <- ex + nn * x / total
			es <- es + nn * s / total
			llf <- llf + nn * log(omega * total) - lgamma(nn+1)
		}
		if (k != l && l != K+1) {
			u <- sapply(indexlist, function(idx) p[idx[1],idx[2]] * ff2(t[k], t[k+1], t[l], t[l+1], lam[idx[1]], a[idx[1]], mu[idx[2]], b[idx[2]]))
			total <- sum(u)
			x <- sapply(indexlist, function(idx) p[idx[1],idx[2]] * ff2x(t[k], t[k+1], t[l], t[l+1], lam[idx[1]], a[idx[1]], mu[idx[2]], b[idx[2]]))
			s <- sapply(indexlist, function(idx) p[idx[1],idx[2]] * ff2s(t[k], t[k+1], t[l], t[l+1], lam[idx[1]], a[idx[1]], mu[idx[2]], b[idx[2]]))
#cat("k=", k, "l=", l, "u=", u, "x=", x, "s=", s, "\n")
			em <- em + nn * u / total
			ex <- ex + nn * x / total
			es <- es + nn * s / total
			llf <- llf + nn * log(omega * total) - lgamma(nn+1)
		}
		if (l == K+1) {
			u <- sapply(indexlist, function(idx) p[idx[1],idx[2]] * ff3(t[k], t[k+1], t[l], lam[idx[1]], a[idx[1]], mu[idx[2]], b[idx[2]]))
			total <- sum(u)
			x <- sapply(indexlist, function(idx) p[idx[1],idx[2]] * ff3x(t[k], t[k+1], t[l], lam[idx[1]], a[idx[1]], mu[idx[2]], b[idx[2]]))
			s <- sapply(indexlist, function(idx) p[idx[1],idx[2]] * ff3s(t[k], t[k+1], t[l], lam[idx[1]], a[idx[1]], mu[idx[2]], b[idx[2]]))
#cat("k=", k, "l=", l, "u=", u, "x=", x, "s=", s, "\n")
			em <- em + nn * u / total
			ex <- ex + nn * x / total
			es <- es + nn * s / total
			llf <- llf + nn * log(omega * total) - lgamma(nn+1)
		}
	}
	u <- sapply(indexlist, function(idx) p[idx[1],idx[2]] * ff4(t[K+1], lam[idx[1]], a[idx[1]], mu[idx[2]], b[idx[2]]))
	x <- sapply(indexlist, function(idx) p[idx[1],idx[2]] * ff4x(t[K+1], lam[idx[1]], a[idx[1]], mu[idx[2]], b[idx[2]]))
	s <- sapply(indexlist, function(idx) p[idx[1],idx[2]] * ff4s(t[K+1], lam[idx[1]], a[idx[1]], mu[idx[2]], b[idx[2]]))
#cat("last", "k=", K+1, "u=", u, "x=", x, "s=", s, "\n")
	total <- sum(u)
	em <- em + omega * u
	ex <- ex + omega * x
	es <- es + omega * s
	llf <- llf - omega * (1 - total)

	list(em=matrix(em, length(a), length(b)), ex=matrix(ex, length(a), length(b)), es=matrix(es, length(a), length(b)), a=a, b=b, llf=llf)
}

mstep <- function(eres) {
	em <- eres$em
	ex <- eres$ex
	es <- eres$es
	a <- eres$a
	b <- eres$b

	omega <- sum(em)
	p <- em / omega
	lam <- a * apply(em, 1, sum) / apply(ex, 1, sum)
	mu <- b * apply(em, 2, sum) / apply(es, 2, sum)

	list(omega=omega, p=p, lam=lam, mu=mu)
}

emest.param <- function(frm, param, atol=1.0e-3, rtol=1.0e-6, maxiter=2000) {
	res <- emest(frm, param$omega, param$p, param$lam, param$a, param$mu, param$b, atol=atol, rtol=rtol, maxiter=maxiter)
	c(res, list(a=param$a, b=param$b))
}

emest <- function(frm, omega, p, lam, a, mu, b, atol=1.0e-3, rtol=1.0e-6, maxiter=2000) {

	if (!is.matrix(p)) {
		p <- matrix(p, length(a), length(b))
	}

	eres <- estep(frm, omega, p, lam, a, mu, b)
#	cat(sprintf("llf %f\n", eres$llf))

	prev <- eres$llf
	param <- mstep(eres)

	iter <- 1
	conv <- FALSE
	repeat {
		eres <- estep(frm, param$omega, param$p, param$lam, a, param$mu, b)
		param <- mstep(eres)

		if (!is.finite(eres$llf)) {
			warning("LLF is NaN")
			eres$llf <- -Inf
			break
		}

		if (eres$llf - prev < 0) {
			warning("LLF decrease")
		}

		aerror <- abs(eres$llf - prev)
		rerror <- abs(aerror / eres$llf)

#		cat(sprintf("llf %f (aerror=%e, rerror=%e)\n", eres$llf, aerror, rerror))
#		cat(sprintf("param: omega=%f lam=%e mu=%e\n", param$omega, param$lam, param$mu))

		if (aerror < atol && rerror < rtol) {
			conv <- TRUE
			break
		}

		if (iter >= maxiter) {
			break
		}

		prev <- eres$llf
		iter <- iter + 1
	}

	c(param, list(llf=eres$llf, aerror=aerror, rerror=rerror, atol=atol, rtol=rtol, iter=iter, maxiter=maxiter))
}

herlang.shape.all <- function(phnum, lbound, ubound) {
  herlang.shape.ipart(1, numeric(phnum), 1, phnum, lbound, ubound, list())
}

herlang.shape.ipart <- function(pos, shape, mini, res, lb, ub, result) {
  if (mini > res) {
    return(result)
  } else {
    for (i in mini:res) {
      shape[pos] <- i
      result <- herlang.shape.ipart(pos+1, shape, i, res-i, lb, ub, result)
    }
    if (lb <= pos && pos <= ub) {
      result <- c(result, list(shape[1:pos]))
    }
    return(result)
  }
}

shape.all <- function(phnum) {
	result <- list()
	for (p in 1:(phnum-1)) {
		dp <- herlang.shape.all(p, 1, p)
		cp <- herlang.shape.all(phnum-p, 1, phnum-p)
		for (xd in dp) {
			for (xc in cp) {
				result <- c(result, list(list(a=xd, b=xc)))
			}
		}
	}
	result
}

mean.frm <- function(frm) {
	md <- sum(frm$t[frm$ij[,1]] * frm$n) / sum(frm$n)
	mc <- sum(frm$t[frm$ij[,2]] * frm$n) / sum(frm$n)
	list(md=md, mc=mc)
}

makeinit <- function(frm, shape) {
	a <- shape$a
	b <- shape$b
	m <- mean.frm(frm)
	lam <- a / m$md
	mu <- b / m$md *0.9
	p <- runif(length(a) * length(b))
	p <- p / sum(p)
	list(omega=sum(frm$n), p=matrix(p, length(a), length(b)), lam=lam, a=a, mu=mu, b=b)
}

emest.all <- function(frm, phnum, atol=1.0e-3, rtol=1.0e-6, maxiter=2000) {
	slist <- shape.all(phnum)
	best.llf <- -Inf
	for (s in slist) {
		print(makeinit(frm, s))
		res <- emest.param(frm, makeinit(frm, s), atol=atol, rtol=rtol, maxiter=maxiter)
		if (res$llf > best.llf) {
			best.llf <- res$llf
#			best.result <- c(res, list(a=s$a, b=s$b))
			best.result <- res
		}
	}
	best.result
}

emest.all.all <- function(frm, phnum, atol=1.0e-3, rtol=1.0e-6, maxiter=2000) {
	result <- list()
	for (ph in 2:phnum) {
		tres <- system.time(res <- emest.all(frm=frm, phnum=ph, atol=atol, rtol=rtol, maxiter=maxiter))
		res <- c(res, list(phsize=ph, aic=-2*res$llf+2*2*ph, ctime=tres[1]))
		result <- c(result, list(res))
	}
	result
}

## inc

herlang.shape.increment <- function(shape, ubound, phsize) {
  ## append
  if (length(shape) == ubound || sum(shape) == phsize) {
    return(list())
  }
  result <- c(1, shape)
  retval <- list(result)

  ## add
  for (i in unique(shape)) {
    tmp <- shape
    l <- which(tmp == i)
    n <- length(l)
    tmp[l[n]] <- tmp[l[n]] + 1
    retval <- c(retval, list(tmp))
  }
  return(retval)
}

herlang.shape.decrement <- function(shape, lbound) {
  if (length(shape)==1 && shape[1]==1) {
    return(list())
  }
  retval <- list()
  ## subtract
  for (i in unique(shape)) {
    tmp <- shape
    l <- which(tmp == i)
    n <- length(l)
    tmp[l[1]] <- tmp[l[1]] - 1
    tmp <- tmp[tmp != 0]
    if (length(tmp) >= lbound)
      retval <- c(retval, list(tmp))
  }
  return(retval)
}

shape.inc <- function(s, phsize) {
	result <- list()
	dp <- herlang.shape.increment(s$a, phsize - sum(s$b), phsize - sum(s$b))
	for (xd in dp) {
		result <- c(result, list(list(a=xd, b=s$b)))
	}
	cp <- herlang.shape.increment(s$b, phsize - sum(s$a), phsize - sum(s$a))
	for (xc in cp) {
		result <- c(result, list(list(a=s$a, b=xc)))
	}
	result
}

shape.dec <- function(s) {
	result <- list()
	dp <- herlang.shape.decrement(s$a, 1)
	for (xd in dp) {
		result <- c(result, list(list(a=xd, b=s$b)))
	}
	cp <- herlang.shape.decrement(s$b, 1)
	for (xc in cp) {
		result <- c(result, list(list(a=s$a, b=xc)))
	}
	result
}

emest.inc <- function(frm, phnum, shape=list(a=c(1), b=c(1)), atol=1.0e-3, rtol=1.0e-6, maxiter=2000) {
	maxllf <- -Inf
  # maxllfv <- rep(-Inf, phsize)
  # shape <- list(a=c(1), b=c(1))
	repeat {
		cat("shape: (", shape$a, "|", shape$b, ")\n")
#		shapelist <- c(shape.inc(shape, phnum), shape.dec(shape))
		shapelist <- shape.inc(shape, phnum)
		shape <- NULL
    for (s in shapelist) {
			cat("(a,b)=(", s$a, "|", s$b, ")")
			res <- emest.param(frm, makeinit(frm, s), atol=atol, rtol=rtol, maxiter=maxiter)
			cat(" llf=", res$llf, "\n")
			if (res$llf > maxllf) {
				maxllf <- res$llf
	#			best.result <- c(res, list(a=s$a, b=s$b))
				shape <- s
				best.result <- res
			}
		}
    if (is.null(shape)) {
##      warning(message="break")
      break
    }
  }
	best.result
}

emest.inc.aic <- function(frm, phnum, shape=list(a=c(1), b=c(1)), atol=1.0e-3, rtol=1.0e-6, maxiter=2000) {
	minaic <- +Inf
  # maxllfv <- rep(-Inf, phsize)
	repeat {
		cat("shape: (", shape$a, "|", shape$b, ")\n")
		shapelist <- c(shape.inc(shape, phnum), shape.dec(shape))
		shape <- NULL
    for (s in shapelist) {
			cat("(a,b)=(", s$a, "|", s$b, ")")
			res <- emest.param(frm, makeinit(frm, s), atol=atol, rtol=rtol, maxiter=maxiter)
			aic=-2*res$llf+2*2*(sum(s$a)+sum(s$b))
			cat("aic=", aic, "\n")
			if (aic < minaic) {
				minaic <- aic
	#			best.result <- c(res, list(a=s$a, b=s$b))
				shape <- s
				best.result <- c(res, list(aic=aic))
			}
		}
    if (is.null(shape)) {
##      warning(message="break")
      break
    }
  }
	best.result
}

###

cor <- function(res) {
	xmean <- sum((res$a / res$lam) * apply(res$p, 1, sum))
	x2mean <- sum(((res$a + res$a^2) / res$lam^2) * apply(res$p, 1, sum))
	xvar <- x2mean - xmean^2
	ymean <- sum((res$b / res$mu) * apply(res$p, 2, sum))
	y2mean <- sum(((res$b + res$b^2) / res$mu^2) * apply(res$p, 2, sum))
	yvar <- y2mean - ymean^2
	xymean <- sum(res$p * matrix(res$a / res$lam, ncol=1) %*% matrix(res$b / res$mu, nrow=1))
	cov <- xymean - xmean * ymean
	rho <- cov / sqrt(xvar * yvar)
	list(xmean=xmean, xvar=xvar, ymean=ymean, yvar=yvar, cov=cov, rho=rho)
}

mvfd <- function(t, res) {
  p <- apply(res$p, 1, sum)
  v <- numeric(length(t))
  for (i in 1:length(res$a)) {
    v <- v + p[i] * pgamma(q=t, shape=res$a[i], rate=res$lam[i])
 }
 res$omega * v
}

mvfc <- function(t, res) {
  v <- 0
  for (i in 1:length(res$a)) {
  	for (j in 1:length(res$b)) {
  erlangQ <- makeErlang(res$lam[i], res$a[i], res$mu[j], res$b[j])
  tmp <- sapply(t, function(x) sum(mexpAv.pade(A=erlangQ$T, v=erlangQ$alpha, t=x, transpose=TRUE)))
  v <- v + res$p[i,j] * (1- tmp)
  	}
 }
 res$omega * v
}

reli <- function(t, s, res) {
  p <- apply(res$p, 1, sum)
  v <- numeric(length(t))
  for (i in 1:length(res$a)) {
    v <- v + p[i] * (pgamma(q=t, shape=res$a[i], rate=res$lam[i]) - pgamma(q=s, shape=res$a[i], rate=res$lam[i]))
 }
 exp(-res$omega * v)
}

ppp <- function(t, s, res) {
	indexlist <- makelist(length(res$a), length(res$b))
	v <- sapply(t, function(tt) sum(sapply(indexlist, function(idx) res$p[idx[1],idx[2]] * ff1(s, tt, res$lam[idx[1]], res$a[idx[1]], res$mu[idx[2]], res$b[idx[2]]))))
	exp(res$omega * v)
}

relic <- function(t, s, res) {
	reli(t, s, res) * ppp(t, s, res)
}

###

logreli <- function(t, s, res) {
  p <- apply(res$p, 1, sum)
  v <- numeric(length(t))
  for (i in 1:length(res$a)) {
    v <- v + p[i] * (pgamma(q=t, shape=res$a[i], rate=res$lam[i]) - pgamma(q=s, shape=res$a[i], rate=res$lam[i]))
 }
 -res$omega * v
}

logppp <- function(t, s, res) {
	indexlist <- makelist(length(res$a), length(res$b))
	v <- sapply(t, function(tt) sum(sapply(indexlist, function(idx) res$p[idx[1],idx[2]] * ff1(s, tt, res$lam[idx[1]], res$a[idx[1]], res$mu[idx[2]], res$b[idx[2]]))))
	res$omega * v
}

logrelic <- function(t, s, res) {
	logreli(t, s, res) + logppp(t, s, res)
}
