
#' @export

emest <- function(data, params, atol=1.0e-3, rtol=1.0e-6, maxiter=2000, gn=15, geps=1.0e-8) {
	t <- cumsum(c(0,data$time))
	i <- data$counts@i
	j <- data$counts@j
	n <- data$counts@x
	Cemest(t, i, j, n, params$omega, params$p, params$lam, params$a, params$mu, params$b, gn, geps, atol, rtol, maxiter)
}

mean.frm <- function(data) {
	t <- cumsum(data$time)
	n1 <- apply(data$counts, 1, sum)
	n2 <- apply(data$counts, 2, sum)
	md <- sum(t * n1) / sum(n1)
	mc <- sum(c(t,max(t)), n2) / sum(n2)
	list(md=md, mc=mc)
}

#' @export

makeinit <- function(data, shape) {
	a <- shape$a
	b <- shape$b
	m <- mean.frm(data)
	lam <- a / m$md
	mu <- b / m$md *0.9
	p <- runif(length(a) * length(b))
	p <- p / sum(p)
	list(omega=sum(data$counts), p=matrix(p, length(a), length(b)), lam=lam, a=a, mu=mu, b=b)
}

#' @export

emest.all <- function(data, phnum, atol=1.0e-3, rtol=1.0e-6, maxiter=2000, gn=15, geps=1.0e-8) {
	slist <- shape.all(phnum)
	best.llf <- -Inf
	for (s in slist) {
		cat("(a,b)=(", s$a, "|", s$b, ")")
		res <- emest(data, makeinit(data, s), atol=atol, rtol=rtol, maxiter=maxiter, gn=gn, geps=geps)
		cat(" llf=", res$llf, "\n")
		if (res$llf > best.llf) {
			best.llf <- res$llf
			best.result <- res
		}
	}
	best.result
}

# emest.all.all <- function(data, phnum, atol=1.0e-3, rtol=1.0e-6, maxiter=2000, gn=15, geps=1.0e-8) {
# 	result <- list()
# 	for (ph in 2:phnum) {
# 		tres <- system.time(res <- emest.all(data=data, phnum=ph, atol=atol, rtol=rtol, maxiter=maxiter, gn=gn, geps=geps))
# 		res <- c(res, list(phsize=ph, aic=-2*res$llf+2*2*ph, ctime=tres[1]))
# 		result <- c(result, list(res))
# 	}
# 	result
# }

#' @export

emest.inc <- function(data, phnum, shape=list(a=c(1), b=c(1)), atol=1.0e-3, rtol=1.0e-6, maxiter=2000, gn=15, geps=1.0e-8) {
	maxllf <- -Inf
	repeat {
		# cat("shape: (", shape$a, "|", shape$b, ")\n")
		shapelist <- shape.inc(shape, phnum)
		shape <- NULL
    for (s in shapelist) {
			cat("(a,b)=(", s$a, "|", s$b, ")")
			res <- emest(data, makeinit(data, s), atol=atol, rtol=rtol, maxiter=maxiter, gn=gn, geps=geps)
			cat(" llf=", res$llf, "\n")
			if (res$llf > maxllf) {
				maxllf <- res$llf
				shape <- s
				best.result <- res
			}
		}
    if (is.null(shape)) {
      break
    }
  }
	best.result
}

emest.inc.aic <- function(data, phnum, shape=list(a=c(1), b=c(1)), atol=1.0e-3, rtol=1.0e-6, maxiter=2000, gn=15, geps=1.0e-8) {
	minaic <- +Inf
  # maxllfv <- rep(-Inf, phsize)
	repeat {
		cat("shape: (", shape$a, "|", shape$b, ")\n")
		shapelist <- c(shape.inc(shape, phnum), shape.dec(shape))
		shape <- NULL
    for (s in shapelist) {
			cat("(a,b)=(", s$a, "|", s$b, ")")
			res <- emest(data, makeinit(data, s), atol=atol, rtol=rtol, maxiter=maxiter, gn=gn, geps=geps)
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
