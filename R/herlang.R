
#' @export

cor <- function(model) {
	xmean <- sum((model$a / model$lam) * apply(model$p, 1, sum))
	x2mean <- sum(((model$a + model$a^2) / model$lam^2) * apply(model$p, 1, sum))
	xvar <- x2mean - xmean^2
	ymean <- sum((model$b / model$mu) * apply(model$p, 2, sum))
	y2mean <- sum(((model$b + model$b^2) / model$mu^2) * apply(model$p, 2, sum))
	yvar <- y2mean - ymean^2
	xymean <- sum(model$p * matrix(model$a / model$lam, ncol=1) %*% matrix(model$b / model$mu, nrow=1))
	cov <- xymean - xmean * ymean
	rho <- cov / sqrt(xvar * yvar)
	list(xmean=xmean, xvar=xvar, ymean=ymean, yvar=yvar, cov=cov, rho=rho)
}

#' @export

mvfd <- function(t, model) {
	Cmvfd(t, model)
}

#' @export

mvfc <- function(t, model, gn = 15L, geps = 1.0e-8) {
	Cmvfc(t, model, gn, geps)
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
