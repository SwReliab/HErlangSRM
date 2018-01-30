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
	Cherlang_mvfd(t, model)
}

#' @export

mvfc <- function(t, model, gn = 15L, geps = 1.0e-8) {
	Cherlang_mvfc(t, model, gn, geps)
}

#' @export

reli <- function(t, s, model) {
	Cherlang_reli(t, s, model)
}

#' @export

pc <- function(t, s, model, gn = 15L, geps = 1.0e-8) {
	Cherlang_pc(t, s, model, gn, geps)
}

#' @export

relic <- function(t, s, model, gn = 15L, geps = 1.0e-8) {
	Cherlang_relic(t, s, model, gn, geps)
}
