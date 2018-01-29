
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

#' @export

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

#' @export

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

#' @export

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
