##
##  minnndist.R
##
## Fast versions of min(nndist(X)), max(nndist(X))
##
##  $Revision: 1.11 $  $Date: 2022/05/21 09:52:11 $

minnndist <- function(X, positive=FALSE, by=NULL) {
  stopifnot(is.ppp(X))
  if(!is.null(by)) {
    stopifnot(length(by) == npoints(X))
    if(positive) {
      retain <- !duplicated(X)
      X <- X[retain]
      by <- by[retain]
    }
    nn <- nndist(X, by=by)
    result <- aggregate(nn, by=list(from=by), min, drop=FALSE)[,-1,drop=FALSE]
    return(result)
  }
  n <- npoints(X)
  if(n <= 1) return(NA)
  x <- X$x
  y <- X$y
  o <- fave.order(y)
  big <- sqrt(.Machine$double.xmax)
  if(positive) {
      z <- .C(SG_minPnnd2,
              n = as.integer(n),
              x = as.double(x[o]),
              y = as.double(y[o]),
              as.double(big),
              result = as.double(numeric(1)),
              PACKAGE="spatstat.geom")
      } else {
      z <- .C(SG_minnnd2,
              n = as.integer(n),
              x = as.double(x[o]),
              y = as.double(y[o]),
              as.double(big),
              result = as.double(numeric(1)),
              PACKAGE="spatstat.geom")
  }
  return(sqrt(z$result))
}

maxnndist <- function(X, positive=FALSE, by=NULL) {
  stopifnot(is.ppp(X))
  if(!is.null(by)) {
    stopifnot(length(by) == npoints(X))
    if(positive) {
      retain <- !duplicated(X)
      X <- X[retain]
      by <- by[retain]
    }
    nn <- nndist(X, by=by)
    result <- aggregate(nn, by=list(from=by), max, drop=FALSE)[,-1,drop=FALSE]
    return(result)
  }
  n <- npoints(X)
  if(n <= 1) return(NA)
  x <- X$x
  y <- X$y
  o <- fave.order(y)
  big <- sqrt(.Machine$double.xmax)
  if(positive) {
      z <- .C(SG_maxPnnd2,
              n = as.integer(n),
              x = as.double(x[o]),
              y = as.double(y[o]),
              as.double(big),
              result = as.double(numeric(1)),
              PACKAGE="spatstat.geom")
      } else {
      z <- .C(SG_maxnnd2,
              n = as.integer(n),
              x = as.double(x[o]),
              y = as.double(y[o]),
              as.double(big),
              result = as.double(numeric(1)),
              PACKAGE="spatstat.geom")
  }
  return(sqrt(z$result))
}

avenndist <- function(X) mean(nndist(unique(X)))

