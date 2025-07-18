#
#   nndist.R
#
#   nearest neighbour distances (nndist) and identifiers (nnwhich)
#
#   $Revision: 1.19 $ $Date: 2025/06/16 05:35:37 $
#

nndist <- function(X, ...) {
  UseMethod("nndist")
}

nndist.ppp <- local({

  nndist.ppp <- function(X, ..., k=1, by=NULL,
                         proper=FALSE, method="C", metric=NULL) {
    verifyclass(X, "ppp")
    trap.extra.arguments(..., .Context="In nndist.ppp")
    if(is.null(by)) {
      if(!proper) {
        ## allow duplicated points
        if(is.null(metric)) {
          ## default case: Euclidean distance
          d <- nndist.default(X$x, X$y, k=k, method=method)
        } else {
          ## some other distance metric
          d <- invoke.metric(metric, "nndist.ppp",
                             X, ..., k=k, method=method)
        }
      } else {
        ## collapse duplicated points
        um <- uniquemap(unmark(X))
        isuniq <- (um == seq_along(um))
        UX <- X[isuniq]
        ## compute nn distances between unique points
        if(is.null(metric)) {
          dU <- nndist.default(UX$x, UX$y, k=k, method=method)
        } else {
          dU <- invoke.metric(metric, "nndist.ppp",
                              UX, ..., k=k, method=method)
        }
        ## remap to original dataset
        d <- marksubset(dU, cumsum(isuniq))
      }
    } else {
      ## split by a factor
      if(is.character(by)) {
        ## Interpret using split.ppp
        Y <- split(X, f=by, drop=FALSE)
        by <- attr(Y, "fgroup")
      }
      idX <- seq_len(npoints(X)) 
      Y <- split(X %mark% idX, f=by, un=FALSE)
      if(proper) {
        ## merge duplicates, but only within a group
        Xid <- lapply(Y, marks)
        Zid <- lapply(lapply(Y, unmark), uniquemap)
        Yid <- mapply("[", x=Xid, i=Zid, SIMPLIFY=FALSE)
        Y <- mapply("marks<-", x=Y, value=Yid, SIMPLIFY=FALSE)
        idX <- unsplit(Yid, f=by)
      }
      distY <- lapply(Y, nndistsub, XX=X, iX=idX, k=k, metric=metric)
      d <- do.call(cbind, distY)
    }
    return(d)
  }
  
  nndistsub <- function(Z, XX, iX, k, metric=NULL) {
    if(is.null(metric)) {
      d <- nncross(XX, Z, iX=iX, iY=marks(Z), k=k, what="dist")
    } else {
      d <- invoke.metric(metric, "nncross.ppp",
                    X=XX, Y=Z, iX=iX, iY=marks(Z), k=k, what="dist")
    }
    return(d)
  }

  nndist.ppp
})

nndist.default <- function(X, Y=NULL, ..., k=1, by=NULL, method="C")
{
  warn.no.metric.support("nndist.default", ...)
	#  computes the vector of nearest-neighbour distances 
	#  for the pattern of points (x[i],y[i])
	#
  xy <- xy.coords(X,Y)[c("x","y")]
  x <- xy$x
  y <- xy$y

  # validate
  n <- length(x)
  if(length(y) != n)
    stop("lengths of x and y do not match")

  method <- match.arg(method, c("C", "interpreted", "test"))
  # other arguments ignored
  trap.extra.arguments(..., .Context="In nndist.default")

  # split by factor ?
  if(!is.null(by)) {
    X <- as.ppp(xy, W=boundingbox)
    return(nndist(X, by=by, k=k))
  }
  
  # k can be a single integer or an integer vector
  if(length(k) == 0)
    stop("k is an empty vector")
  else if(length(k) == 1) {
    if(k != round(k) || k <= 0)
      stop("k is not a positive integer")
  } else {
    if(any(k != round(k)) || any(k <= 0))
      stop(paste("some entries of the vector",
           sQuote("k"), "are not positive integers"))
  }
  k <- as.integer(k)
  kmax <- max(k)

  # trivial cases
  if(n <= 1) {
    # empty pattern => return numeric(0)
    # or pattern with only 1 point => return Inf
    nnd <- matrix(Inf, nrow=n, ncol=kmax)
    nnd <- nnd[,k, drop=TRUE]
    return(nnd)
  }
  
  # number of neighbours that are well-defined
  kmaxcalc <- min(n-1, kmax)

  # calculate k-nn distances for k <= kmaxcalc
  
  if(kmaxcalc == 1) {
    # calculate nearest neighbour distance only
    switch(method,
         test = ,
         interpreted={
           #  matrix of squared distances between all pairs of points
           sq <- function(a, b) { (a-b)^2 }
           squd <-  outer(x, x, sq) + outer(y, y, sq)
           #  reset diagonal to a large value so it is excluded from minimum
           diag(squd) <- Inf
           #  nearest neighbour distances
           nnd <- sqrt(apply(squd,1,min))
         },
         C={
           nnd<-numeric(n)
           o <- fave.order(y)
           big <- sqrt(.Machine$double.xmax)
           z <- .C(SG_nndistsort,
                  n= as.integer(n),
                  x= as.double(x[o]), y= as.double(y[o]), nnd= as.double(nnd),
                  as.double(big),
                  PACKAGE="spatstat.geom")
           nnd[o] <- z$nnd
         },
         stop(paste("Unrecognised method", sQuote(method)))
         )
  } else {
    # case kmaxcalc > 1
    switch(method,
           test = ,
           interpreted={
             if(n <= 1000 && method == "interpreted") {
               # form n x n matrix of squared distances
               D2 <- pairdist.default(x, y, method=method, squared=TRUE)
               # find k'th smallest squared distance
               diag(D2) <- Inf
               NND2 <- t(apply(D2, 1, sort))[, 1:kmaxcalc]
               nnd <- sqrt(NND2)
             } else {
               # avoid creating huge matrix
               # handle one row of D at a time
               NND2 <- matrix(numeric(n * kmaxcalc), nrow=n, ncol=kmaxcalc)
               for(i in seq_len(n)) {
                 D2i <- (x - x[i])^2 + (y - y[i])^2
                 D2i[i] <- Inf
                 NND2[i,] <- orderstats(D2i, k=1:kmaxcalc)
               }
               nnd <- sqrt(NND2)
             }
           },
           C={
             nnd<-numeric(n * kmaxcalc)
             o <- fave.order(y)
             big <- sqrt(.Machine$double.xmax)
             z <- .C(SG_knndsort,
                    n    = as.integer(n),
                    kmax = as.integer(kmaxcalc),
                    x    = as.double(x[o]),
                    y    = as.double(y[o]),
                    nnd  = as.double(nnd),
                    huge = as.double(big),
                    PACKAGE="spatstat.geom")
             nnd <- matrix(nnd, nrow=n, ncol=kmaxcalc)
             nnd[o, ] <- matrix(z$nnd, nrow=n, ncol=kmaxcalc, byrow=TRUE)
           },
           stop(paste("Unrecognised method", sQuote(method)))
           )
  }

  # post-processing
  if(kmax > kmaxcalc) {
    # add columns of Inf
    infs <- matrix(Inf, nrow=n, ncol=kmax-kmaxcalc)
    nnd <- cbind(nnd, infs)
  }

  if(kmax > 1)
    colnames(nnd) <- paste0("dist.", 1:kmax)
  
  if(length(k) < kmax) {
    # select only the specified columns
    nnd <- nnd[, k, drop=TRUE]
  }
  
  return(nnd)
}


nnwhich <- function(X, ...) {
  UseMethod("nnwhich")
}

nnwhich.ppp <- local({

  nnwhich.ppp <- function(X, ..., k=1, by=NULL,
                          proper=FALSE, method="C", metric=NULL) {
    verifyclass(X, "ppp")
    trap.extra.arguments(..., .Context="In nnwhich.ppp")
    if(is.null(by)) {
      if(!proper) {
        ## allow duplicated points
        if(is.null(metric)) {
          ## default case: Euclidean distance
          nw <- nnwhich.default(X$x, X$y, k=k, method=method)
        } else {
          ## some other distance metric
          nw <- invoke.metric(metric, "nnwhich.ppp",
                             X, ..., k=k, method=method)
        }
      } else {
        ## collapse duplicated points
        um <- uniquemap(unmark(X))
        isuniq <- (um == seq_along(um))
        UX <- X[isuniq]
        ## compute nn distances between unique points
        if(is.null(metric)) {
          nwUU <- nnwhich.default(UX$x, UX$y, k=k, method=method)
        } else {
          nwUU <- invoke.metric(metric, "nnwhich.ppp",
                                UX, ..., k=k, method=method)
        }
        ## remap to original dataset
        UtoX <- which(isuniq)
        nwUX <- nwUU
        nwUX[] <- UtoX[nwUU[]]
        nw <- marksubset(nwUX, um)
      }
    } else {
      ## split by a factor
      if(is.character(by)) {
        ## Interpret using split.ppp
        Y <- split(X, f=by, drop=FALSE)
        by <- attr(Y, "fgroup")
      }
      idX <- seq_len(npoints(X)) 
      Y <- split(X %mark% idX, f=by, un=FALSE)
      if(proper) {
        ## merge duplicates, but only within a group
        Xid <- lapply(Y, marks)
        Zid <- lapply(lapply(Y, unmark), uniquemap)
        Yid <- mapply("[", x=Xid, i=Zid, SIMPLIFY=FALSE)
        Y <- mapply("marks<-", x=Y, value=Yid, SIMPLIFY=FALSE)
        idX <- unsplit(Yid, f=by)
      }
      nwY <- lapply(Y, nnwhichsub, XX=X, iX=idX, k=k, metric=metric)
      nw <- do.call(cbind, nwY)
    }
    return(nw)
  }

  nnwhichsub <- function(Z, XX, iX, k, metric=NULL) {
    # marks(Z) gives original serial numbers of subset Z
    iY <- marks(Z)
    Zid <- nncross(XX, Z, iX=iX, iY=iY, k=k, what="which", metric=metric)
    nk <- length(k)
    if(nk == 1) {
      Yid <- iY[Zid]
    } else {
      Zid <- as.vector(as.matrix(Zid))
      Yid <- iY[Zid]
      Yid <- data.frame(which=matrix(Yid, ncol=nk))
    }
    return(Yid)
  }

  nnwhich.ppp
})


nnwhich.default <-
  function(X, Y=NULL, ..., k=1, by=NULL, method="C")
{
  warn.no.metric.support("nnwhich.default", ...)
	#  identifies nearest neighbour of each point in
	#  the pattern of points (x[i],y[i])
	#
  xy <- xy.coords(X,Y)[c("x","y")]
  x <- xy$x
  y <- xy$y

  # validate
  n <- length(x)
  if(length(y) != n)
    stop("lengths of x and y do not match")
  
  method <- match.arg(method, c("C", "interpreted", "test"))
  # other arguments ignored
  trap.extra.arguments(..., .Context="In nnwhich.default")

  # split by factor ?
  if(!is.null(by)) {
    X <- as.ppp(xy, W=boundingbox)
    return(nnwhich(X, by=by, k=k))
  }
  
  # k can be a single integer or an integer vector
  if(length(k) == 0)
    stop("k is an empty vector")
  else if(length(k) == 1) {
    if(k != round(k) || k <= 0)
      stop("k is not a positive integer")
  } else {
    if(any(k != round(k)) || any(k <= 0))
      stop(paste("some entries of the vector",
           sQuote("k"), "are not positive integers"))
  }
  k <- as.integer(k)
  kmax <- max(k)

  # special cases
  if(n <= 1) {
    # empty pattern => return integer(0)
    # or pattern with only 1 point => return NA
    nnw <- matrix(as.integer(NA), nrow=n, ncol=kmax)
    nnw <- nnw[,k, drop=TRUE]
    return(nnw)
  }

  # number of neighbours that are well-defined
  kmaxcalc <- min(n-1, kmax)

  # identify k-nn for k <= kmaxcalc

  if(kmaxcalc == 1) {
    # identify nearest neighbour only
    switch(method,
           test = ,
           interpreted={
             #  matrix of squared distances between all pairs of points
             sq <- function(a, b) { (a-b)^2 }
             squd <-  outer(x, x, sq) + outer(y, y, sq)
             #  reset diagonal to a large value so it is excluded from minimum
             diag(squd) <- Inf
             #  nearest neighbours
             nnw <- apply(squd,1,which.min)
           },
           C={
             nnw <- integer(n)
             o <- fave.order(y)
             big <- sqrt(.Machine$double.xmax)
             z <- .C(SG_nnwhichsort,
                    n = as.integer(n),
                    x = as.double(x[o]),
                    y = as.double(y[o]),
                    nnwhich = as.integer(nnw),
                    huge = as.double(big),
                    PACKAGE="spatstat.geom")
             witch <- z$nnwhich # sic 
             if(any(witch <= 0))
               stop("Internal error: non-positive index returned from C code")
             if(any(witch > n))
               stop("Internal error: index returned from C code exceeds n")
             nnw[o] <- o[witch]
           },
           stop(paste("Unrecognised method", sQuote(method)))
           )
  } else {
    # case kmaxcalc > 1
    switch(method,
           test = ,
           interpreted={
             if(n <= 1000 && method == "interpreted") {
               # form n x n matrix of squared distances
               D2 <- pairdist.default(x, y, method=method, squared=TRUE)
               # find k'th smallest squared distance
               diag(D2) <- Inf
               nnw <- t(apply(D2, 1, fave.order))[, 1:kmaxcalc]
             } else {
               # avoid creating huge matrix
               # handle one row of D at a time
               nnw <- matrix(as.integer(NA), nrow=n, ncol=kmaxcalc)
               for(i in seq_len(n)) {
                 D2i <- (x - x[i])^2 + (y - y[i])^2
                 D2i[i] <- Inf
                 nnw[i,] <- fave.order(D2i)[1:kmaxcalc]
               }      
             }
           },
           C={
             nnw <- matrix(integer(n * kmaxcalc), nrow=n, ncol=kmaxcalc)
             o <- fave.order(y)
             big <- sqrt(.Machine$double.xmax)
             z <- .C(SG_knnwhich,
                    n = as.integer(n),
                    kmax = as.integer(kmaxcalc),
                    x = as.double(x[o]),
                    y = as.double(y[o]),
                    nnwhich = as.integer(nnw),
                    huge = as.double(big),
                    PACKAGE="spatstat.geom")
             witch <- z$nnwhich # sic
             witch <- matrix(witch, nrow=n, ncol=kmaxcalc, byrow=TRUE)
             if(any(witch <= 0))
               stop("Internal error: non-positive index returned from C code")
             if(any(witch > n))
               stop("Internal error: index returned from C code exceeds n")
             # convert back to original ordering
             nnw[o,] <- matrix(o[witch], nrow=n, ncol=kmaxcalc)
           },
           stop(paste("Unrecognised method", sQuote(method)))
           )
  }
  
  # post-processing
  if(kmax > kmaxcalc) {
    # add columns of NA's
    nas <- matrix(as.numeric(NA), nrow=n, ncol=kmax-kmaxcalc)
    nnw <- cbind(nnw, nas)
  }

  if(kmax > 1)
    colnames(nnw) <- paste0("which.", 1:kmax)

  if(length(k) < kmax) {
    # select only the specified columns
    nnw <- nnw[, k, drop=TRUE]
  }
  return(nnw)
}
