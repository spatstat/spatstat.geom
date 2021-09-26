#'
#'  nnutils.R
#'
#'  Utilities for extracting nndist/nncross from distance matrices
#'
#'  $Revision: 1.3 $  $Date: 2021/09/26 07:32:16 $


PDtoNN <- function(d, what=c("dist", "which"), k=1L, ...) {
  ## Given a matrix of pairwise distances,
  ## determine the nearest neighbours
  ## and return in standard format
  stopifnot(is.matrix(d))
  stopifnot(nrow(d) == ncol(d))
  nX <- nrow(d)

  what   <- match.arg(what, several.ok=TRUE)
  want.dist  <- "dist" %in% what 
  want.which <- "which" %in% what
  want.both  <- want.dist && want.which

  if(!missing(k)) {
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
  }
  k <- as.integer(k)
  kmax <- max(k)
  nk <- length(k)
  kmaxcalc <- min(nX, kmax) # number of neighbours that are well-defined

  ## deal with null cases
  if(nX == 0) 
    return(as.data.frame(list(dist=matrix(0, nrow=0, ncol=nk),
                              which=matrix(0L, nrow=0, ncol=nk)))[,what])

  ##
  diag(d) <- Inf
  NND <- NNW <- NULL
  if(kmax == 1L) {
    if(want.dist)  NND <- apply(d, 1, min) 
    if(want.which) NNW <- apply(d, 1, which.min) 
  } else {
    kuse <- k[k <= kmaxcalc]
    nkuse <- length(kuse)
    kmap <- match(kuse, k)
    if(want.dist) {
      NND <- apply(d, 1, orderstats, k=kuse)
      if(nkuse > 1) NND <- t(NND)
      if(nk > nkuse) {
        NNDfull <- matrix(Inf, nrow=nX, ncol=nk)
        NNDfull[, kmap] <- NND
        NND <- NNDfull
      }
    }
    if(want.which) {
      NNW <- apply(d, 1, orderwhich, k=kuse)
      if(nkuse > 1) NNW <- t(NNW)
      if(nk > nkuse) {
        NNWfull <- matrix(NA_integer_, nrow=nX, ncol=nk)
        NNWfull[, kmap] <- NNW
        NNW <- NNWfull
      }
    }
  }
  result <- packupNNdata(NND, NNW, what, k)
  return(result)
}

XDtoNN <- function(d, what=c("dist", "which"),
                   iX=NULL, iY=NULL, k=1L, ...) {
  ## Given a matrix of cross-pairwise distances,
  ## determine the nearest neighbours
  ## and return in standard format
  stopifnot(is.matrix(d))
  nX <- nrow(d)
  nY <- ncol(d)

  what   <- match.arg(what, several.ok=TRUE)
  want.dist  <- "dist" %in% what 
  want.which <- "which" %in% what
  want.both  <- want.dist && want.which

  if(!missing(k)) {
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
  }
  k <- as.integer(k)
  kmax <- max(k)
  nk <- length(k)
  kmaxcalc <- min(nY, kmax) # number of neighbours that are well-defined

  ## deal with null cases
  if(nX == 0) 
    return(as.data.frame(list(dist=matrix(0, nrow=0, ncol=nk),
                              which=matrix(0L, nrow=0, ncol=nk)))[,what])
  if(nY == 0)
    return(as.data.frame(list(dist=matrix(Inf, nrow=nX, ncol=nk),
                              which=matrix(NA_integer_, nrow=nX, ncol=nk))[what]))

  ## exclusion of identical pairs
  if(is.null(iX) != is.null(iY))
    stop("If one of iX, iY is given, then both must be given")
  exclude <- (!is.null(iX) || !is.null(iY))
  if(exclude) {
    stopifnot(is.integer(iX) && is.integer(iY))
    if(length(iX) != nX)
      stop("length of iX does not match the number of points in X")
    if(length(iY) != nY)
      stop("length of iY does not match the number of points in Y")
    d[cbind(iX, iY)] <- Inf
    need.dist <- want.which
  }

  ##
  NND <- NNW <- NULL
  if(kmax == 1L) {
    if(want.dist || need.dist)  NND <- apply(d, 1, min) 
    if(want.which) NNW <- apply(d, 1, which.min) 
  } else {
    kuse <- k[k <= kmaxcalc]
    nkuse <- length(kuse)
    kmap <- match(kuse, k)
    if(want.dist || need.dist) {
      NND <- apply(d, 1, orderstats, k=kuse)
      if(nkuse > 1) NND <- t(NND)
      if(nk > nkuse) {
        NNDfull <- matrix(Inf, nrow=nX, ncol=nk)
        NNDfull[, kmap] <- NND
        NND <- NNDfull
      }
    }
    if(want.which) {
      NNW <- apply(d, 1, orderwhich, k=kuse)
      if(nkuse > 1) NNW <- t(NNW)
      if(nk > nkuse) {
        NNWfull <- matrix(NA_integer_, nrow=nX, ncol=nk)
        NNWfull[, kmap] <- NNW
        NNW <- NNWfull
      }
    }
  }
  ## 
  if(want.which && exclude) {
    if(any(nope <- is.infinite(NND)))
      NNW[nope] <- NA
  }
  ## 
  result <- packupNNdata(NND, NNW, what, k)
  return(result)
}


packupNNdata <- function(NND, NNW, what, k) {
  result <- as.data.frame(list(dist=NND, which=NNW)[what])
  colnames(result) <-
    if(max(k) == 1L) {
      c(if("dist" %in% what) "dist" else NULL,
        if("which" %in% what) "which" else NULL)
    } else {
      c(if("dist" %in% what) paste0("dist.", k) else NULL,
        if("which" %in% what) paste0("which.",k) else NULL)
    }
  if(ncol(result) == 1L)
    result <- result[, , drop=TRUE]
  return(result)
}
