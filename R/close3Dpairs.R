#
# close3Dpairs.R
#
#   $Revision: 1.18 $   $Date: 2022/06/15 01:29:04 $
#
#  extract the r-close pairs from a 3D dataset
# 
#
closepairs.pp3 <- local({

  closepairs.pp3 <- function(X, rmax, twice=TRUE,
                             what=c("all", "indices", "ijd"),
                             distinct=TRUE, neat=TRUE, ...) {
    verifyclass(X, "pp3")
    what <- match.arg(what)
    stopifnot(is.numeric(rmax) && length(rmax) == 1L)
    stopifnot(is.finite(rmax))
    stopifnot(rmax >= 0)
    ordered <- list(...)$ordered
    if(missing(twice) && !is.null(ordered)) {
      warning("Obsolete argument 'ordered' has been replaced by 'twice'")
      twice <- ordered
    }
    npts <- npoints(X)
    nama <- switch(what,
                   all = c("i", "j",
                           "xi", "yi", "zi",
                           "xj", "yj", "zj",
                           "dx", "dy", "dz",
                           "d"),
                   indices = c("i", "j"),
                   ijd     = c("i", "j", "d"))
    names(nama) <- nama
    if(npts == 0) {
      null.answer <- lapply(nama, nuttink)
      return(null.answer)
    }
    ## sort points by increasing x coordinate
    oo <- fave.order(coords(X)$x)
    Xsort <- X[oo]
    ## First make an OVERESTIMATE of the number of pairs
    nsize <- list(...)$nsize # secret option to test overflow code
    if(!is.null(nsize)) {
      splat("Using nsize =", nsize)
    } else {
      #' normal usage
      npairs <- as.double(npts)^2
      if(npairs <= 1024) {
        nsize <- 1024
      } else {
        catchfraction <- (4/3) * pi * (rmax^3)/volume(as.box3(X))
        nsize <- ceiling(4 * catchfraction * npairs)
        nsize <- min(nsize, npairs)
        nsize <- max(1024, nsize)
        if(nsize > .Machine$integer.max) {
          warning(
            "Estimated number of close pairs exceeds maximum possible integer",
            call.=FALSE)
          nsize <- .Machine$integer.max
        }
      }
    }
    ## Now extract pairs
    XsortC <- coords(Xsort)
    x <- XsortC$x
    y <- XsortC$y
    z <- XsortC$z
    r <- rmax
    ng <- nsize
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    storage.mode(z) <- "double"
    storage.mode(r) <- "double"
    storage.mode(ng) <- "integer"
    ## go
    a <- switch(what,
                all = {
                  .Call(SG_close3pairs,
                        xx=x, yy=y, zz=z, rr=r, nguess=ng,
                        PACKAGE="spatstat.geom")
                },
                  indices = {
                  .Call(SG_close3IJpairs,
                        xx=x, yy=y, zz=z, rr=r, nguess=ng,
                        PACKAGE="spatstat.geom")
                },
                    ijd = {
                  .Call(SG_close3IJDpairs,
                        xx=x, yy=y, zz=z, rr=r, nguess=ng,
                        PACKAGE="spatstat.geom")
                })
    names(a) <- nama
    ## convert i,j indices to original sequence
    a$i <- oo[a$i]
    a$j <- oo[a$j]
    ## handle options
    if(twice) {
      ## both (i, j) and (j, i) should be returned
      a <- as.data.frame(a)
      a <- as.list(rbind(a, swapdata(a, what)))
    } else if(neat) {
      ## enforce i < j
      swap <- with(a, (j < i))
      if(any(swap)) {
        a <- as.data.frame(a)
        a[swap,] <- swapdata(a[swap, ,drop=FALSE], what)
        a <- as.list(a)
      }
    }
    ## add pairs of identical points?
    if(!distinct) {
      ii <- seq_len(npts)
      xtra <- switch(what,
                     indices = {
                       data.frame(i = ii, j=ii)
                     },
                     ijd= {
                       data.frame(i = ii, j=ii, d=0)
                     },
                     all = {
                       cooi <- cooj <- coords(X)[, c("x","y","z")]
                       names(cooi) <- c("xi", "yi", "zi")
                       names(cooj) <- c("xj", "yj", "zj")
                       zero <- numeric(npts)
                       cbind(data.frame(i=ii, j=ii),
                             cooi,
                             cooj,
                             data.frame(dx=zero, dy=zero, dz=zero, d=zero))
                     })
      a <- as.list(rbind(as.data.frame(a), xtra))
    }
    ## done
    return(a)
  }

  swapdata <- function(a, what) {
    switch(what,
           all = {
             with(a, data.frame(i  =  j,
                                j  =  i,
                                xi =  xj,
                                yi =  yj,
                                zi =  zj,
                                xj =  xi,
                                yj =  yi,
                                zj =  zi,
                                dx = -dx,
                                dy = -dy,
                                dz = -dz,
                                d  =  d))
           },
           indices = {
             with(a, data.frame(i=j,
                                j=i))
           },
           ijd = {
             with(a, data.frame(i=j,
                                j=i,
                                d=d))
           })
  }
  
  nuttink <- function(x) numeric(0)

  closepairs.pp3
})

#######################

crosspairs.pp3 <- local({

  crosspairs.pp3 <- function(X, Y, rmax, what=c("all", "indices", "ijd"), ...) {
    verifyclass(X, "pp3")
    verifyclass(Y, "pp3")
    what <- match.arg(what)
    stopifnot(is.numeric(rmax) && length(rmax) == 1L && rmax >= 0)
    nama <- switch(what,
                   all = c("i", "j",
                           "xi", "yi", "zi",
                           "xj", "yj", "zj",
                           "dx", "dy", "dz",
                           "d"),
                   indices = c("i", "j"),
                   ijd = c("i", "j", "d"))
    names(nama) <- nama
    nX <- npoints(X)
    nY <- npoints(Y)
    if(nX == 0 || nY == 0) {
      null.answer <- lapply(nama, nuttink)
      return(null.answer)
    }
    ## order patterns by increasing x coordinate
    ooX <- fave.order(coords(X)$x)
    Xsort <- X[ooX]
    ooY <- fave.order(coords(Y)$x)
    Ysort <- Y[ooY]
    ## First (over)estimate the number of pairs
    nsize <- list(...)$nsize # secret option to test overflow code
    if(!is.null(nsize)) {
      splat("Using nsize =", nsize)
    } else {
      #' normal usage
      nXY <- as.double(nX) * as.double(nY)
      if(nXY <= 1024) {
        nsize <- 1024
      } else {
        catchfraction <- (4/3) * pi * (rmax^3)/volume(as.box3(Y))
        nsize <- ceiling(4 * catchfraction * nXY)
        nsize <- min(nXY, nsize)
        nsize <- max(1024, nsize)
        if(nsize > .Machine$integer.max) {
          warning(
            "Estimated number of close pairs exceeds maximum possible integer",
            call.=FALSE)
          nsize <- .Machine$integer.max
        }
      }
    }
    ## .Call
    XsortC <- coords(Xsort)
    YsortC <- coords(Ysort)
    Xx <- XsortC$x
    Xy <- XsortC$y
    Xz <- XsortC$z
    Yx <- YsortC$x
    Yy <- YsortC$y
    Yz <- YsortC$z
    r <- rmax
    ng <- nsize
    storage.mode(Xx) <- storage.mode(Xy) <- storage.mode(Xz) <- "double"
    storage.mode(Yx) <- storage.mode(Yy) <- storage.mode(Yz) <- "double"
    storage.mode(r) <- "double"
    storage.mode(ng) <- "integer"
    ## go
    a <- switch(what,
                all = {
                  .Call(SG_cross3pairs,
                        xx1=Xx, yy1=Xy, zz1=Xz,
                        xx2=Yx, yy2=Yy, zz2=Yz,
                        rr=r, nguess=ng,
                        PACKAGE="spatstat.geom")
                },
                  indices = {
                  .Call(SG_cross3IJpairs,
                        xx1=Xx, yy1=Xy, zz1=Xz,
                        xx2=Yx, yy2=Yy, zz2=Yz,
                        rr=r, nguess=ng,
                        PACKAGE="spatstat.geom")
                },
                    ijd = {
                  .Call(SG_cross3IJDpairs,
                        xx1=Xx, yy1=Xy, zz1=Xz,
                        xx2=Yx, yy2=Yy, zz2=Yz,
                        rr=r, nguess=ng,
                        PACKAGE="spatstat.geom")
                })
    names(a) <- nama
    ## convert i,j indices to original sequence
    a$i <- ooX[a$i]
    a$j <- ooY[a$j]
    return(a)
  }

  nuttink <- function(x) numeric(0)
  
  crosspairs.pp3
})


