#'      convexdist.R
#'
#'  Distance metric whose unit ball is a given, symmetric, convex polygon.
#'
#' $Revision: 1.21 $  $Date: 2026/01/21 06:26:39 $


convexmetric <- local({
  
  #' ..........  utilities ......................

  sptvex <- function(K, origin=c(0,0), show=FALSE) {
    ## find the support vectors of convex polygon K
    if(!missing(origin)) K <- shift(K, origin=origin)
    v <- vertices(convexhull(K))
    x <- v$x
    y <- v$y
    dx <- diff(c(x, x[1]))
    dy <- diff(c(y, y[1]))
    ll <- sqrt(dx^2+dy^2)
    co <- dy/ll
    si <- -dx/ll
    p  <- co * x + si * y
    if(any(bad <- !(is.finite(co) & is.finite(si)))) {
      ## very short segment - direction cannot be determined - use midpoint
      xmid <- x + dx/2
      ymid <- y + dy/2
      ll[bad] <- sqrt(xmid[bad]^2 + ymid[bad]^2)
      co[bad] <- ymid[bad]/ll[bad]
      si[bad] <- -xmid[bad]/ll[bad]
      p[bad] <- co[bad] * x[bad] + si[bad] * y[bad]
      if(any(verybad <- !(is.finite(co) & is.finite(si)))) {
        ## very short segment, very close to origin - remove it
        retain <- !verybad
        co <- co[retain]
        si <- si[retain]
        p  <- p[retain]
      }
    }
    if(show) {
      B <- boundingbox(Frame(K), bounding.box.xy(p*co, p*si))
      plot(B, type="n", main="")
      plot(K, add=TRUE)
      points(0,0,pch=3)
      points(p * co, p * si)
      plot(infline(p=p, theta=atan2(si,co)), lty=3)
    }
    return(data.frame(co=co, si=si, p=p))
  }

  ## find support vectors after enforcing convexity and exact symmetry

  sptvexsym <- function(K, standardise=TRUE) {
    K <- convexhull(K)
    sp <- sptvex(K, show=FALSE)
    spr <- sptvex(reflect(K), show=FALSE)
    sp <- rbind(sp, spr)
    sp <- sp[!duplicated(sp), , drop=FALSE]
    if(standardise) {
      sp$co <- sp$co/sp$p
      sp$si <- sp$si/sp$p
    }
    return(sp)
  }

  #' ..........  main 'engine' functions  ......................

  convexpairdist <- function(X, sp) {
    nX <- npoints(X)
    if(nX <= 1) return(matrix(0, nX, nX))
    a <- coords(X)
    dx <- outer(a$x, a$x, "-")
    dy <- outer(a$y, a$y, "-")
    ex <- sp$co
    ey <- sp$si
    for(i in seq_len(nrow(sp))) {
      ri <- dx * ex[i] + dy * ey[i]
      if(i == 1) r <- ri else r[] <- pmax(r, ri)
    }
    return(r)
  }

  convexnndist <- function(X, sp, k=1L) {
    d <- convexpairdist(X, sp)
    diag(d) <- Inf
    nn <- PDtoNN(d, "dist", k=k)
    return(nn)
  }

  convexnnwhich <- function(X, sp, k=1L) {
    d <- convexpairdist(X, sp)
    diag(d) <- Inf
    nw <- PDtoNN(d, "which", k=k)
    return(nw)
  }

  convexcrossdist <- function(X, Y, sp) {
    if(npoints(X) == 0 || npoints(Y) == 0)
      return(matrix(0, npoints(X), npoints(Y)))
    ex <- sp$co
    ey <- sp$si
    a <- coords(X)
    b <- coords(Y)
    dx <- outer(a$x, b$x, "-")
    dy <- outer(a$y, b$y, "-")
    for(i in seq_len(nrow(sp))) {
      ri <- dx * ex[i] + dy * ey[i]
      if(i == 1) r <- ri else r[] <- pmax(r, ri)
    }
    return(r)
  }

  convexnncross <- function(X, Y, sp, ve,
                            iX=NULL, iY=NULL, what=c("dist", "which"),
                            k=1L) {
    #' X is a point pattern
    #' Y is a point pattern or segment pattern
    what <- match.arg(what, several.ok=TRUE)
    nX <- npoints(X)
    nY <- nobjects(Y)
    if(nX == 0 || nY == 0) {
      d <- matrix(Inf, nX, nY)
    } else if(is.ppp(Y)) {
      d <- convexcrossdist(X, Y, sp)
    } else if(is.psp(Y)) {
      d <- convexPxS(X, Y, sp, ve)
    } else stop("Y should be a point pattern or line segment pattern")
    result <- XDtoNN(d, what=what, iX=iX, iY=iY, k=k)
    return(result)
  }

    convexPxS <- function(X, Y, sp, ve) {
      #' convex distance from each point of X to each segment in Y
      #' requires vertices as well as support vectors
      stopifnot(is.ppp(X))
      stopifnot(is.psp(Y))
      nX <- npoints(X)
      nY <- nsegments(Y)
      if(nX == 0 || nY == 0)
        return(matrix(, nX, nY))
      #' vertices - distance from origin
      vl <- with(ve, sqrt(x^2+y^2))
      nv <- length(vl)
      #' distances from points of X to endpoints of Y
      D1 <- convexcrossdist(X, endpoints.psp(Y, "first"), sp)
      D2 <- convexcrossdist(X, endpoints.psp(Y, "second"), sp)
      D <- matrix(pmin(D1, D2), nX, nY)
      dmax <- apply(D, 1, max)
      #' distances from points of X to locations on segments
      B <- boundingbox(Frame(X), Frame(Y))
      B <- grow.rectangle(B, max(dmax) * max(vl))
      coX <- coords(X)
      xx <- coX[, "x"]
      yy <- coX[, "y"]
      for(i in 1:nrow(coX)) {
        #' construct segments from X[i] along expansion line of each vertex
        Zi <- psp(rep(xx[i], nv),
                  rep(yy[i], nv),
                  xx[i] + dmax[i] * ve$x,
                  yy[i] + dmax[i] * ve$y,
                  window=B)
        #' intersect with target segments
        V <- crossing.psp(Zi, Y, details=TRUE)
        if(npoints(V) > 0) {
          marv <- marks(V)
          #' crossing point with which target segment?
          jj <- marv$jB
          #' crossing point of extension of which vertex?
          kk <- marv$iA
          #' Euclidean distance from X[i] to crossing point
          dE <- crossdist(X[i], V)
          #' metric distance
          dd <- dE/vl[kk]
          #' minimise over each target segment
          oo <- order(jj, dd) 
          jj <- jj[oo]
          dd <- dd[oo]
          ok <- !duplicated(jj)
          jj <- jj[ok]
          dd <- dd[ok]
          #' minimise
          D[i, jj] <- pmin(D[i, jj], dd)
        }
      }
      return(D)
    }

  convexdistmapmask <- function(w, sp, npasses=5, verbose=FALSE) {
    stopifnot(is.mask(w))
    check.1.integer(npasses)
    ## get support vectors
    sx <- sp$co
    sy <- sp$si
    ns <- length(sx)
    ## pad out mask
    nr <- w$dim[1L]
    nc <- w$dim[2L]
    xcol <- w$xcol
    yrow <- w$yrow
    #' input image will be padded out with a margin of width 2 on all sides
    mr <- mc <- 2L
    #' full dimensions of padded image
    Nnr <- nr + 2 * mr
    Nnc <- nc + 2 * mc
    N <- Nnr * Nnc
    #' output image (subset): rows & columns (R indexing)
    rmin <- mr + 1L
    rmax <- Nnr - mr
    cmin <- mc + 1L
    cmax <- Nnc - mc
    #' do padding
    x <- matrix(FALSE, nrow=Nnr, ncol=Nnc)
    x[rmin:rmax, cmin:cmax] <- w$m
    #' compute distmap
    res <- .C(SG_mdtPconv,
              as.double(xcol[1L]),
              as.double(yrow[1L]),
              as.double(xcol[nc]),
              as.double(yrow[nr]),
              nr = as.integer(nr),
              nc = as.integer(nc),
              mr = as.integer(mr),
              mc = as.integer(mc),
              inp = as.integer(t(x)),
              ns = as.integer(ns),
              sx = as.double(sx),
              sy = as.double(sy),
              npasses = as.integer(npasses),
              distances = as.double (double(N)),
              rows      = as.integer(integer(N)),
              cols      = as.integer(integer(N)),
              PACKAGE="spatstat.geom")
    dist <- matrix(res$distances,
                   ncol=Nnc, nrow=Nnr, byrow = TRUE)[rmin:rmax, cmin:cmax]
    result <- as.im(dist, w)
    edge <- TRUE
    if(edge) {
      #' calculate distance transform to boundary
      y <- x
      y[] <- TRUE
      y[rmin:rmax, cmin:cmax] <- FALSE
      y[rmin, ] <- TRUE
      y[rmax, ] <- TRUE
      y[, cmin] <- TRUE
      y[, cmax] <- TRUE
      #' compute distmap
      bres <- .C(SG_mdtPconv,
                 as.double(xcol[1L]),
                 as.double(yrow[1L]),
                 as.double(xcol[nc]),
                 as.double(yrow[nr]),
                 nr = as.integer(nr),
                 nc = as.integer(nc),
                 mr = as.integer(mr),
                 mc = as.integer(mc),
                 inp = as.integer(t(y)),
                 ns = as.integer(ns),
                 sx = as.double(sx),
                 sy = as.double(sy),
                 npasses = as.integer(npasses),
                 distances = as.double (double(N)),
                 rows      = as.integer(integer(N)),
                 cols      = as.integer(integer(N)),
                 PACKAGE="spatstat.geom")
      bdist <- matrix(bres$distances,
                      ncol=Nnc, nrow=Nnr, byrow = TRUE)[rmin:rmax, cmin:cmax]
      bdist <- as.im(bdist, w)
      attr(result, "bdist") <- bdist
      }
      return(result)
    }


    #' >>>>>>>>>>>>>> Function to create a metric <<<<<<<<<<<<<<<<<<<<<<
    
    convexmetric <- function(K) {
      stopifnot(is.owin(K))
      stopifnot(is.convex(K))
      if(!inside.owin(0, 0, K))
        stop("The origin (0,0) must be inside the set K")
      if(bdist.points(ppp(0, 0, window=K)) < sqrt(.Machine$double.eps))
        stop("The origin (0,0) must lie in the interior of K")
      ## vertices of K
      vertK <- vertices(K)
      ## support vectors of K
      spK <- sptvexsym(K, standardise=TRUE)
      if(!all(is.finite(unlist(spK))))
        stop("Support vectors are singular (infinite or undefined)",
             call.=FALSE)
      ## build object in this environment so that K, spK, vertK are accessible
      result <- list(
        pairdist.ppp=function(X, ..., squared=FALSE) {
          warn.unsupported.args(list(periodic=FALSE, method="C"), ...)
          y <- convexpairdist(X, spK)
          return(if(squared) y^2 else y)
        },
        nndist.ppp=function(X, ..., k=1L) {
          warn.unsupported.args(list(by=NULL, method="C"), ...)
          convexnndist(X, spK, k=k)
        },
        nnwhich.ppp=function(X, ..., k=1L) {
          warn.unsupported.args(list(by=NULL, method="C"), ...)
          convexnnwhich(X, spK, k=k)
        },
        crossdist.ppp=function(X, Y, ..., squared=FALSE) {
          warn.unsupported.args(list(periodic=FALSE, method="C"), ...)
          y <- convexcrossdist(X,Y,spK)
          return(if(squared) y^2 else y)
        },
        nncross.ppp=function(X,Y,iX=NULL, iY=NULL, what=c("dist","which"),
                             ..., k=1L) {
          warn.unsupported.args(list(sortby=c("range", "var", "x", "y"),
                                     is.sorted.X=FALSE, is.sorted.Y=FALSE),
                                ...)
          convexnncross(X,Y,spK,vertK, iX, iY, what, k)
        },
        distmap.ppp=function(X, ...) {
          warn.unsupported.args(list(clip=FALSE), ...)
          w <- pixellate(X, ..., preserve=TRUE)
          w <- solutionset(w > 0)
          convexdistmapmask(w, spK)
        },
        distmap.owin=function(X, ...) {
          warn.unsupported.args(list(discretise=FALSE, invert=FALSE), ...)
          w <- do.call.matched(as.mask, list(w=quote(X), ...))
          convexdistmapmask(w, spK)
        },
        distmap.psp=function(X, ...) {
          warn.unsupported.args(list(extras=TRUE, clip=FALSE), ...)
          w <- do.call.matched(psp2mask, list(x=quote(X), ...),
                               extrargs=names(formals(as.mask))[-1])
          convexdistmapmask(w, spK)
        },
        disc=function(radius=1, centre=c(0,0), ..., mask=FALSE) {
          warn.unsupported.args(list(npoly=128, delta=NULL), ...)
          check.1.real(radius)
          stopifnot(radius > 0)
          centre <- as2vector(centre)
          B <- shift(scalardilate(K, radius), vec=centre)
          if(mask) B <- as.mask(B, ...)
          return(B)
        },
        print=function(...) {
          splat("Distance metric defined by the convex set:") 
          print(K, prefix="\t")
          invisible(NULL)
        }
      )
      class(result) <- "metric"
      return(result)
    }

    
    class(convexmetric) <- "metricfun"
    attr(convexmetric, "explain") <-
      "Creates a distance metric based on a convex set K"

    convexmetric

})
  
