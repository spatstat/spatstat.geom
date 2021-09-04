#'      convexdist.R
#'
#'  Distance metric whose unit ball is a given, symmetric, convex polygon.
#'
#' $Revision: 1.7 $  $Date: 2021/09/04 10:38:39 $


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

  #' ..........  main functions  ......................

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

  convexnndist <- function(X, sp) {
    nX <- npoints(X)
    if(nX <= 1) return(numeric(0))
    if(nX == 1) return(Inf)
    d <- convexpairdist(X, sp)
    diag(d) <- Inf
    nn <- apply(d, 1, min)
    return(nn)
  }

  convexnnwhich <- function(X, sp) {
    nX <- npoints(X)
    if(nX == 0) return(integer(0))
    if(nX == 1) return(NA_integer_)
    d <- convexpairdist(X, sp)
    diag(d) <- Inf
    nw <- apply(d, 1, which.min)
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

  convexnncross <- function(X, Y, sp, what=c("dist", "which")) {
    what <- match.arg(what, several.ok=TRUE)
    nX <- npoints(X)
    nY <- npoints(Y)
    ## trivial cases
    if(nX == 0 || nY == 0) {
      nnd <- rep(Inf, nX)
      nnw <- rep(NA_integer_, nX)
    } else {
      d <- convexcrossdist(X, Y, sp)
      nnd <- if("dist" %in% what) apply(d, 1, min) else NULL
      nnw <- if("which" %in% what) apply(d, 1, which.min) else NULL
    }
    result <- list(dist=nnd, which=nnw)[what]
    result <- as.data.frame(result)[,,drop=TRUE]
    return(result)
  }

  
  convexdistmap <- function(X, sp, npasses=5, verbose=FALSE) {
    ## converts X to a set
    w <- as.mask(X)
    check.1.integer(npasses)
    ## get support vectors
    sx <- sp$co
    sy <- sp$si
    ns <- length(sx)
    ## pad out mask
    nr <- w$dim[1L]
    nc <- w$dim[2L]
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
              as.double(w$xrange[1L]),
              as.double(w$yrange[1L]),
              as.double(w$xrange[2L]),
              as.double(w$yrange[2L]),
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
                 as.double(w$xrange[1L]),
                 as.double(w$yrange[1L]),
                 as.double(w$xrange[2L]),
                 as.double(w$yrange[2L]),
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

    convexmetric <- function(K) {
      stopifnot(is.owin(K))
      stopifnot(is.convex(K))
      if(!inside.owin(0, 0, K))
        stop("The origin (0,0) must be inside the set K")
      if(bdist.points(ppp(0, 0, window=K)) < sqrt(.Machine$double.eps))
        stop("The origin (0,0) must lie in the interior of K")
      spK <- sptvexsym(K, standardise=TRUE)
      if(!all(is.finite(unlist(spK))))
        stop("Support vectors are singular (infinite or undefined)",
             call.=FALSE)
      result <- list(
        tasks=list(
          pairdist.ppp="pairdist",
          nndist.ppp="nndist",
          nnwhich.ppp="nnwhich",
          crossdist.ppp="crossdist",
          nncross.ppp="nncross",
          distmap.ppp="distmap",
          distmap.owin="distmap",
          distmap.psp="distmap",
          disc="ball",
          print="print"
        ),
        functions=list(
          pairdist=function(X) { convexpairdist(X, spK) },
          nndist=function(X) { convexnndist(X, spK) },
          nnwhich=function(X) { convexnnwhich(X, spK) },
          crossdist=function(X, Y) { convexcrossdist(X,Y,spK) },
          nncross=function(X,Y,what=c("dist","which")) {
            convexnncross(X,Y,spK,what)
          },
          distmap=function(X, ...) {
            ## handles all spatial objects by discretising
            X <- if(is.psp(X)) as.mask.psp(X, ...) else as.mask(X, ...)
            convexdistmap(X, spK)
          },
          ball=function(radius=1, centre=c(0,0), ...) {
            check.1.real(radius)
            stopifnot(radius > 0)
            centre <- as2vector(centre)
            B <- shift(scalardilate(K, radius), vec=centre)
            return(B)
          },
          print=function(...) {
            splat("Distance metric defined by the convex set:") 
            print(K)
            invisible(NULL)
          })
        )
      class(result) <- "metric"
      return(result)
    }

    class(convexmetric) <- "metricfun"
    attr(convexmetric, "explain") <-
      "Creates a distance metric based on a convex set K"

    convexmetric

})
  
