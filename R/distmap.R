#
#
#      distmap.R
#
#      $Revision: 1.34 $     $Date: 2023/08/28 06:38:54 $
#
#
#     Distance transforms
#
#
distmap <- function(X, ...) {
  UseMethod("distmap")
}

distmap.ppp <- function(X, ..., clip=FALSE, metric=NULL) {
  verifyclass(X, "ppp")
  if(!is.null(metric)) {
    ans <- invoke.metric(metric, "distmap.ppp", X=X, ..., clip=clip)
    return(ans)
  }
  e <- exactdt(X, ...)
  W <- e$w
  uni <- unitname(W)
  dmat <- e$d
  imat <- e$i
  V <- im(dmat, W$xcol, W$yrow, unitname=uni)
  I <- im(imat, W$xcol, W$yrow, unitname=uni)
  if(X$window$type == "rectangle") {
    # distance to frame boundary
    bmat <- e$b
    B <- im(bmat, W$xcol, W$yrow, unitname=uni)
  } else {
    ## distance to window boundary, not frame boundary
    bmat <- bdist.pixels(W, style="matrix")
    B <- im(bmat, W$xcol, W$yrow, unitname=uni)
    if(clip) {
      ## clip all to window
      V <- V[W, drop=FALSE]
      I <- I[W, drop=FALSE]
      B <- B[W, drop=FALSE]
    }
  }
  attr(V, "index") <- I
  attr(V, "bdry")  <- B
  return(V)
}

distmap.owin <- function(X, ..., discretise=FALSE, invert=FALSE,
                         connect=8, metric=NULL) {
  verifyclass(X, "owin")
  uni <- unitname(X)
  if(!is.null(metric)) {
    ans <- invoke.metric(metric, "distmap.owin", X=X, ...,
                         discretise=discretise, invert=invert)
    return(ans)
  }
  if(is.empty(X)) {
    ## handle empty window
    Dist <- as.im(Inf, X)
    attr(Dist, "bdry") <- framedist.pixels(X, ...)
    return(Dist)
  }
  if(X$type == "rectangle") {
    M <- as.mask(X, ...)
    Bdry <- im(bdist.pixels(M, style="matrix"),
               M$xcol, M$yrow, unitname=uni)
    if(!invert)
      Dist <- as.im(M, value=0)
    else 
      Dist <- Bdry
  } else if(X$type == "polygonal" && !discretise) {
    Edges <- edges(X)
    Dist <- distmap(Edges, ...)
    Bdry <- attr(Dist, "bdry")
    if(!invert) 
      Dist[X] <- 0
    else {
      bb <- as.rectangle(X)
      bigbox <- grow.rectangle(bb, diameter(bb)/4)
      Dist[complement.owin(X, bigbox)] <- 0
    }
  } else {
    check.1.integer(connect)
    if(!(connect %in% c(8, 24)))
      stop("Argument 'connect' must equal 8 or 24", call.=FALSE)
    X <- as.mask(X, ...)
    if(invert)
      X <- complement.owin(X)
    xc <- X$xcol
    yr <- X$yrow
    nr <- X$dim[1L]
    nc <- X$dim[2L]
    ## pad out the input image with a margin of width 1 on all sides
    mat <- X$m
    mat <- cbind(FALSE, mat, FALSE)
    mat <- rbind(FALSE, mat, FALSE)
    ## call C routine
    res <- .C(SG_distmapbin,
              connect=as.integer(connect),
              xmin=as.double(xc[1L]),
              ymin=as.double(yr[1L]),
              xmax=as.double(xc[nc]),
              ymax=as.double(yr[nr]),
              nr = as.integer(nr),
              nc = as.integer(nc),
              inp = as.integer(as.logical(t(mat))),
              distances = as.double(matrix(0, ncol = nc + 2, nrow = nr + 2)),
              boundary = as.double(matrix(0, ncol = nc + 2, nrow = nr + 2)),
              PACKAGE="spatstat.geom")
  # strip off margins again
    dist <- matrix(res$distances,
                   ncol = nc + 2, byrow = TRUE)[2:(nr + 1), 2:(nc +1)]
    bdist <- matrix(res$boundary,
                    ncol = nc + 2, byrow = TRUE)[2:(nr + 1), 2:(nc +1)]
  # cast as image objects
    Dist <- im(dist,  xc, yr, unitname=uni)
    Bdry <- im(bdist, xc, yr, unitname=uni)
  }
  attr(Dist, "bdry")  <- Bdry
  return(Dist)
}

distmap.psp <- function(X, ..., extras=TRUE, clip=FALSE, metric=NULL) {
  verifyclass(X, "psp")
  if(!is.null(metric)) {
    ans <- invoke.metric(metric, "distmap.psp", X=X, ...,
                         extras=extras, clip=clip)
    return(ans)
  }
  W <- Window(X)
  uni <- unitname(W)
  M <- as.mask(W, ...)
  ## handle empty pattern
  if(nsegments(X) == 0) {
    Dist <- as.im(Inf, W)
    if(extras) {
      Indx <- as.im(NA, W)
      Bdry <- bdist.pixels(M)
      if(clip) {
       Indx <- Indx[M, drop=FALSE]
       Bdry <- Bdry[M, drop=FALSE]
      }
      attr(Dist, "index") <- Indx
      attr(Dist, "bdry")  <- Bdry
    }
    return(Dist)
  }
  rxy <- rasterxy.mask(M)
  xp <- rxy$x
  yp <- rxy$y
  E <- X$ends
  big <- 2 * diameter(Frame(W))^2
  z <- NNdist2segments(xp, yp, E$x0, E$y0, E$x1, E$y1, big, wantindex=extras)
  xc <- M$xcol
  yr <- M$yrow
  Dist <- im(array(sqrt(z$dist2), dim=M$dim), xc, yr, unitname=uni)
  if(clip <- clip && !is.rectangle(W))
    Dist <- Dist[M, drop=FALSE]
  if(extras) {
    Indx <- im(array(z$index, dim=M$dim), xc, yr, unitname=uni)
    Bdry <- im(bdist.pixels(M, style="matrix"), xc, yr, unitname=uni)
    if(clip) {
       Indx <- Indx[M, drop=FALSE]
       Bdry <- Bdry[M, drop=FALSE]
    }
    attr(Dist, "index") <- Indx
    attr(Dist, "bdry")  <- Bdry
  }
  return(Dist)
}

