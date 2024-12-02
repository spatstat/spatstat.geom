##
## persp.im.R
##
##  'persp' method for image objects
##      plus annotation
##  
##  $Revision: 1.29 $ $Date: 2023/02/28 01:53:02 $
##

persp.im <- function(x, ...,
                     colmap=NULL, colin=x, apron=FALSE,
                     visible=FALSE) {
  xname <- short.deparse(substitute(x))
  xinfo <- summary(x)
  if(xinfo$type == "factor")
    stop("Perspective plot is inappropriate for factor-valued image")
  ## check whether 'col' was specified when 'colmap' was intended
  Col <- list(...)$col
  if(is.null(colmap) && !is.null(Col) && !is.matrix(Col) && length(Col) != 1)
    warning("Argument col is not a matrix. Did you mean colmap?")
  if(!missing(colin)) {
    ## separate image to determine colours
    verifyclass(colin, "im")
    if(!compatible(colin, x)) {
      ## resample 'colin' onto grid of 'x'
      colin <- as.im(colin, W=x)
    }
    if(is.null(colmap))
      colmap <- spatstat.options("image.colfun")(128)
  }
  pop <- spatstat.options("par.persp")
  ##
  if(is.function(colmap) && !inherits(colmap, "colourmap")) {
    ## coerce to a 'colourmap' if possible
    clim <- range(colin, finite=TRUE)
    if(names(formals(colmap))[1] == "n") {
      colval <- colmap(128)
      colmap <- colourmap(colval, range=clim)
    } else {
      ## colour map determined by a rule (e.g. 'beachcolours')
      colmap <- invokeColourmapRule(colmap, colin,
                                    zlim=clim, colargs=list(...))
      if(is.null(colmap))
        stop("Unrecognised syntax for colour function")
    }
  }
  ## colour map?
  if(is.null(colmap)) {
    colinfo <- list(col=NULL)
  } else if(inherits(colmap, "colourmap")) {
    ## colour map object
    ## apply colour function to image data
    colval <- eval.im(colmap(colin))
    colval <- t(as.matrix(colval))
    ## strip one row and column for input to persp.default
    colval <- colval[-1, -1]
    ## replace NA by arbitrary value
    isna <- is.na(colval)
    if(any(isna)) {
      stuff <- attr(colmap, "stuff")
      colvalues <- stuff$outputs
      colval[isna] <- colvalues[1]
    }
    ## pass colour matrix (and suppress lines)
    colinfo <- list(col=colval, border=NA)
  } else {
    ## interpret 'colmap' as colour map
    if(is.list(colmap) && all(c("breaks", "col") %in% names(colmap))) {
      breaks <- colmap$breaks
      colvalues <- colmap$col
    } else if(is.vector(colmap)) {
      colvalues <- colmap
      breaks <- quantile(colin,
                         seq(from=0,to=1,length.out=length(colvalues)+1))
      if(!all(ok <- !duplicated(breaks))) {
        breaks <- breaks[ok]
        colvalues <- colvalues[ok[-1]]
      }
    } else warning("Unrecognised format for colour map")
    ## apply colour map to image values
    colid <- cut.im(colin, breaks=breaks, include.lowest=TRUE)
    colval <- eval.im(colvalues[unclass(colid)])
    colval <- t(as.matrix(colval))
    ## strip one row and column for input to persp.default
    colval <- colval[-1, -1]
    colval[is.na(colval)] <- colvalues[1]
    ## pass colour matrix (and suppress lines)
    colinfo <- list(col=colval, border=NA)
  }

  if(apron) {
    ## add an 'apron'
    zlim <- list(...)$zlim
    bottom <- if(!is.null(zlim)) zlim[1] else min(x)
    x <- na.handle.im(x, na.replace=bottom)
    x <- padimage(x, bottom)
    xinfo <- summary(x)
    if(is.matrix(colval <- colinfo$col)) {
      colval <- matrix(col2hex(colval), nrow(colval), ncol(colval))
      grijs <- col2hex("lightgrey")
      colval <- cbind(grijs, rbind(grijs, colval, grijs), grijs)
      colinfo$col <- colval
    }
  }

  if(spatstat.options("monochrome"))
    colinfo$col <- to.grey(colinfo$col)
  
  ## get reasonable z scale while fixing x:y aspect ratio
  if(xinfo$type %in% c("integer", "real")) {
    zrange <- xinfo$range
    if(diff(zrange) > 0) {
      xbox <- as.rectangle(x)
      zscale <- 0.5 * mean(diff(xbox$xrange), diff(xbox$yrange))/diff(zrange)
      zlim <- zrange
    } else {
      zscale <- NULL
      mx <- xinfo$mean
      zlim <- mx + c(-1,1) * if(mx == 0) 0.1 else min(abs(mx), 1)
    }
  } else zscale <- zlim <- NULL

  dotargs <- list(...)
  if(spatstat.options("monochrome"))
    dotargs <- col.args.to.grey(dotargs)

  ## catch argument 'adj.main' and convert to recognised argument 'adj'
  if(!is.na(k <- match("adj.main", names(dotargs))))
    names(dotargs)[k] <- "adj"
    
  xcol <- x$xcol
  yrow <- x$yrow
  zmat <- t(x$v)
  dont.complain.about(xcol, yrow, zmat)
  yargh <- resolve.defaults(list(x=quote(xcol), y=quote(yrow), z=quote(zmat)),
                            dotargs,
                            pop,
                            colinfo,
                            list(xlab="x", ylab="y", zlab=xname),
                            list(scale=FALSE, expand=zscale,
                                 zlim=zlim),
                            list(main=xname),
                            .StripNull=TRUE)

  jawab <- do.call.matched(persp, yargh, 
                           funargs=graphicsPars("persp"))

  attr(jawab, "expand") <- yargh$expand
  
  if(visible)
    attr(jawab, "visible") <- perspVisible(x=x, M=jawab)
    
  return(invisible(jawab))
}

perspVisible <- function(x, y, z, M) {
  if(!is.matrix(M)) stop("M should be a matrix")
  if(!all(dim(M) == c(4,4))) stop("M should be a 4x4 matrix")
  ## handle all options available in persp.default
  xgiven <- !missing(x)
  ygiven <- !missing(y)
  zgiven <- !missing(z)
  if(xgiven && ygiven && zgiven) {
    xmargin <- x
    ymargin <- y
    values <- z
  } else if(xgiven && !ygiven && !zgiven) {
    values <- x
    xmargin <- ymargin <- NULL
  } else if(!xgiven && !ygiven && zgiven) {
    values <- z
    xmargin <- ymargin <- NULL
  } else stop("x or z must be given")
  ## extract data as matrix 'Xmat' and data frame 'xyz'
  if(is.im(values)) {
    ## spatstat image convention x = col, y = row
    X <- values
    Xmat <- as.matrix(X)
    xyz <- as.matrix(as.data.frame(X)) # drops NA entries
    xstep <- X$xstep
    ystep <- X$ystep
  } else if(is.matrix(x)) {
    ## base graphics image convention x = row, y = col
    ## convert to spatstat convention
    Xmat <- t(values)
    if(is.null(xmargin)) xmargin <- seq(0, 1, length.out=ncol(Xmat))
    if(is.null(ymargin)) ymargin <- seq(0, 1, length.out=nrow(Xmat))
    xyz <- cbind(x=xmargin[col(Xmat)],
                 y=ymargin[row(Xmat)],
                 z=as.numeric(Xmat))
    xyz <- xyz[complete.cases(xyz), ,drop=FALSE]
    xstep <- mean(diff(xmargin))
    ystep <- mean(diff(ymargin))
  } else stop("format is not understood")
  ## project the coordinates
  ## onto (x,y) plane of plot and z axis pointing out of it
  v <- cbind(xyz, 1) %*% M
  px <- v[,1]/v[,4]
  py <- v[,2]/v[,4]
  pz <- v[,3]/v[,4]
  pw <- v[,4]
  ## determine greatest possible difference in 'depth' in one pixel step
  PZ <- Xmat
  ok <- !is.na(PZ)
  PZ[ok] <- pz
  maxslipx <- max(0, abs(apply(PZ, 1, diff)), na.rm=TRUE)
  maxslipy <- max(0, abs(apply(PZ, 2, diff)), na.rm=TRUE)
  ## First, determine which pixels are in front
  d <- ceiling(dim(Xmat)/2)
  jx <- cut(px, breaks=d[2])
  iy <- cut(py, breaks=d[1])
  zmax <- tapply(pz, list(iy,jx), max)
  infront <- (pz > zmax[cbind(iy,jx)] - maxslipx - maxslipy)
  ## Second, determine whether outward normal to surface is pointing to viewer
  dzdx <- cbind(0, t(apply(Xmat, 1, diff)))/xstep
  dzdy <- rbind(0, apply(Xmat, 2, diff))/ystep
  dzdx <- as.vector(dzdx[ok])
  dzdy <- as.vector(dzdy[ok])
  ## unscaled normal vector
  normalx <- -dzdx
  normaly <- -dzdy
  normalz <- 1
  ## derivative of projected depth with respect to 3D input position
  dDdx <- (M[1,3] - M[1,4]/pz)/pw
  dDdy <- (M[2,3] - M[2,4]/pz)/pw
  dDdz <- (M[3,3] - M[3,4]/pz)/pw
  ## inner product = derivative of projected depth along outward normal vector
  dotprod <- normalx * dDdx + normaly * dDdy + normalz * dDdz
  ## Visible?
  isvis <- infront & (dotprod < 0)
  if(!anyNA(Xmat)) {
    answer <- isvis
  } else {
    answer <- !is.na(as.vector(Xmat))
    answer[answer] <- isvis
  }
  Vmat <- matrix(answer, nrow(Xmat), ncol(Xmat))
  if(!is.im(values)) return(t(Vmat))
  V <- (X > 0)
  V[drop=FALSE] <- Vmat
  return(V)
}

perspPoints <- function(x, y=NULL, ..., Z, M, occluded=TRUE) {
  xy <- xy.coords(x, y)
  stopifnot(is.im(Z))
  X <- as.ppp(xy, W=Frame(Z))
  if(!(is.matrix(M) && all(dim(M) == 4)))
    stop("M should be a 4 x 4 matrix, returned from persp()")
  if(occluded) {
    V <- attr(M, "visible")
    if(is.null(V)) {
      warning(paste("M does not contain visibility information;",
                    "it should be recomputed by persp() with visible=TRUE"))
    } else {
      ## restrict to visible points
      VX <- V[X, drop=FALSE]
      VX[is.na(VX)] <- FALSE
      X <- X[VX]
    }
  }
  #' determine heights
  ZX <- Z[X, drop=FALSE] # may contain NA
  #' transform and plot
  points(trans3d(X$x, X$y, ZX, M), ...)
}

perspSegments <- local({
  perspSegments <- function(x0, y0=NULL, x1=NULL, y1=NULL, ...,
                            Z, M, occluded=TRUE) {
    stopifnot(is.im(Z))
    if(!(is.matrix(M) && all(dim(M) == 4)))
      stop("M should be a 4 x 4 matrix, returned from persp()")
    if(occluded) {
      V <- attr(M, "visible")
      if(is.null(V))
        warning(paste("M does not contain visibility information;",
                      "it should be recomputed by persp() with visible=TRUE"))
    }
    
    if(is.psp(X <- x0) && is.null(y0) && is.null(x1) && is.null(y1)) {
      eX <- X$ends
#      nX <- nrow(eX)
    } else {
#      nX <- length(x0)
      check.nvector(x0, naok=TRUE, vname="x0")
      check.nvector(y0, naok=TRUE, vname="y0")
      check.nvector(x1, naok=TRUE, vname="x1")
      check.nvector(y1, naok=TRUE, vname="y1")
      eX <- cbind(x0, y0, x1, y1)
    }
    if(!occluded || is.null(V)) {
      ## no segments will be occluded
      Y <- eX
    } else {
      ## chop each segment to length of single pixel along either axis
      xstep <- Z$xstep
      ystep <- Z$ystep
      Y <- do.call(rbind, lapply(as.data.frame(t(eX)), chopsegment,
                                 eps1=xstep, eps2=ystep))
      ## determine which segments are visible
      yleft  <- list(x=Y[,1], y=Y[,2])
      yright <- list(x=Y[,3], y=Y[,4])
      ok <- V[yleft, drop=FALSE] & V[yright, drop=FALSE]
      ok[is.na(ok)] <- FALSE
      Y <- Y[ok, ,drop=FALSE]
    }
    if(nrow(Y) == 0) return(invisible(NULL))
    ## map to projected plane
    x0y0 <- trans3d(Y[,1], Y[,2], Z[list(x=Y[,1],y=Y[,2]), drop=FALSE], M)
    x1y1 <- trans3d(Y[,3], Y[,4], Z[list(x=Y[,3],y=Y[,4]), drop=FALSE], M)
    segments(x0y0$x, x0y0$y, x1y1$x, x1y1$y, ...)
  }

  chopsegment <- function(x, eps1, eps2) {
    n1 <- ceiling(abs(x[3] - x[1])/eps1)
    n2 <- ceiling(abs(x[4] - x[2])/eps2)
    n <- max(1, n1, n2)
    b <- (1:n)/n
    a <- (0:(n-1))/n
    return(cbind(x[1] + a * (x[3]-x[1]),
                 x[2] + a * (x[4]-x[2]),
                 x[1] + b * (x[3]-x[1]),
                 x[2] + b * (x[4]-x[2])))
  }
      
  perspSegments
})

perspLines <- function(x, y=NULL, ..., Z, M, occluded=TRUE) {
  xy <- xy.coords(x, y)
  n <- length(xy$x)
  perspSegments(x[-n], y[-n], x[-1], y[-1], Z=Z, M=M, ..., occluded=occluded)
}

perspContour <- function(Z, M, ...,
                         nlevels=10, levels=pretty(range(Z), nlevels),
                         occluded=TRUE) {
  cl <- contourLines(x=Z$xcol,
                     y=Z$yrow,
                     z=t(Z$v),
                     nlevels=nlevels, levels=levels)
  for(i in seq_along(cl)) {
    cli <- cl[[i]]
    perspLines(cli$x, cli$y, ..., Z=Z, M=M, occluded=occluded)
  }
  invisible(NULL)
}

