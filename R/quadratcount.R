#
#  quadratcount.R
#
#  $Revision: 1.68 $  $Date: 2025/04/05 06:28:43 $
#

quadratcount <- function(X, ...) {
  UseMethod("quadratcount")
}

quadratcount.splitppp <- function(X, ...) {
  solapply(X, quadratcount, ...)
}

quadratcount.ppp <- function(X, nx=5, ny=nx, ...,
                             xbreaks=NULL, ybreaks=NULL, left.open=TRUE, 
                             tess=NULL)  {
  verifyclass(X, "ppp")
  W <- X$window

  if(is.null(tess)) {
    ## rectangular boundaries 
    if(!is.numeric(nx))
      stop("nx should be numeric")
    ## start with rectangular tessellation
    tess <- quadrats(as.rectangle(W),
                     nx=nx, ny=ny, xbreaks=xbreaks, ybreaks=ybreaks)
    ## fast code for counting points in rectangular grid
    Xcount <- rectquadrat.countEngine(X$x, X$y, tess$xgrid, tess$ygrid,
                                      left.open=left.open)
    ## 
    if(W$type != "rectangle") {
      # intersections of rectangles with window including empty intersections
      tess <- quadrats(X,
                       nx=nx, ny=ny, xbreaks=xbreaks, ybreaks=ybreaks,
                       keepempty=TRUE)
      nonempty <- !tiles.empty(tess)
      if(!any(nonempty))
        stop("All tiles are empty")
      if(!all(nonempty)) {
        ## reshape counts as 1-dim table corresponding to tiles
        Xcount <- as.integer(t(Xcount))
        ## retain nonempty tiles
        tess   <- tess[nonempty]
        ## retain corresponding counts 
        Xcount <- Xcount[nonempty]
        ## attach tile names to counts
        Xcount <- array(Xcount,
                        dimnames=list(tile=tilenames(tess)))
        class(Xcount) <- "table"
      }
    }
  } else {
    # user-supplied tessellation
    if(!inherits(tess, "tess")) {
      tess <- try(as.tess(tess), silent=TRUE)
      if(inherits(tess, "try-error"))
        stop("The argument tess should be a tessellation", call.=FALSE)
    }
    if(tess$type == "rect") {
      # fast code for counting points in rectangular grid
      Xcount <- rectquadrat.countEngine(X$x, X$y, tess$xgrid, tess$ygrid,
                                        left.open=left.open)
    } else {
      # quadrats are another type of tessellation
      Y <- cut(X, tess)
      if(anyNA(marks(Y)))
        warning("Tessellation does not contain all the points of X")
      Xcount <- table(tile=marks(Y))
    }
  }
  attr(Xcount, "tess") <- tess
  class(Xcount) <- c("quadratcount", class(Xcount))
  return(Xcount)
}

plot.quadratcount <- function(x, ...,
                              add=FALSE, entries=as.integer(t(x)),
                              dx=0, dy=0, show.tiles=TRUE,
                              textargs = list()) {
  xname <- short.deparse(substitute(x))
  tess <- attr(x, "tess")
  # add=FALSE, show.tiles=TRUE  => plot tiles + numbers
  # add=FALSE, show.tiles=FALSE => plot window (add=FALSE) + numbers
  # add=TRUE,  show.tiles=TRUE  => plot tiles  (add=TRUE) + numbers
  # add=TRUE,  show.tiles=FALSE => plot numbers
  if(show.tiles || !add) {
    context <- if(show.tiles) tess else as.owin(tess)
    dont.complain.about(context)
    do.call(plot,
            resolve.defaults(list(quote(context), add=add),
                             list(...),
                             list(main=xname),
                             .StripNull=TRUE))
  }
  if(!is.null(entries)) {
    labels <- paste(as.vector(entries))
    til <- tiles(tess)
    incircles <- lapply(til, incircle)
    x0 <- sapply(incircles, getElement, name="x")
    y0 <- sapply(incircles, getElement, name="y")
    ra <- sapply(incircles, getElement, name="r")
    xx <- x0 + dx *ra
    yy <- y0 + dy *ra
    dont.complain.about(xx, yy, labels)
    do.call.matched(text.default,
                    resolve.defaults(list(x=quote(xx),
                                          y = quote(yy),
                                          labels=quote(labels)),
                                     textargs, 
                                     list(...)),
                    funargs=graphicsPars("text"))
  }
  return(invisible(NULL))
}

rectquadrat.breaks <- function(xr, yr, nx=5, ny=nx, xbreaks=NULL, ybreaks=NULL) {
  if(is.null(xbreaks))
    xbreaks <- seq(from=xr[1], to=xr[2], length.out=nx+1)
  else if(min(xbreaks) > xr[1] || max(xbreaks) < xr[2])
    stop("xbreaks do not span the range of x coordinates in the window")
  if(is.null(ybreaks))
    ybreaks <- seq(from=yr[1], to=yr[2], length.out=ny+1)
  else if(min(ybreaks) > yr[1] || max(ybreaks) < yr[2])
    stop("ybreaks do not span the range of y coordinates in the window")
  return(list(xbreaks=xbreaks, ybreaks=ybreaks))
}

rectquadrat.countEngine <- function(x, y, xbreaks, ybreaks, weights,
                                    left.open=TRUE) {
  if(length(x) > 0) {
    # check validity of breaks
    if(!all(inside.range(range(x), range(xbreaks))))
      stop("xbreaks do not span the actual range of x coordinates in data")
    if(!all(inside.range(range(y), range(ybreaks))))
      stop("ybreaks do not span the actual range of y coordinates in data")
  }
  # WAS: 
  # xg <- cut(x, breaks=xbreaks, include.lowest=TRUE)
  # yg <- cut(y, breaks=ybreaks, include.lowest=TRUE)
  xg <- fastFindInterval(x, xbreaks, labels=TRUE, left.open=left.open)
  yg <- fastFindInterval(y, ybreaks, labels=TRUE, left.open=left.open)
  if(missing(weights)) {
    sumz <- table(list(y=yg, x=xg))
  } else {
    # was: 
    # sumz <- tapply(weights, list(y=yg, x=xg), sum)
    # if(any(nbg <- is.na(sumz)))
    #  sumz[nbg] <- 0
    sumz <- tapplysum(weights, list(y=yg, x=xg), do.names=TRUE)
  }
  # reverse order of y 
  sumz <- sumz[rev(seq_len(nrow(sumz))), ]
  sumz <- as.table(sumz)
  #
  attr(sumz, "xbreaks") <- xbreaks
  attr(sumz, "ybreaks") <- ybreaks
  return(sumz)
}

quadrats <- function(X, nx=5, ny=nx, xbreaks = NULL, ybreaks = NULL,
                     keepempty=FALSE) {
  W <- as.owin(X)
  xr <- W$xrange
  yr <- W$yrange
  b <- rectquadrat.breaks(xr, yr, nx, ny, xbreaks, ybreaks)
  # rectangular tiles
  Z <- tess(xgrid=b$xbreaks, ygrid=b$ybreaks, unitname=unitname(W))
  if(W$type != "rectangle") {
    # intersect rectangular tiles with window W
    if(!keepempty) {
      Z <- intersect.tess(Z, W)
    } else {
      til <- tiles(Z)
      for(i in seq_along(til))
        til[[i]] <- intersect.owin(til[[i]], W)
      Z <- tess(tiles=til, window=W, keepempty=TRUE)
    }
  }
  return(Z)
}

as.tess.quadratcount <- function(X) {
  Y <- attr(X, "tess")
  m <- as.integer(t(X)) ## counts in order corresponding to tiles
  marks(Y) <- m
  return(Y)
}

as.owin.quadratcount <- function(W, ..., fatal=TRUE) {
  return(as.owin(as.tess(W), ..., fatal=fatal))
}

domain.quadratcount <- Window.quadratcount <- function(X, ...) { as.owin(X) }

intensity.quadratcount <- function(X, ..., image=FALSE) {
  Y <- as.tess(X)  # marks are counts
  lambda <- marks(Y)[,1]/tile.areas(Y) 
  if(image) {
    ## make an image
    tileid <- as.im(Y, ...)  # values are tile index
    result <- eval.im(lambda[tileid])
  } else {
    ## save as a table corresponding to X
    result <- X
    flip <- (length(dim(X)) == 2)
    if(flip) result <- t(result)
    result[] <- lambda
    if(flip) result <- t(result)
    class(result) <- "table"
    attr(result, "tess") <- NULL
  }
  return(result)
}

## The shift method is undocumented.
## It is only needed in plot.listof / plot.solist / plot.layered

shift.quadratcount <- function(X, ...) {
  attr(X, "tess") <- te <- shift(attr(X, "tess"), ...)
  attr(X, "lastshift") <- getlastshift(te)
  return(X)
}

