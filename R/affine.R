#
#	affine.R
#
#	$Revision: 1.58 $	$Date: 2024/04/26 02:29:54 $
#

affinexy <- function(X, mat=diag(c(1,1)), vec=c(0,0), invert=FALSE) {
  if(length(X$x) == 0 && length(X$y) == 0)
    return(list(x=numeric(0),y=numeric(0)))
  if(invert) {
    mat <- invmat <- solve(mat)
    vec <- - as.numeric(invmat %*% vec)
  }
  # Y = M X + V
  ans <- mat %*% rbind(X$x, X$y) + matrix(vec, nrow=2L, ncol=length(X$x))
  return(list(x = ans[1L,],
              y = ans[2L,]))
}

affinexypolygon <- function(p, mat=diag(c(1,1)), vec=c(0,0),
                             detmat=det(mat)) {
  # transform (x,y)
  p[c("x","y")] <- affinexy(p, mat=mat, vec=vec)
  # transform area
  if(!is.null(p$area))
    p$area <- p$area * detmat
  # if map has negative sign, cyclic order was reversed; correct it
  if(detmat < 0)
    p <- reverse.xypolygon(p, adjust=TRUE)
  return(p)
}
       
"affine" <- function(X, ...) {
  UseMethod("affine")
}

"affine.owin" <- function(X,  mat=diag(c(1,1)), vec=c(0,0), ...,
                          rescue=TRUE) {
  verifyclass(X, "owin")
  vec <- as2vector(vec)
  if(!is.matrix(mat) || any(dim(mat) != c(2,2)))
    stop(paste(sQuote("mat"), "should be a 2 x 2 matrix"))
  diagonalmatrix <- all(mat == diag(diag(mat)))
  scaletransform <- diagonalmatrix && (length(unique(diag(mat))) == 1)
  newunits <- if(scaletransform) unitname(X) else as.unitname(NULL)
  #
  switch(X$type,
         rectangle={
           if(diagonalmatrix) {
             # result is a rectangle
             Y <- owinInternalRect(range(mat[1L,1L] * X$xrange + vec[1L]),
                       range(mat[2L,2L] * X$yrange + vec[2L]))
             unitname(Y) <- newunits
             return(Y)
           } else {
             # convert rectangle to polygon
             P <- as.polygonal(X)
             # call polygonal case
             return(affine.owin(P, mat, vec, rescue=rescue))
           }
         },
         polygonal={
           # Transform the polygonal boundaries
           bdry <- lapply(X$bdry, affinexypolygon, mat=mat, vec=vec,
                          detmat=det(mat))
           # Compile result
           W <- owin(poly=bdry, check=FALSE, unitname=newunits)
           # Result might be a rectangle: if so, convert to rectangle type
           if(rescue)
             W <- rescue.rectangle(W)
           return(W)
         },
         mask={
           #' binary mask
           if(sqrt(abs(det(mat))) < .Machine$double.eps)
             stop("Matrix of linear transformation is singular")
           newpixels <- (length(list(...)) > 0) &&
                        any(dim(X) != rev(spatstat.options("npixel")))
           if(diagonalmatrix && !newpixels) {
             #' diagonal matrix: apply map to row and column locations
             m      <- X$m
             d      <- X$dim
             newbox <- affine(Frame(X), mat=mat, vec=vec)
             xscale <- diag(mat)[1L]
             yscale <- diag(mat)[2L]
             xcol <- xscale * X$xcol + vec[1L]
             yrow <- yscale * X$yrow + vec[2L]
             if(xscale < 0) {
               #' x scale is negative
               xcol <- rev(xcol)
               m <- m[, (d[2L]:1)]
             }
             if(yscale < 0) {
               ## y scale is negative
               yrow <- rev(yrow)
               m <- m[(d[1L]:1), ]
             }
             W <- owin(mask=m,
                       xy=list(x=xcol, y=yrow),
                       xrange=newbox$xrange, yrange=newbox$yrange,
                       unitname=newunits)
           } else {
             #' general case
             #' create box containing transformed window
             newframe <- boundingbox(affinexy(corners(X), mat, vec))
             if(newpixels) {
               W <- as.mask(newframe, ...)
             } else {
               #' determine pixel size
               mp <- mat %*% cbind(c(X$xstep, 0), c(0, X$ystep))
               len <- sqrt(colSums(mp^2))
               xcos <- abs(mp[1,1])/len[1]
               ycos <- abs(mp[1,2])/len[2]
               if(xcos > ycos) {
                 newxstep <- len[1]
                 newystep <- len[2]
               } else {
                 newxstep <- len[2]
                 newystep <- len[1]
               }
               W <- as.mask(newframe, eps=c(newxstep, newystep))
             }
             pixelxy <- rasterxy.mask(W)
             xybefore <- affinexy(pixelxy, mat, vec, invert=TRUE)
             W$m[] <- with(xybefore, inside.owin(x, y, X))
             W <- intersect.owin(W, boundingbox(W))
           }
           if(rescue)
             W <- rescue.rectangle(W)
           return(W)
         },
         stop("Unrecognised window type")
         )
}

"affine.ppp" <- function(X, mat=diag(c(1,1)), vec=c(0,0), ...) {
  verifyclass(X, "ppp")
  vec <- as2vector(vec)
  r <- affinexy(X, mat, vec)
  w <- affine.owin(X$window, mat, vec, ...)
  return(ppp(r$x, r$y, window=w, marks=marks(X, dfok=TRUE), check=FALSE))
}

"affine.im" <- function(X,  mat=diag(c(1,1)), vec=c(0,0), ...) {
  verifyclass(X, "im")
  vec <- as2vector(vec)
  if(!is.matrix(mat) || any(dim(mat) != c(2,2)))
    stop(paste(sQuote("mat"), "should be a 2 x 2 matrix"))
  # Inspect the determinant
  detmat <- det(mat)
  if(sqrt(abs(detmat)) < .Machine$double.eps)
    stop("Matrix of linear transformation is singular")
  #
  diagonalmatrix <- all(mat == diag(diag(mat)))
  scaletransform <- diagonalmatrix && (length(unique(diag(mat))) == 1L)
  newunits <- if(scaletransform) unitname(X) else as.unitname(NULL)
  newpixels <- (length(list(...)) > 0) &&
               any(dim(X) != rev(spatstat.options("npixel")))
  #
  if(diagonalmatrix && !newpixels) {
    # diagonal matrix: apply map to row and column locations
    v      <- X$v
    d      <- X$dim
    newbox <- affine(as.rectangle(X), mat=mat, vec=vec)
    xscale <- diag(mat)[1L]
    yscale <- diag(mat)[2L]
    xcol <- xscale * X$xcol + vec[1L]
    yrow <- yscale * X$yrow + vec[2L]
    if(xscale < 0) {
      # x scale is negative
      xcol <- rev(xcol)
      v <- v[, (d[2L]:1)]
    }
    if(yscale < 0) {
      # y scale is negative
      yrow <- rev(yrow)
      v <- v[(d[1L]:1), ]
    }
    Y <- im(v, xcol=xcol, yrow=yrow,
            xrange=newbox$xrange, yrange=newbox$yrange,
            unitname=newunits)
  } else {
    # general case
    # create box containing transformed image
    newframe <- boundingbox(affinexy(corners(X), mat, vec))
    if(newpixels) {
      W <- as.mask(newframe, ...)
    } else {
      #' determine pixel size
      mp <- mat %*% cbind(c(X$xstep, 0), c(0, X$ystep))
      len <- sqrt(colSums(mp^2))
      xcos <- abs(mp[1,1])/len[1]
      ycos <- abs(mp[1,2])/len[2]
      if(xcos > ycos) {
        newxstep <- len[1]
        newystep <- len[2]
      } else {
        newxstep <- len[2]
        newystep <- len[1]
      }
      W <- as.mask(newframe, eps=c(newxstep, newystep))
    }
    unitname(W) <- newunits
    # raster for transformed image
    naval <- switch(X$type,
                    factor= , 
                    integer = NA_integer_,
                    logical = as.logical(NA_integer_),
                    real = NA_real_,
                    complex = NA_complex_, 
                    character = NA_character_,
                    NA)
    Y <- as.im(W, value=naval)
    # preimages of pixels of transformed image
    xx <- as.vector(rasterx.im(Y))
    yy <- as.vector(rastery.im(Y))
    pre <- affinexy(list(x=xx, y=yy), mat, vec, invert=TRUE)
    # sample original image
    if(X$type != "factor") {
      Y$v[] <- lookup.im(X, pre$x, pre$y, naok=TRUE)
    } else {
      lab <- levels(X)
      lev <- seq_along(lab)
      Y$v[] <- lookup.im(eval.im(as.integer(X)), pre$x, pre$y, naok=TRUE)
      Y <- eval.im(factor(Y, levels=lev, labels=lab))
    }
  }
  return(Y)
}


### ---------------------- reflect ----------------------------------

reflect <- function(X) {
  UseMethod("reflect")
}

reflect.default <- function(X) { affine(X, mat=diag(c(-1,-1))) }

reflect.im <- function(X) {
  stopifnot(is.im(X))
  out <- with(X,
              list(v      = v[dim[1L]:1, dim[2L]:1],
                   dim    = dim,
                   xrange = rev(-xrange),
                   yrange = rev(-yrange),
                   xstep  = xstep,
                   ystep  = ystep,
                   xcol   = rev(-xcol),
                   yrow   = rev(-yrow),
                   type   = type,
                   units  = units))
  class(out) <- "im"
  return(out)
}

### ---------------------- shift ----------------------------------

"shift" <- function(X, ...) {
  UseMethod("shift")
}

shiftxy <- function(X, vec=c(0,0)) {
  vec <- as.numeric(vec)
  n <- length(vec)
  if(n == 0) {
    warning("Null displacement vector; treated as zero")
    return(X)
  } else if(n != 2) {
    stop(paste("Displacement vector has length", n, "!= 2"),
         call.=FALSE)
  }
  list(x = X$x + vec[1L],
       y = X$y + vec[2L])
}

shiftxypolygon <- function(p, vec=c(0,0)) {
  # transform (x,y), retaining other data
  p[c("x","y")] <- shiftxy(p, vec=vec)
  return(p)
}

shift.owin <- function(X,  vec=c(0,0), ..., origin=NULL) {
  verifyclass(X, "owin")
  if(!is.null(origin)) {
    if(!missing(vec))
      warning("argument vec ignored; argument origin has precedence")
    locn <- interpretAsOrigin(origin, X)
    vec <- -locn
  }
  vec <- as2vector(vec)
  # Shift the bounding box
  X$xrange <- X$xrange + vec[1L]
  X$yrange <- X$yrange + vec[2L]
  switch(X$type,
         rectangle={
         },
         polygonal={
           # Shift the polygonal boundaries
           X$bdry <- lapply(X$bdry, shiftxypolygon, vec=vec)
         },
         mask={
           # Shift the pixel coordinates
           X$xcol <- X$xcol + vec[1L]
           X$yrow <- X$yrow + vec[2L]
           # That's all --- the mask entries are unchanged
         },
         stop("Unrecognised window type")
         )
  # tack on shift vector
  attr(X, "lastshift") <- vec
  # units are unchanged
  return(X)
}

shift.ppp <- function(X, vec=c(0,0), ..., origin=NULL) {
  verifyclass(X, "ppp")
  if(!is.null(origin)) {
    if(!missing(vec))
      warning("argument vec ignored; argument origin has precedence")
    locn <- interpretAsOrigin(origin, Window(X))
    vec <- -locn
  }
  vec <- as2vector(vec)
  # perform shift
  r <- shiftxy(X, vec)
  w <- shift.owin(X$window, vec)
  Y <- ppp(r$x, r$y, window=w, marks=marks(X, dfok=TRUE), check=FALSE)
  # tack on shift vector
  attr(Y, "lastshift") <- vec
  return(Y)
}

getlastshift <- function(X) {
  v <- attr(X, "lastshift")
  if(is.null(v))
    stop(paste("Internal error: shifted object of class",
               sQuote(as.character(class(X))[1L]),
               "does not have \"lastshift\" attribute"),
         call.=FALSE)
  if(!(is.numeric(v) && length(v) == 2L))
    stop("Internal error: \"lastshift\" attribute is not a vector",
         call.=FALSE)
  return(v)
}

putlastshift <- function(X, vec) {
  attr(X, "lastshift") <- vec
  return(X)
}

interpretAsOrigin <- function(x, W) {
  if(is.character(x)) {
    x <- paste(x, collapse="")
    x <- match.arg(x,
                   c("centroid", "midpoint",
                     "left", "right", "top", "bottom",
                     "bottomleft", "bottomright",
                     "topleft", "topright"))
    W <- as.owin(W)
    xr <- W$xrange
    yr <- W$yrange
    x <- switch(x,
                centroid    = { unlist(centroid.owin(W)) },
                midpoint    = { c(mean(xr), mean(yr)) },
                left        = { c(xr[1L],   mean(yr)) },
                right       = { c(xr[2L],   mean(yr)) },
                top         = { c(mean(xr), yr[2L]) },
                bottom      = { c(mean(xr), yr[1L]) },
                bottomleft  = { c(xr[1L],   yr[1L]) },
                bottomright = { c(xr[2L],   yr[1L]) },
                topleft     = { c(xr[1L],   yr[2L]) },
                topright    = { c(xr[2L],   yr[2L]) },
                stop(paste("Unrecognised option",sQuote(x)),
                     call.=FALSE))
  } 
  return(as2vector(x))
}

### ---------------------- scalar dilation ---------------------------------

scalardilate <- function(X, f, ...) {
  UseMethod("scalardilate")
}

scalardilate.default <- function(X, f, ...) {
  trap.extra.arguments(..., .Context="In scalardilate(X,f)")
  check.1.real(f, "In scalardilate(X,f)")
  stopifnot(is.finite(f) && f > 0)
  Y <- affine(X, mat=diag(c(f,f)))
  return(Y)
}

scalardilate.breakpts <- function(X, f, ...) {
  out <- with(X,
              list(val    = f*val,
                   max    = f*max,
                   ncells = ncells,
                   r      = f*r,
                   even   = even,
                   npos   = npos,
                   step   = if(is.null(step)) NULL else (f*step)))
  class(out) <- "breakpts"
  out
}  
                            
scalardilate.im <- scalardilate.owin <- scalardilate.psp <- scalardilate.ppp <-
  function(X, f, ..., origin=NULL) {
  trap.extra.arguments(..., .Context="In scalardilate(X,f)")
  check.1.real(f, "In scalardilate(X,f)")
  stopifnot(is.finite(f) && f > 0)
  if(!is.null(origin)) {
    X <- shift(X, origin=origin)
    negorig <- getlastshift(X)
  } else negorig <- c(0,0)
  Y <- affine(X, mat=diag(c(f, f)), vec = -negorig)
  return(Y)
}
  


