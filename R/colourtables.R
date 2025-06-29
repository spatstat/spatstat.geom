#
# colourtables.R
#
# support for colour maps and other lookup tables
#
# $Revision: 1.70 $ $Date: 2025/06/28 02:23:09 $
#

colourmap <- function(col, ..., range=NULL, breaks=NULL, inputs=NULL, gamma=1,
                      compress=NULL, decompress=NULL) {
  if(nargs() == 0) {
    ## null colour map
    f <- lut()
  } else {
    ## validate colour data 
    col2hex(col)
    ## store without conversion
    f <- lut(col, ..., range=range, breaks=breaks, inputs=inputs, gamma=gamma,
             compress=compress, decompress=decompress)
  }
  class(f) <- c("colourmap", class(f))
  f
}

lut <- function(outputs, ..., range=NULL, breaks=NULL, inputs=NULL,
                gamma=1, compress=NULL, decompress=NULL) {
  if(nargs() == 0) {
    ## null lookup table
    f <- function(x, what="value"){NULL}
    class(f) <- c("lut", class(f))
    attr(f, "stuff") <- list(n=0)
    return(f)
  }
  if(is.null(gamma)) gamma <- 1
  if(is.null(compress)) {
    decompress <- NULL
  } else {
    stopifnot(is.function(compress))
    if(!is.null(decompress)) {
      stopifnot(is.function(decompress))
    } else {
      ## Argument decompress is missing
      ## Try to construct it from 'compress'
      if(is.primitive(compress) && samefunction(compress, log10)) {
        ## logarithmic lookup table
        decompress <- TenPower
      } else if(inherits(compress, c("ecdf", "ewcdf", "interpolatedCDF"))) {
        ## histogram-equalised lookup table
        decompress <- quantilefun(compress)
      } else {
        ## not recognised
        stop("Argument 'decompress' is required", call.=FALSE)
      }
    } 
  }
  n <- length(outputs)
  given <- c(!is.null(range), !is.null(breaks), !is.null(inputs))
  names(given) <- nama <- c("range", "breaks", "inputs")
  ngiven <- sum(given)
  if(ngiven == 0L)
    stop(paste("One of the arguments",
               commasep(sQuote(nama), "or"),
               "should be given"))
  if(ngiven > 1L) {
    offending <- nama[given]
    stop(paste("The arguments",
               commasep(sQuote(offending)),
               "are incompatible"))
  }
  if(!is.null(inputs)) {
    #' discrete set of input values mapped to output values
    if(n == 1L) {
      #' constant output
      n <- length(inputs)
      outputs <- rep(outputs, n)
    } else stopifnot(length(inputs) == length(outputs))
    stuff <- list(n=n, discrete=TRUE, inputs=inputs, outputs=outputs,
                  compress=compress, decompress=decompress)
    f <- function(x, what="value") {
      inputs <- stuff$inputs
      if(is.function(compress <- stuff$compress)) {
        x <- compress(x)
        inputs <- compress(inputs)
      }
      m <- match(x, inputs)
      if(what == "index")
        return(m)
      cout <- stuff$outputs[m]
      return(cout)
    }
  } else {
    #' range of numbers, or date/time interval, mapped to colours
    #' determine type of domain
    timeclasses <- c("Date", "POSIXt")
    is.time <- inherits(range, timeclasses) || inherits(breaks, timeclasses)
    if(is.null(breaks)) {
      #' determine breaks
      if(is.null(compress)) {
        breaks <- gammabreaks(range, n + 1L, gamma)
      } else {
        breaks <- decompress(gammabreaks(compress(range), n + 1L, gamma))
      }
      gamma.used <- gamma
    } else {
      #' check user-specified breaks
      stopifnot(length(breaks) >= 2)
      if(length(outputs) == 1L) {
        n <- length(breaks) - 1L
        outputs <- rep(outputs, n)
      } else stopifnot(length(breaks) == length(outputs) + 1L)
      if(!all(diff(breaks) > 0))
        stop("breaks must be increasing")
      gamma.used <- NULL
    }
    stuff <- list(n=n, discrete=FALSE, breaks=breaks, outputs=outputs,
                  gamma=gamma.used, compress=compress, decompress=decompress)
    #' use appropriate function
    if(is.time) {
      f <- function(x, what="value") {
        breaks <- stuff$breaks
        if(is.function(compress <- stuff$compress)) {
          x <- compress(x)
          breaks <- compress(breaks)
        }
        x <- as.vector(as.numeric(x))
        z <- findInterval(x, breaks, rightmost.closed=TRUE)
        oo <- stuff$outputs
        z[z <= 0 | z > length(oo)] <- NA
        if(what == "index")
          return(z)
        cout <- oo[z]
        return(cout)
      }
    } else {
      f <- function(x, what="value") {
        breaks <- stuff$breaks
        if(is.function(compress <- stuff$compress)) {
          x <- compress(x)
          breaks <- compress(breaks)
        }
        stopifnot(is.numeric(x))
        x <- as.vector(x)
        z <- findInterval(x, breaks, rightmost.closed=TRUE)
        oo <- stuff$outputs
        z[z <= 0 | z > length(oo)] <- NA
        if(what == "index")
          return(z)
        cout <- stuff$outputs[z]
        return(cout)
      }
    }
  }
  attr(f, "stuff") <- stuff
  class(f) <- c("lut", class(f))
  f
}

print.lut <- function(x, ...) {
  if(inherits(x, "colourmap")) {
    tablename <- "Colour map"
    outputname <- "colour"
  } else {
    tablename  <- "Lookup table"
    outputname <- "output"
  }
  stuff <- attr(x, "stuff")
  n <- stuff$n
  if(n == 0) {
    ## Null map
    cat(paste("Null", tablename, "\n"))
    return(invisible(NULL))
  }
  if(stuff$discrete) {
    cat(paste(tablename, "for discrete set of input values\n"))
    out <- data.frame(input=stuff$inputs, output=stuff$outputs)
  } else {
    b <- stuff$breaks
    cat(paste(tablename, "for the range", prange(b[c(1L,n+1L)]), "\n"))
    leftend  <- rep("[", n)
    rightend <- c(rep(")", n-1), "]")
    inames <- paste(leftend, b[-(n+1L)], ", ", b[-1L], rightend, sep="")
    out <- data.frame(interval=inames, output=stuff$outputs)
  }
  colnames(out)[2L] <- outputname
  print(out)
  if(!is.null(gamma <- stuff$gamma) && gamma != 1)
    splat("Generated using gamma =", gamma)
  if(!is.null(compress <- stuff$compress)) {
    if(samefunction(compress, log10)) {
      splat("Logarithmic", tolower(tablename))
    } else {
      splat("Input compression map:")
      print(stuff$compress)
    }
  }
  invisible(NULL)
}

print.colourmap <- function(x, ...) {
  NextMethod("print")
}

summary.lut <- function(object, ...) {
  s <- attr(object, "stuff")
  if(inherits(object, "colourmap")) {
    s$tablename <- "Colour map"
    s$outputname <- "colour"
  } else {
    s$tablename  <- "Lookup table"
    s$outputname <- "output"
  }
  class(s) <- "summary.lut"
  return(s)
}

print.summary.lut <- function(x, ...) {
  n <- x$n
  if(n == 0) {
    cat(paste("Null", x$tablename, "\n"))
    return(invisible(NULL))
  }
  if(x$discrete) {
    cat(paste(x$tablename, "for discrete set of input values\n"))
    out <- data.frame(input=x$inputs, output=x$outputs)
  } else {
    b <- x$breaks
    cat(paste(x$tablename, "for the range", prange(b[c(1L,n+1L)]), "\n"))
    leftend  <- rep("[", n)
    rightend <- c(rep(")", n-1L), "]")
    inames <- paste(leftend, b[-(n+1L)], ", ", b[-1L], rightend, sep="")
    out <- data.frame(interval=inames, output=x$outputs)
  }
  colnames(out)[2L] <- x$outputname
  print(out)  
  if(!is.null(compress <- x$compress)) {
    if(samefunction(compress, log10)) {
      splat("Logarithmic", tolower(x$tablename))
    } else {
      splat("Input compression map:")
      print(compress)
    }
  }
  return(invisible(NULL))
}

plot.colourmap <- local({

  # recognised additional arguments to image.default() and axis()
  
  imageparams <- c("main", "asp", "sub", "axes", "ann",
                   "cex", "font", 
                   "cex.axis", "cex.lab", "cex.main", "cex.sub",
                   "col.axis", "col.lab", "col.main", "col.sub",
                   "font.axis", "font.lab", "font.main", "font.sub")
  axisparams <- c("cex", 
                  "cex.axis", "cex.lab",
                  "col.axis", "col.lab",
                  "font.axis", "font.lab",
                  "las", "mgp", "xaxp", "yaxp",
                  "tck", "tcl", "xpd")

  linmap <- function(x, from, to) {
    dFrom <- as.numeric(diff(from))
    dTo <- as.numeric(diff(to))
    b <- dTo/dFrom
    if(is.nan(b)) b <- 0
    if(!is.finite(b)) stop("Internal error: Cannot map zero width interval")
    to[1L] + b * (x - from[1L])
  }

  ensurenumeric <- function(x) { if(is.numeric(x)) x else as.numeric(x) }

  # rules to determine the ribbon dimensions when one dimension is given
  widthrule <- function(heightrange, separate, n, gap) {
    dh <- diff(heightrange)
    if(separate || dh == 0) 1 else dh/10
  }
  heightrule <- function(widthrange, separate, n, gap) {
    dw <- diff(widthrange)
    if(dw == 0) 1 else (dw * (if(separate) (n + (n-1)*gap) else 10))
  }

  Ticks <- function(usr, log=FALSE, nint=NULL, ..., clip=TRUE) {
    #' modification of grDevices::axisTicks
    #'      constrains ticks to be inside the specified range 'usr' if clip=TRUE
    #'      accepts nint=NULL as if it were missing
    z <- if(is.null(nint)) axisTicks(usr=usr, log=log, ...) else
         axisTicks(usr=usr, log=log, nint=nint, ...) 
    if(clip) {
      zlimits <- if(log) 10^usr else usr
      z <- z[inside.range(z, zlimits)]
    }
    return(z)
  }

    
  plot.colourmap <- function(x, ..., main,
                             xlim=NULL, ylim=NULL,
                             vertical=FALSE,
                             axis=TRUE,
                             side = if(vertical) "right" else "bottom",
                             labelmap=NULL, gap=0.25, add=FALSE,
                             increasing=NULL, nticks=5, at=NULL, box=NULL) {
    if(missing(main))
      main <- short.deparse(substitute(x))
    if(missing(vertical) && !missing(side)) 
      vertical <- (sideCode(side) %in% c(2, 4))
    
    dotargs <- list(...)
    if(inherits(dotargs$col, "colourmap"))
      dotargs <- dotargs[names(dotargs) != "col"]
    
    stuff <- attr(x, "stuff")
    col <- stuff$outputs
    n   <- stuff$n
    if(n == 0) {
      ## Null map
      return(invisible(NULL))
    }
    discrete <- stuff$discrete
    if(discrete) {
      check.1.real(gap, "In plot.colourmap")
      explain.ifnot(gap >= 0, "In plot.colourmap")
    }
    separate <- discrete && (gap > 0)

    compress <- stuff$compress %orifnull% identity
    decompress <- stuff$decompress %orifnull% identity
    is.log <- samefunction(compress, log10)

    if(is.null(labelmap)) {
      labelmap <- identity
    } else if(is.numeric(labelmap) && length(labelmap) == 1L && !discrete) {
      labscal <- labelmap
      labelmap <- function(x) { x * labscal }
    } else stopifnot(is.function(labelmap))

    #' map values back to original scale
    Labelmap <- function(x) { labelmap(decompress(x)) }

    if(is.null(increasing))
      increasing <- !(discrete && vertical)
    reverse <- !increasing

    #' determine pixel entries 'v' and colour map breakpoints 'bks'
    #' to be passed to 'image.default'
    trivial <- FALSE
    if(!discrete) {
      # real numbers: continuous ribbon
      bks <- compress(stuff$breaks)
      rr <- range(bks)
      trivial <- (diff(rr) == 0)
      v <- if(trivial) rr[1] else
           seq(from=rr[1L], to=rr[2L], length.out=max(n+1L, 1024))
    } else if(!separate) {
      # discrete values: blocks of colour, run together
      v <- (1:n) - 0.5
      bks <- 0:n
      rr <- c(0,n)
    } else {
      # discrete values: separate blocks of colour
      vleft <- (1+gap) * (0:(n-1L))
      vright <- vleft + 1
      v <- vleft + 0.5
      rr <- c(0, n + (n-1)*gap)
    }
    # determine position of ribbon or blocks of colour
    if(is.null(xlim) && is.null(ylim)) {
      u <- widthrule(rr, separate, n, gap)
      if(!vertical) {
        xlim <- rr
        ylim <- c(0,u)
      } else {
        xlim <- c(0,u)
        ylim <- rr
      }
    } else if(is.null(ylim)) {
      if(!vertical) 
        ylim <- c(0, widthrule(xlim, separate, n, gap))
      else 
        ylim <- c(0, heightrule(xlim, separate, n, gap))
    } else if(is.null(xlim)) {
      if(!vertical) 
        xlim <- c(0, heightrule(ylim, separate, n, gap))
      else 
        xlim <- c(0, widthrule(ylim, separate, n, gap))
    } 

    # .......... initialise plot ...............................
    if(!add)
      do.call.matched(plot.default,
                      resolve.defaults(list(x=xlim, y=ylim,
                                            type="n", main=main,
                                            axes=FALSE, xlab="", ylab="",
                                            asp=1.0),
                                       dotargs))

    if(separate) {
      # ................ plot separate blocks of colour .................
      if(reverse) 
        col <- rev(col)
      if(!vertical) {
        # horizontal arrangement of blocks
        xleft <- linmap(vleft, rr, xlim)
        xright <- linmap(vright, rr, xlim)
        y <- ylim
        z <- matrix(1, 1L, 1L)
        for(i in 1:n) {
          x <- c(xleft[i], xright[i])
          do.call.matched(image.default,
                          resolve.defaults(list(x=ensurenumeric(x),
                                                y=ensurenumeric(y),
                                                z=z,
                                                add=TRUE,
                                                col=col[i]),
                                           dotargs),
                          extrargs=imageparams)
        }
      } else {
        # vertical arrangement of blocks
        x <- xlim 
        ylow <- linmap(vleft, rr, ylim)
        yupp <- linmap(vright, rr, ylim)
        z <- matrix(1, 1L, 1L)
        for(i in 1:n) {
          y <- c(ylow[i], yupp[i])
          do.call.matched(image.default,
                          resolve.defaults(list(x=ensurenumeric(x),
                                                y=ensurenumeric(y),
                                                z=z,
                                                add=TRUE,
                                                col=col[i]),
                                           dotargs),
                          extrargs=imageparams)
        }
      }
    } else {
      # ................... plot ribbon image .............................
      if(!vertical) {
        # horizontal colour ribbon
        x <- linmap(v, rr, xlim)
        y <- ylim
        z <- matrix(v, ncol=1L)
      } else {
        # vertical colour ribbon
        y <- linmap(v, rr, ylim)
        z <- matrix(v, nrow=1L)
        x <- xlim
      }
      #' deal with Date or integer values
      x <- ensurenumeric(x)
      if(!trivial) {
        if(any(diff(x) == 0)) 
          x <- seq(from=x[1L], to=x[length(x)], length.out=length(x))
        y <- ensurenumeric(y)
        if(any(diff(y) == 0)) 
          y <- seq(from=y[1L], to=y[length(y)], length.out=length(y))
        bks <- ensurenumeric(bks)
        if(any(diff(bks) <= 0)) {
          ok <- (diff(bks) > 0)
          bks <- bks[ok]
          col <- col[ok]
        }
      }
      if(reverse)
        col <- rev(col)
      do.call.matched(image.default,
                      resolve.defaults(list(x=x, y=y, z=z,
                                            add=TRUE,
                                            breaks=ensurenumeric(bks),
                                            col=col),
                                       dotargs),
                      extrargs=imageparams)
    }
    #' draw box around colours?
    #' default is TRUE unless drawing blocks of colour with gaps between.
    if(is.null(box)) box <- !separate
    if(!isFALSE(box))
      rect(xlim[1], ylim[1], xlim[2], ylim[2])

    if(axis) {
      # ................. draw annotation ..................
      if(!vertical) {
          # add horizontal axis/annotation
        if(discrete) {
          la <- paste(Labelmap(stuff$inputs))
          at <- linmap(v, rr, xlim)
        } else {
          atpos <- compress(at %orifnull% Ticks(rr, log=is.log, nint=nticks))
          at <- linmap(atpos, rr, xlim)
          la <- Labelmap(atpos)
        }
        if(reverse)
          at <- rev(at)
        # default axis position is below the ribbon (side=1)
        sidecode <- sideCode(side)
        if(!(sidecode %in% c(1L,3L)))
          warning(paste("side =", 
                        if(is.character(side)) sQuote(side) else side,
                        "is not consistent with horizontal orientation"))
        pos <- c(ylim[1L], xlim[1L], ylim[2L], xlim[2L])[sidecode]
        # don't draw axis lines if plotting separate blocks
        lwd0 <- if(separate) 0 else 1
        # draw axis
        do.call.matched(graphics::axis,
                        resolve.defaults(dotargs,
                                         list(side = sidecode,
                                              pos = pos,
                                              at = ensurenumeric(at),
                                              labels=la, lwd=lwd0)),
                        extrargs=axisparams)
      } else {
        # add vertical axis
        if(discrete) {
          la <- paste(Labelmap(stuff$inputs))
          at <- linmap(v, rr, ylim)
        } else {
          atpos <- compress(at %orifnull% Ticks(rr, log=is.log, nint=nticks))
          at <- linmap(atpos, rr, ylim)
          la <- Labelmap(atpos)
        }
        if(reverse)
          at <- rev(at)
        # default axis position is to the right of ribbon (side=4)
        sidecode <- sideCode(side)
        if(!(sidecode %in% c(2L,4L)))
          warning(paste("side =",
                        if(is.character(side)) sQuote(side) else side,
                        "is not consistent with vertical orientation"))
        pos <- c(ylim[1L], xlim[1L], ylim[2L], xlim[2L])[sidecode]
        # don't draw axis lines if plotting separate blocks
        lwd0 <- if(separate) 0 else 1
        # draw labels horizontally if plotting separate blocks
        las0 <- if(separate) 1 else 0
        # draw axis
        do.call.matched(graphics::axis,
                        resolve.defaults(dotargs,
                                         list(side=sidecode,
                                              pos=pos,
                                              at=ensurenumeric(at),
                                              labels=la, lwd=lwd0, las=las0)),
                        extrargs=axisparams)
      }
    }
    invisible(NULL)
  }

  plot.colourmap
})


# Interpolate a colourmap or lookup table defined on real numbers

interp.colourmap <- function(m, n=512) {
  if(!inherits(m, "colourmap"))
    stop("m should be a colourmap")
  st <- attr(m, "stuff")
  if(st$discrete) {
    # discrete set of input values mapped to colours
    xknots <- st$inputs
    # Ensure the inputs are real numbers
    if(!is.numeric(xknots))
      stop("Cannot interpolate: inputs are not numerical values")
  } else {
    # interval of real line, chopped into intervals, mapped to colours
    # Find midpoints of intervals
    bks <- st$breaks
    nb <- length(bks)
    xknots <- (bks[-1L] + bks[-nb])/2
  }
  # corresponding colours in hsv coordinates
  yknots.hsv <- rgb2hsva(col2rgb(st$outputs, alpha=TRUE))
  # transform 'hue' from polar to cartesian coordinate
  # divide domain into n equal intervals
  xrange <- range(xknots)
  xbreaks <- seq(xrange[1L], xrange[2L], length=n+1L)
  xx <- (xbreaks[-1L] + xbreaks[-(n+1L)])/2
  # interpolate saturation and value in hsv coordinates
  yy.sat <- approx(x=xknots, y=yknots.hsv["s", ], xout=xx)$y
  yy.val <- approx(x=xknots, y=yknots.hsv["v", ], xout=xx)$y
  # interpolate hue by first transforming polar to cartesian coordinate
  yknots.hue <- 2 * pi * yknots.hsv["h", ]
  yy.huex <- approx(x=xknots, y=cos(yknots.hue), xout=xx)$y
  yy.huey <- approx(x=xknots, y=sin(yknots.hue), xout=xx)$y
  yy.hue <- (atan2(yy.huey, yy.huex)/(2 * pi)) %% 1
  # handle transparency
  yknots.alpha <- yknots.hsv["alpha", ]
  if(all(yknots.alpha == 1)) {
    ## opaque colours: form using hue, sat, val
    yy <- hsv(yy.hue, yy.sat, yy.val)
  } else {
    ## transparent colours: interpolate alpha
    yy.alpha <- approx(x=xknots, y=yknots.alpha, xout=xx)$y
    ## form colours using hue, sat, val, alpha
    yy <- hsv(yy.hue, yy.sat, yy.val, yy.alpha)    
  }
  # done
  f <- colourmap(yy, breaks=xbreaks,
                 compress=st$compress, decompress=st$decompress)
  return(f)
}

interp.colours <- function(x, length.out=512) {
  y <- colourmap(x, range=c(0,1))
  z <- interp.colourmap(y, length.out)
  oo <- attr(z, "stuff")$outputs
  return(oo)
}

tweak.colourmap <- local({

  is.hex <- function(z) {
    is.character(z) &&
    all(nchar(z, keepNA=TRUE) %in% c(7L,9L)) &&
    identical(substr(z, 1L, 7L), substr(col2hex(z), 1L, 7L))
  }

  tweak.colourmap <- function(m, col, ..., inputs=NULL, range=NULL) {
    if(!inherits(m, "colourmap"))
      stop("m should be a colourmap")
    if(is.null(inputs) && is.null(range)) 
      stop("Specify either inputs or range")
    if(!is.null(inputs) && !is.null(range))
      stop("Do not specify both inputs and range")
    ## determine indices of colours to be changed
    if(!is.null(inputs)) {
      ix <- m(inputs, what="index")
    } else {
      if(!(is.numeric(range) && length(range) == 2 && diff(range) > 0))
        stop("range should be a numeric vector of length 2 giving (min, max)")
      if(length(col2hex(col)) != 1L)
        stop("When range is given, col should be a single colour value")
      ixr <- m(range, what="index")
      ix <- (ixr[1L]):(ixr[2L])
    }
    ## reassign colours
    st <- attr(m, "stuff")
    outputs <- st$outputs
    result.hex <- FALSE
    if(is.hex(outputs)) {
      ## convert replacement data to hex
      col <- col2hex(col)
      result.hex <- TRUE
    } else if(is.hex(col)) {
      ## convert existing data to hex
      outputs <- col2hex(outputs)
      result.hex <- TRUE
    } else if(!(is.character(outputs) && is.character(col))) {
      ## unrecognised format - convert both to hex
      outputs <- col2hex(outputs)
      col <- col2hex(col)
      result.hex <- TRUE
    }
    if(result.hex) {
      ## hex codes may be 7 or 9 characters
      outlen <- nchar(outputs)
      collen <- nchar(col)
      if(length(unique(c(outlen, collen))) > 1L) {
        ## convert all to 9 characters
        if(any(bad <- (outlen == 7))) 
          outputs[bad] <- paste0(outputs[bad], "FF")
        if(any(bad <- (collen == 7))) 
          col[bad] <- paste0(col[bad], "FF")
      }
    }
    ## Finally, replace
    outputs[ix] <- col
    st$outputs <- outputs
    attr(m, "stuff") <- st
    assign("stuff", st, envir=environment(m))
    return(m)
  }

  tweak.colourmap
})

colouroutputs <- function(x) {
  stopifnot(inherits(x, "colourmap"))
  attr(x, "stuff")$outputs
}

"colouroutputs<-" <- function(x, value) {
  stopifnot(inherits(x, "colourmap"))
  st <- attr(x, "stuff")
  col2hex(value) # validates colours
  st$outputs[] <- value
  attr(x, "stuff") <- st
  assign("stuff", st, envir=environment(x))
  return(x)
}

rev.colourmap <- function(x) {
  colouroutputs(x) <- rev(colouroutputs(x))
  return(x)
}

restrict.colourmap <- function(x, ..., range=NULL, breaks=NULL, inputs=NULL) {
  stopifnot(inherits(x, "colourmap"))
  given <- c(!is.null(range), !is.null(breaks), !is.null(inputs))
  names(given) <- nama <- c("range", "breaks", "inputs")
  ngiven <- sum(given)
  if(ngiven == 0L)
    return(x)
  if(ngiven > 1L) {
    offending <- nama[given]
    stop(paste("The arguments",
               commasep(sQuote(offending)),
               "are incompatible"))
  }
  stuff <- attr(x, "stuff")
  if(!is.null(inputs)) {
    ## discrete colour map
    if(!stuff$discrete) 
      stop("Cannot update 'inputs'; the existing colour map is not discrete",
           call.=FALSE)
    oldinputs <- stuff$inputs
    oldoutputs <- stuff$outputs
    m <- match(inputs, oldinputs)
    if(any(is.na(m)))
      stop("New inputs are not a subset of the old inputs", call.=FALSE)
    result <- colourmap(oldoutputs[m], inputs=inputs,
                        compress=stuff$compress,
                        decompress=stuff$decompress)
  } else if(!is.null(range)) {
    ## colour map for continuous domain
    ## range specified
    if(stuff$discrete) 
      stop("Cannot update 'range'; the existing colour map is discrete",
           call.=FALSE)
    check.range(range)
    oldbreaks <- stuff$breaks
    if(!all(inside.range(range, range(oldbreaks))))
      stop("new range of values is not a subset of current range")
    ## restrict existing breaks to new range
    newbreaks <- pmax(range[1], pmin(range[2], oldbreaks))
    newbreaks <- unique(newbreaks)
    ## evaluate current colour at midpoint of each new interval
    newmid <- newbreaks[-length(newbreaks)] + diff(newbreaks)/2
    newout <- x(newmid)
    result <- colourmap(newout, breaks=newbreaks,
                        compress=stuff$compress,
                        decompress=stuff$decompress)
  } else {
    ## colour map for continuous domain
    ## breaks specified
    if(stuff$discrete) 
      stop("Cannot update 'breaks'; the existing colour map is discrete",
           call.=FALSE)
    oldbreaks <- stuff$breaks
    if(!all(inside.range(breaks, range(oldbreaks))))
      stop("new range of 'breaks' is not a subset of current range of 'breaks'",
           call.=FALSE)
    newmid <- breaks[-length(breaks)] + diff(breaks)/2
    newout <- x(newmid)
    result <- colourmap(newout, breaks=breaks,
                        compress=stuff$compress,
                        decompress=stuff$decompress)
  }
  return(result)
}

as.colourmap <- function(x, ...) {
  UseMethod("as.colourmap")
}

as.colourmap.colourmap <- function(x, ...) {
  #' remove attributes which are not part of class 'colourmap'
  atr <- attributes(x)
  attributes(x) <- atr[names(atr) %in% c("stuff", "class")]
  return(x)
}

