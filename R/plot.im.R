#
#   plot.im.R
#
#  $Revision: 1.176 $   $Date: 2025/06/30 04:12:26 $
#
#  Plotting code for pixel images
#
#  plot.im
#  image.im
#  contour.im
#
###########################################################################

plot.im <- local({

  ## auxiliary functions

  image.doit <- function(imagedata, ...,
                         add=FALSE, show.all=!add,
                         extrargs=graphicsPars("image"), W,
                         addcontour=FALSE, contourargs=list(),
                         args.contour=list(), # legacy - undocumented
                         workaround=FALSE) {
    aarg <- resolve.defaults(..., list(add=add, show.all=show.all))
    ##
    if(add && show.all) {
      ## set up the window space *with* the main title
      ## using the same code as plot.owin, for consistency
      force(W)
      do.call.matched(plot.owin,
                      resolve.defaults(list(x=quote(W), type="n"), aarg), 
                      extrargs=graphicsPars("owin"))
    }
    if(workaround && isTRUE(aarg$useRaster)) {
      #' workaround for bug 16035
      #' detect reversed coordinates
      usr <- par('usr')
      xrev <- (diff(usr[1:2]) < 0) 
      yrev <- (diff(usr[3:4]) < 0)
      if(xrev || yrev) {
        #' flip matrix of pixel values, because the device driver does not
        z <- imagedata$z
        d <- dim(z) # z is in the orientation expected for image.default
        if(xrev) z <- z[d[1]:1,       , drop=FALSE]
        if(yrev) z <- z[      , d[2]:1, drop=FALSE]
        imagedata$z <- z
      }
    }
    extrargs <- setdiff(extrargs, c("claim.title.space", "box"))
    if(!is.na(k <- match("adj.main", names(aarg))))
      names(aarg)[k] <- "adj"
    z <- do.call.matched(image.default,
                         append(imagedata, aarg),
                         extrargs=extrargs)
    if(addcontour) {
      do.call(do.contour,
              resolve.defaults(imagedata,
                               list(add=TRUE),
                               contourargs,
                               args.contour,
                               list(col=par('fg')),
                               aarg,
                               .StripNull=TRUE))
    }
    
    return(z)
  }

  do.contour <- function(x, y, z, ...,
                         nlevels=10, levels=NULL, labels=NULL, drawlabels=TRUE, 
                         values.are.log=FALSE) {
    nx <- length(x)
    ny <- length(y)
    nz <- dim(z)
    if(nx > nz[1]) {
      if(nz[1] == 1) {
        z <- rbind(z, z)
        nz <- dim(z)
        drawlabels <- FALSE
      } else {
        x <- (x[-1] + x[-nx])/2
        nx <- nx-1
      }
    }
    if(ny > nz[2]) {
      if(nz[2] == 1) {
        z <- cbind(z, z)
        nz <- dim(z)
        drawlabels <- FALSE
      } else {
        y <- (y[-1] + y[-ny])/2
        ny <- ny-1
      }
    }
    if(values.are.log) {
      ## ................  logarithmic case ........................
      ## z is log10 of actual value
      if(!is.null(levels)) {
        labels <- paste(levels)
        levels <- log10(levels)
      } else {
        logra <- range(z, finite=TRUE)
        ## default levels commensurate with logarithmic colour scale
        dlr <- diff(logra)
        if(dlr > 1.5) {
          ## usual case - data ranges over several powers of 10
          wholepowers <- 10^(floor(logra[1]):ceiling(logra[2]))
          levelsperdecade <- nlevels/max(1, length(wholepowers)-1)
          if(levelsperdecade >= 1) {
            ## At least one contour level for every power of 10
            ## Decide on leading digits
            if(levelsperdecade < 1.5) {
              leadingdigits <- 1
            } else if(levelsperdecade < 2.5) {
              leadingdigits <- c(1,3)
            } else if(levelsperdecade < 3.5) {
              leadingdigits <- c(1,2,5)
            } else {
              ## use fractional powers of 10, equally spaced on log scale
              leadingdigits <- 10^seq(0, 1, length.out=ceiling(levelsperdecade)+1)
              leadingdigits <- leadingdigits[-length(leadingdigits)]
              leadingdigits <- round(leadingdigits, max(2, ceiling(log10(levelsperdecade))))
            }
            explevels <- sort(unique(as.numeric(outer(wholepowers, leadingdigits, "*"))))
          } else {
            ## more than one power of 10 between successive contour levels
            explevels <- wholepowers
            thinperiod <- max(1L, floor(1/levelsperdecade))
            if(thinperiod > 1) {
              i <- floor((thinperiod+1)/2)
              explevels <- explevels[seq_along(explevels) %% thinperiod == i]
            }
          }
        } else {
          ## Small range: use standard (non-logarithmic scale) values
          explevels <- pretty(10^logra, nlevels)
        }
        ## restrict to actual range
        explevels <- explevels[inside.range(explevels, 10^logra)]
        if(length(explevels) == 0) explevels <- 10^mean(logra)
        ## finally define levels and labels
        labels <- paste(explevels)
        levels <- log10(explevels)
      }
      ## ................ end logarithmic case ........................
    }
    do.call.matched(contour.default,
                    resolve.defaults(list(x=x, y=y, z=z,
                                          ...,
                                          nlevels=nlevels,
                                          levels=levels,
                                          labels=labels,
                                          drawlabels=drawlabels),
                                     .StripNull=TRUE))
  }
                 
  do.box.etc <- function(bb, add, argh) {
    do.call(box.etc, append(list(bb=bb, add=add), argh))
  }
  
  box.etc <- function(bb, ..., add=FALSE, box=!add,
                      axes=FALSE, ann=FALSE, xlab="", ylab="") {
    # axes for image
    xr <- bb$xrange
    yr <- bb$yrange
    if(box)
      rect(xr[1], yr[1], xr[2], yr[2])
    if(axes) {
      px <- prettyinside(xr)
      py <- prettyinside(yr)
      do.call.plotfun(graphics::axis,
                      resolve.defaults(
                                       list(side=1, at=px), 
                                       list(...),
                                       list(pos=yr[1])),
                      extrargs=graphicsPars("axis"))
      do.call.plotfun(graphics::axis,
                      resolve.defaults(
                                       list(side=2, at=py), 
                                       list(...),
                                       list(pos=xr[1])),
                      extrargs=graphicsPars("axis"))
    }
    ## axis labels xlab, ylab
    if(ann) {
      dox <- any(nzchar(xlab))
      doy <- any(nzchar(ylab))
      line0 <- if(axes) 1 else 0
      if(dox || doy) {
        mtargs <- resolve.defaults(list(...), list(line=line0))
        if(dox)
          do.call.matched(mtext, append(list(text=xlab, side=1), mtargs))
        if(doy)
          do.call.matched(mtext, append(list(text=ylab, side=2), mtargs))
      }
    }
  }
  
  
  clamp <- function(x, v, tol=0.02 * diff(v)) {
    ok <- (x >= v[1] - tol) & (x <= v[2] + tol)
    x[ok]
  }
  
  cellbreaks <- function(x, dx) {
    nx <- length(x)
    seq(x[1] - dx/2, x[nx] + dx/2, length.out=nx+1)
  }

  log10orNA <- function(x) {
    y <- rep(NA_real_, length(x))
    ok <- !is.na(x) & (x > 0)
    y[ok] <- log10(x[ok])
    return(y)
  }

  Ticks <- function(usr, log=FALSE, nint=NULL, ..., clip=TRUE, deco=identity) {
    #' modification of grDevices::axisTicks
    #'      constrains ticks to be inside the specified range if clip=TRUE
    #'      accepts nint=NULL as if it were missing
    z <- if(is.null(nint)) axisTicks(usr=usr, log=log, ...) else
         axisTicks(usr=usr, log=log, nint=nint, ...) 
    if(clip) {
      zlimits <- if(log) 10^usr else usr
      z <- z[inside.range(z, zlimits)]
    }
    if(!log) z <- deco(z)
    return(unique(z))
  }

  numericalRange <- function(x, zlim=NULL) {
    xr <- suppressWarnings(range(x, finite=TRUE))
    if(!all(is.finite(xr)))
      warning("All pixel values are NA", call.=FALSE)
    if(!is.null(zlim)) 
      xr <- suppressWarnings(range(xr, zlim, finite=TRUE))
    if(!all(is.finite(xr))) {
      warning("Cannot determine range of values for colour map",
              call.=FALSE)
      xr <- c(0,0)
    }
    return(xr)
  }
  
  # main function
  PlotIm <- function(x, ...,
                     main, 
                     add=FALSE, clipwin=NULL,
                     col=NULL, reverse.col=FALSE,
                     valuesAreColours=NULL, log=FALSE,
                     ncolours=256, gamma=1, 
                     ribbon=show.all, show.all=!add,
                     drop.ribbon=FALSE,
                     ribside=c("right", "left", "bottom", "top"),
                     ribsep=0.15, ribwid=0.05, ribn=1024,
                     ribscale=1, ribargs=list(), riblab=NULL, colargs=list(),
                     useRaster=NULL, workaround=FALSE, zap=1,
                     do.plot=TRUE,
                     addcontour=FALSE, contourargs=list(),
                     background=NULL, clip.background=FALSE) {
    if(missing(main)) main <- short.deparse(substitute(x))
    verifyclass(x, "im")
    force(show.all)
    force(ribbon)
    if(x$type == "complex") {
      cl <- match.call()
      cl$x <- solist(Re=Re(x), Im=Im(x), Mod=Mod(x), Arg=Arg(x))
      cl[[1]] <- as.name('plot')
      cl$main <- main
      out <- eval(cl, parent.frame())
      return(invisible(out))
    }
    ribside <- match.arg(ribside)
    col.given <- !is.null(col)
    dotargs <- list(...)

    stopifnot(is.list(ribargs))
    user.ticks <- ribargs$at
    user.nint <- ribargs$nint
    user.ribbonlabels <- ribargs$labels
    
    if(!is.null(clipwin)) {
      x <- x[as.rectangle(clipwin)]
      if(!is.rectangle(clipwin)) x <- x[clipwin, drop=FALSE]
    }

    zlim <- dotargs$zlim

    x <- repair.image.xycoords(x)

    xtype <- x$type
    xbox <- as.rectangle(x)
    
    do.log <- identical(log, TRUE)
    if(do.log && !(x$type %in% c("real", "integer")))
      stop(paste("Log transform is undefined for an image of type",
                 sQuote(xtype)))

    ## secret interface for handling log-transformed data
    already.log <- identical(log, "already")
    
    # determine whether pixel values are to be treated as colours
    if(!is.null(valuesAreColours)) {
      # argument given - validate
      stopifnot(is.logical(valuesAreColours))
      if(valuesAreColours) {
        ## pixel values must be factor or character
        if(!xtype %in% c("factor", "character")) {
          if(do.plot)
            warning(paste("Pixel values of type", sQuote(xtype),
                          "are not interpretable as colours"))
          valuesAreColours <- FALSE
        } else if(col.given) {
          ## colour info provided: contradictory
          if(do.plot)
            warning(paste("Pixel values are taken to be colour values,",
                          "because valuesAreColours=TRUE;", 
                          "the colour map (argument col) is ignored"),
                    call.=FALSE)
          col <- NULL
        }
        if(do.log && do.plot) 
          warning(paste("Pixel values are taken to be colour values,",
                        "because valuesAreColours=TRUE;", 
                        "the argument log=TRUE is ignored"),
                  call.=FALSE)
      }
    } else if(col.given) {
      # argument 'col' controls colours
      valuesAreColours <- FALSE
    } else if(spatstat.options("monochrome")) {
      valuesAreColours <- FALSE
    } else {
      ## default : determine whether pixel values are colours
      strings <- switch(xtype,
                        character = { as.vector(x$v) },
                        factor    = { levels(x) },
                        { NULL })
      valuesAreColours <- is.character(strings) && 
      !inherits(try(col2rgb(strings), silent=TRUE), "try-error")
      if(valuesAreColours && do.plot)
        splat("Interpreting pixel values as colours",
              "(valuesAreColours=TRUE)")
    }
    # 
    if(valuesAreColours) {
      # colour-valued images are plotted using the code for factor images
      # with the colour map equal to the levels of the factor
      switch(xtype,
             factor = {
               col <- levels(x)
             },
             character = {
               x <- eval.im(factor(x))
               xtype <- "factor"
               col <- levels(x)
             },
             {
               if(do.plot)
                 warning(paste("Pixel values of type", sQuote(xtype),
                               "are not interpretable as colours"))
             })
      # colours not suitable for ribbon
      ribbon <- FALSE
    } 
    
    # transform pixel values to log scale?
    if(do.log) {
      rx <- range(x, finite=TRUE)
      if(all(rx > 0)) {
        x <- eval.im(log10(x))
      } else {
        if(do.plot && any(rx < 0)) 
          warning(paste("Negative pixel values",
                        "omitted from logarithmic colour map;",
                        "range of values =", prange(rx)),
                  call.=FALSE)
        if(do.plot && !all(rx > 0))
          warning("Zero pixel values omitted from logarithmic colour map",
                  call.=FALSE)
        x <- eval.im(log10orNA(x))
      } 
      xtype <- x$type
      values.are.log <- TRUE
      ## functions 'Log' and 'Exp' are used to determine tick marks and labels
      Log <- log10
      Exp <- TenPower
      if(!is.null(zlim))
        dotargs$zlim <- zlim <- log10(zlim)
    } else if(already.log) {
      values.are.log <- TRUE
      Log <- log10
      Exp <- TenPower
    } else {
      values.are.log <- FALSE
      Log <- Exp <- identity
    }
    compress <- decompress <- identity
    
    imagebreaks <- NULL
#    ribbonvalues <- ribbonbreaks <- NULL
    ribbonvalues <- NULL

    ## NOW DETERMINE THE COLOUR MAP
    colfun <- colmap <- NULL
    if(valuesAreColours) {
      ## pixel values are colours; set of colours was determined earlier
      colmap <- colourmap(col=col, inputs=col)
    } else if(!col.given) {
      ## no colour information given: use default
      colfun <- spatstat.options("image.colfun")
    } else if(inherits(col, "colourmap")) {
      ## Bob's your uncle
      colmap <- col
    } else if(is.function(col)) {
      ## Some kind of function determining a colour map
      if(names(formals(col))[1] == "n") {
        ## function(n) -> colour values
        colfun <- col
      } else {
        ## colour map determined by a rule (e.g. 'beachcolours')
        colmap <- invokeColourmapRule(col, x, zlim=zlim, colargs=colargs)
        if(is.null(colmap))
          stop("Unrecognised syntax for colour function")
      }
    }

    switch(xtype,
           real    = {
             vrange <- numericalRange(x, zlim)
             if(!is.null(colmap)) {
               # explicit colour map
               s <- summary(colmap)
               col <- s$outputs
               if(s$discrete)
                 stop("Discrete colour map is not applicable to real values")
               imagebreaks <- s$breaks
               if(is.function(s$compress)) {
                 ## remember these transformations
                 compress <- s$compress
                 decompress <- s$decompress
                 ## transform pixel values to the compressed scale
                 x <- eval.im(compress(x))
                 imagebreaks <- compress(imagebreaks)
                 values.are.log <- samefunction(compress, log10)
               }
               vrange <- range(imagebreaks)
             }
             trivial <- (diff(vrange) <= zap * .Machine$double.eps)
             #' ribbonvalues: a sequence of pixel values, mapped to colours
             #' ribbonrange:  (min, max) of pixel values mapped by ribbon
             #' nominalrange: range of (scaled) values shown on ribbon 
             #' nominalmarks: (scaled) values shown on ribbon at tick marks
             #' ribbonticks: pixel values corresponding to tick marks 
             #' ribbonlabels: text displayed at tick marks
             #' reusableticks: reusable value of user.ticks
             if(trivial) {
               ribbonvalues <- mean(vrange)
               nominalmarks <- compress(Log(ribscale * Exp(decompress(ribbonvalues))))
             } else {
               ribbonvalues <- seq(from=vrange[1L], to=vrange[2L],
                                   length.out=ribn)
               ribbonrange <- vrange
               nominalrange <- compress(Log(ribscale * Exp(decompress(ribbonrange))))
               nominalmarks <- user.ticks %orifnull% Ticks(nominalrange,
                                                           log=values.are.log,
                                                           nint=user.nint,
                                                           deco=decompress)
             }
             ribbonticks <- compress(Log(nominalmarks/ribscale))
             ribbonlabels <- user.ribbonlabels %orifnull% paste(nominalmarks)
             reusableticks <- nominalmarks
           },
           integer = {
             vrange <- numericalRange(x, zlim)
             if(!is.null(colmap)) {
               # explicit colour map
               s <- summary(colmap)
               col <- s$outputs
               if(s$discrete) {
                 imagebreaks <- c(s$inputs[1] - 0.5, s$inputs + 0.5)
               } else {
                 imagebreaks <- s$breaks
                 if(is.function(s$compress)) {
                   ## remember these transformations
                   compress <- s$compress
                   decompress <- s$decompress
                   ## transform pixel values to the compressed scale
                   x <- eval.im(compress(x))
                   imagebreaks <- compress(imagebreaks)
                   vrange <- range(imagebreaks)
                   values.are.log <- samefunction(compress, log10)
                 }
               }
             } 
             trivial <- (diff(vrange) < sqrt(.Machine$double.eps))
             nominalrange <- Log(ribscale * Exp(vrange))
             if(!is.null(user.ticks)) {
               nominalmarks <- user.ticks
             } else {
               nominalmarks <- Ticks(nominalrange,
                                     log=do.log,
                                     nint = user.nint,
                                     deco = decompress)
               nominalmarks <- nominalmarks[nominalmarks %% 1 == 0]
               nominalmarks <- decompress(nominalmarks)
             }
             reusableticks <- nominalmarks
             ribbonticks <- compress(Log(nominalmarks/ribscale))
             ribbonlabels <- user.ribbonlabels %orifnull% paste(nominalmarks)
             if(!do.log && isTRUE(all.equal(ribbonticks,
                                            vrange[1]:vrange[2]))) {
               #' each possible pixel value will appear in ribbon
               ribbonvalues <- vrange[1]:vrange[2]
               imagebreaks <- c(ribbonvalues - 0.5, vrange[2] + 0.5)
               ribbonrange <- range(imagebreaks)
               ribbonticks <- ribbonvalues
               ribbonlabels <- user.ribbonlabels %orifnull% paste(ribbonticks * ribscale)
             } else {
               ## not all possible values will appear in ribbon
               ribn <- min(ribn, diff(vrange)+1)
               ribbonvalues <- seq(from=vrange[1], to=vrange[2],
                                   length.out=ribn)
               ribbonrange <- vrange
             }
           },
           logical = {
             vrange <- c(0,1)
             trivial <- FALSE
             imagebreaks <- c(-0.5, 0.5, 1.5)
             ribbonvalues <- c(0,1)
             ribbonrange <- range(imagebreaks)
#             ribbonbreaks <- imagebreaks
             ribbonticks <- user.ticks %orifnull% ribbonvalues
             ribbonlabels <- user.ribbonlabels %orifnull% c("FALSE", "TRUE")
             reusableticks <- ribbonticks
             if(!is.null(colmap)) 
               col <- colmap(c(FALSE,TRUE))
           },
           factor  = {
             lev <- levels(x)
             nvalues <- length(lev)
             trivial <- (nvalues < 2)
             # ensure all factor levels plotted separately
             fac <- factor(lev, levels=lev)
             intlev <- as.integer(fac)
             imagebreaks <- c(intlev - 0.5, max(intlev) + 0.5)
             ribbonvalues <- intlev
             ribbonrange <- range(imagebreaks)
#             ribbonbreaks <- imagebreaks
             ribbonticks <- user.ticks %orifnull% ribbonvalues
             ribbonlabels <- user.ribbonlabels %orifnull% paste(lev)
             reusableticks <- ribbonticks
             vrange <- range(intlev)
             if(!is.null(colmap) && !valuesAreColours) 
               col <- colmap(fac)
           },
           character  = {
             x <- eval.im(factor(x))
             lev <- levels(x)
             nvalues <- length(lev)
             trivial <- (nvalues < 2)
             # ensure all factor levels plotted separately
             fac <- factor(lev, levels=lev)
             intlev <- as.integer(fac)
             imagebreaks <- c(intlev - 0.5, max(intlev) + 0.5)
             ribbonvalues <- intlev
             ribbonrange <- range(imagebreaks)
#             ribbonbreaks <- imagebreaks
             ribbonticks <- user.ticks %orifnull% ribbonvalues
             ribbonlabels <- user.ribbonlabels %orifnull% paste(lev)
             reusableticks <- ribbonticks
             vrange <- range(intlev)
             if(!is.null(colmap)) 
               col <- colmap(fac)
           },
           stop(paste("Do not know how to plot image of type", sQuote(xtype)))
           )
  
    ## Compute colour values to be passed to image.default
    if(!is.null(colmap)) {
      ## Explicit colour map object
      colourinfo <- list(breaks=imagebreaks, col=col)
    } else if(!is.null(colfun)) {
      ## Function colfun(n)
      if(trivial) ncolours <- 1
      colourinfo <-
        if(is.null(imagebreaks)) list(col=colfun(ncolours)) else
        list(breaks=imagebreaks, col=colfun(length(imagebreaks) - 1L))
    } else if(col.given) {
      ## Colour values
      if(inherits(try(col2rgb(col), silent=TRUE), "try-error"))
        stop("Unable to interpret argument col as colour values")
      if(is.null(imagebreaks)) {
        colourinfo <- list(col=col)
      } else {
        nintervals <- length(imagebreaks) - 1
        colourinfo <- list(breaks=imagebreaks, col=col)
        if(length(col) != nintervals)
          stop(paste("Length of argument", dQuote("col"),
                     paren(paste(length(col))),
                     "does not match the number of distinct values",
                     paren(paste(nintervals))))
      }
    } else stop("Internal error: unable to determine colour values")

    if(spatstat.options("monochrome")) {
      ## transform to grey scale
      colourinfo$col <- to.grey(colourinfo$col)
    }

    if(isTRUE(reverse.col) && !valuesAreColours) {
      ## reverse the colour sequence (using rev.colourmap or rev.default)
      colourinfo$col <- rev(colourinfo$col)
    }
    
    # colour map to be returned (invisibly)
    i.col <- colourinfo$col
    i.bks <- colourinfo$breaks
    output.colmap <-
      if(is.null(i.col)) NULL else
      if(inherits(i.col, "colourmap")) i.col else
      if(inherits(colmap, "colourmap")) colmap else
      if(valuesAreColours) colourmap(col=i.col, inputs=i.col) else
      switch(xtype,
             integer=,
             real= {
               if(!do.log) {
                 ## linear colour map 
                 if(!is.null(i.bks)) {
                   ## possibly uneven breaks
                   colourmap(col=i.col, breaks=i.bks)
                 } else {
                   colourmap(col=i.col, range=vrange, gamma=gamma)
                 }
               } else {
                 ## logarithmic colour map
                 if(!is.null(i.bks)) {
                   colourmap(col=i.col, breaks=TenPower(i.bks),
                             compress=log10, decompress=TenPower)
                 } else {
                   colourmap(col=i.col, range=TenPower(vrange),
                             gamma=gamma,
                             compress=log10, decompress=TenPower)
                 }
               }
             },
             logical={
               colourmap(col=i.col, inputs=c(FALSE,TRUE))
             },
             character=,
             factor={
               colourmap(col=i.col, inputs=lev)
             },
             NULL)

    ## save tickmark values
    attr(output.colmap, "at") <- reusableticks
    
    ## gamma correction
    soc <- summary(output.colmap)
    if(!is.null(gamma <- soc$gamma) && gamma != 1)
      colourinfo$breaks <- soc$breaks

    ##  ........ decide whether to use rasterImage .........

    if(!isFALSE(useRaster)) {
      ## get device capabilities
      ##      (this will start a graphics device if none is active)
      rasterable <- safeDevCapabilities()$rasterImage
      if(is.null(rasterable)) rasterable <- "no"
      ##
      can.use.raster <-
        switch(rasterable,
               yes=TRUE,
               no=FALSE,
               "non-missing"=!anyNA(x$v),
               FALSE)
      if(is.null(useRaster)) {
        useRaster <- can.use.raster
      } else if(useRaster && !can.use.raster) {
        whinge <- "useRaster=TRUE is not supported by the graphics device"
        if(rasterable == "non-missing")
          whinge <- paste(whinge, "for images with NA values")
        warning(whinge, call.=FALSE)
      }
    }

    ## ........ catch old usage (undocumented ) ................
    contourargs <- resolve.defaults(contourargs, dotargs$args.contour)

    ## ........ background object ..............................
    if(!is.null(background)) {
      if(isTRUE(clip.background)) {
        bkg <- try(background[as.rectangle(x), drop=FALSE], silent=TRUE)
        if(inherits(bkg, "try-error")) {
          warning("Unable to clip the background object", call.=FALSE)
        } else {
          background <- bkg
        }
      }
      backbox <- Frame(background)
    } else {
      backbox <- NULL
    }
    
    ## ........ start plotting .................

    if(!isTRUE(ribbon) || (trivial && isTRUE(drop.ribbon))) {
      ## no ribbon wanted

      attr(output.colmap, "bbox") <- boundingbox(as.rectangle(x), backbox)
      if(!do.plot)
        return(output.colmap)

      ## plot background if specified
      if(!is.null(background)) {
        plot(background, main="")
        add <- TRUE
      }

      ## plot title centred over main image area
      if(show.all && sum(nzchar(main))) {
        mainargnames <- c("cex.main", "adj.main", "col.main")
        mainargs <- dotargs[names(dotargs) %in% mainargnames]
        do.call.plotfun(plot.owin,
                        resolve.defaults(list(x=quote(xbox),
                                              type="n",
                                              main=main,
                                              add=add,
                                              show.all=TRUE),
                                         mainargs,
                                         list(claim.title.space=TRUE)),
                        extrargs=graphicsPars("owin"))
        main <- ""
        add <- TRUE
      }

      ## plot image without ribbon
      image.doit(imagedata=list(x=cellbreaks(x$xcol, x$xstep),
                                y=cellbreaks(x$yrow, x$ystep),
                                z=t(x$v)),
                 ## formal arguments
                 add=add, show.all=show.all,
                 W=xbox,
                 addcontour=addcontour, contourargs=contourargs, 
                 workaround=workaround,
                 ## argument lists 
                 list(axes=FALSE, xlab="",ylab=""),
                 dotargs,
                 list(useRaster=useRaster),
                 colourinfo,
                 list(zlim=vrange, asp = 1, main = main),
                 list(values.are.log=values.are.log))
##      if(add && show.all)
##        fakemaintitle(x, main, dotargs)

      do.box.etc(Frame(x), add, dotargs)
      
      return(invisible(output.colmap))
    }
    
    # determine plot region
    bb <- owinInternalRect(x$xrange, x$yrange)
    Width <- diff(bb$xrange)
    Height <- diff(bb$yrange)
    Size <- max(Width, Height)
    switch(ribside,
           right={
             # ribbon to right of image
             bb.rib <- owinInternalRect(bb$xrange[2] + c(ribsep, ribsep+ribwid) * Size,
                            bb$yrange)
             rib.iside <- 4
           },
           left={
             # ribbon to left of image
             bb.rib <- owinInternalRect(bb$xrange[1] - c(ribsep+ribwid, ribsep) * Size,
                            bb$yrange)
             rib.iside <- 2
           },
           top={
             # ribbon above image
             bb.rib <- owinInternalRect(bb$xrange,
                            bb$yrange[2] + c(ribsep, ribsep+ribwid) * Size)
             rib.iside <- 3
           },
           bottom={
             # ribbon below image
             bb.rib <- owinInternalRect(bb$xrange,
                            bb$yrange[1] - c(ribsep+ribwid, ribsep) * Size)
             rib.iside <- 1
           })
    bb.all <- boundingbox(bb.rib, bb, backbox)

    attr(output.colmap, "bbox") <- bb.all
    attr(output.colmap, "bbox.legend") <- bb.rib
    attr(output.colmap, "side.legend") <- rib.iside
    if(!do.plot)
      return(output.colmap)

    pt <- prepareTitle(main)
    
    if(!add) {
      ## establish coordinate system
      do.call.plotfun(plot.owin,
                      resolve.defaults(list(x=quote(bb.all),
                                            type="n",
                                            main=pt$blank),
                                       dotargs),
                      extrargs=graphicsPars("owin"))
      add <- TRUE
    }
    if(show.all) {
      ## plot title centred over main image area 'bb'
      do.call.plotfun(plot.owin,
                      resolve.defaults(list(x=quote(bb),
                                            type="n",
                                            main=main,
                                            add=add,
                                            show.all=TRUE),
                                       dotargs,
                                       list(claim.title.space=TRUE)),
                      extrargs=graphicsPars("owin"))
      main <- ""
      add <- TRUE
    }
    if(!is.null(background)) {
      ## plot background
      plot(background, add=TRUE)
    }
    # plot image
    image.doit(imagedata=list(x=cellbreaks(x$xcol, x$xstep),
                              y=cellbreaks(x$yrow, x$ystep),
                              z=t(x$v)),
               ## formal arguments               
               add=TRUE, show.all=show.all,
               W=xbox,
               addcontour=addcontour, contourargs=contourargs,
               workaround=workaround,
               ## argument lists
               list(axes=FALSE, xlab="", ylab=""),
               dotargs,
               list(useRaster=useRaster),
               colourinfo,
               list(zlim=vrange, asp = 1, main = main),
               list(values.are.log=values.are.log))

##    if(add && show.all)
##      fakemaintitle(bb.all, main, ...)
    
    # box or axes for image
    do.box.etc(bb, add, dotargs)

    # plot ribbon image containing the range of image values
    rib.npixel <- length(ribbonvalues) + 1
    switch(ribside,
           left=,
           right={
             # vertical ribbon
             rib.xcoords <- bb.rib$xrange
             rib.ycoords <- seq(from=bb.rib$yrange[1],
                                to=bb.rib$yrange[2],
                                length.out=rib.npixel)
             rib.z <- matrix(ribbonvalues, ncol=1)
             rib.useRaster <- useRaster
           },
           top=,
           bottom={
             # horizontal ribbon
             rib.ycoords <- bb.rib$yrange
             rib.xcoords <- seq(from=bb.rib$xrange[1],
                                to=bb.rib$xrange[2],
                                length.out=rib.npixel)
             rib.z <- matrix(ribbonvalues, nrow=1)
             # bug workaround
             rib.useRaster <- FALSE 
           })
    image.doit(imagedata=list(x=rib.xcoords,
                              y=rib.ycoords,
                              z=t(rib.z)),
               ## formal arguments
               add=TRUE, show.all=show.all,
               W=bb.rib,
               addcontour=addcontour, contourargs=contourargs,
               workaround=workaround,
               ## argument lists
               ribargs,
               list(useRaster=rib.useRaster),
               list(main="", sub="", xlab="", ylab=""),
               dotargs,
               colourinfo,
               list(values.are.log=values.are.log))
    # box around ribbon?
    resol <- resolve.defaults(ribargs, dotargs)
    if(!identical(resol$box, FALSE))
      plot(as.owin(bb.rib), add=TRUE)
    # scale axis for ribbon image
    ribaxis <- !(identical(resol$axes, FALSE) || identical(resol$ann, FALSE))
    if(ribaxis) {
      ribaxis.iside <- rib.iside
      ## check for user-supplied xlim, ylim with reverse order
      ll <- resolve.defaults(ribargs, dotargs, list(xlim=NULL, ylim=NULL))
      xlimflip <- is.numeric(ll$xlim) && (diff(ll$xlim) < 0)
      ylimflip <- is.numeric(ll$ylim) && (diff(ll$ylim) < 0)
      if(xlimflip) ribaxis.iside <- c(1, 4, 3, 2)[ribaxis.iside] 
      if(ylimflip) ribaxis.iside <- c(3, 2, 1, 4)[ribaxis.iside]
      ##
      axisargs <- list(side=ribaxis.iside, labels=ribbonlabels)
      switch(ribside,
             right={
               if(trivial) {
                 at <- mean(bb.rib$yrange)
               } else {
                 scal <- diff(bb.rib$yrange)/diff(ribbonrange)
                 at <- bb.rib$yrange[1] + scal * (ribbonticks - ribbonrange[1])
               }
               axisargs <- append(axisargs, list(at=at))
               posargs <- list(pos=bb.rib$xrange[2],
                               yaxp=c(bb.rib$yrange, length(ribbonticks)))
             },
             left={
               if(trivial) {
                 at <- mean(bb.rib$yrange)
               } else {
                 scal <- diff(bb.rib$yrange)/diff(ribbonrange)
                 at <- bb.rib$yrange[1] + scal * (ribbonticks - ribbonrange[1])
               }
               axisargs <- append(axisargs, list(at=at))
               posargs <- list(pos=bb.rib$xrange[1],
                               yaxp=c(bb.rib$yrange, length(ribbonticks)))
             },
             top={
               if(trivial) {
                 at <- mean(bb.rib$xrange)
               } else {
                 scal <- diff(bb.rib$xrange)/diff(ribbonrange)
                 at <- bb.rib$xrange[1] + scal * (ribbonticks - ribbonrange[1])
               }
               axisargs <- append(axisargs, list(at=at))
               posargs <- list(pos=bb.rib$yrange[2],
                               xaxp=c(bb.rib$xrange, length(ribbonticks)))
             },
             bottom={
               if(trivial) {
                 at <- mean(bb.rib$xrange)
               } else {
                 scal <- diff(bb.rib$xrange)/diff(ribbonrange)
                 at <- bb.rib$xrange[1] + scal * (ribbonticks - ribbonrange[1])
               }
               axisargs <- append(axisargs, list(at=at))
               posargs <- list(pos=bb.rib$yrange[1],
                               xaxp=c(bb.rib$xrange, length(ribbonticks)))
             })
      do.call.plotfun(graphics::axis,
                      resolve.defaults(axisargs, ribargs, dotargs, posargs),
                      extrargs=graphicsPars("axis"))
    }
    ## label next to ribbon
    if(!is.null(riblab)) {
      ## assemble arguments
      if(!is.list(riblab)) riblab <- list(text=riblab)
      riblab <- resolve.defaults(riblab, list(side=rib.iside))
      ## temporarily suppress clipping 
      opa <- par(xpd=NA)
      on.exit(par(opa))
      if(sideCode(riblab$side) == sideCode(rib.iside)) {
        ## use 'mtext'
        do.call(mtext, riblab)
      } else {
        ## will use 'text'
        rls <- sideCode(riblab$side, "word")
        ## conform to formal arguments of 'text'
        riblab <- riblab[names(riblab) != "side"]
        names(riblab) <- sub("text", "labels", names(riblab))
        ## determine spatial position arguments for 'text'
        switch(rls,
               bottom = {
                 x0 <- mean(bb.rib$xrange)
                 y0 <- bb.rib$yrange[1L]
                 srt <- 0
               },
               left = {
                 x0 <- bb.rib$xrange[1L]
                 y0 <- mean(bb.rib$yrange)
                 srt <- 90
               },
               top = {
                 x0 <- mean(bb.rib$xrange)
                 y0 <- bb.rib$yrange[2L]
                 srt <- 0
               },
               right = {
                 x0 <- bb.rib$xrange[2L]
                 y0 <- mean(bb.rib$yrange)
                 srt <- -90
               })
        rlpos <- sideCode(rls)
        do.call(text,
                resolve.defaults(riblab,
                                 list(x=x0, y=y0, pos=rlpos, srt=srt)))
      }
    }
    #
    return(invisible(output.colmap))
  }

  PlotIm
})

invokeColourmapRule <- function(colfun, x, ..., zlim=NULL, colargs=list()) {
  ## utility for handling special functions that generate colour maps
  ## either 
  ##        function(... range) -> colourmap
  ##        function(... inputs) -> colourmap
  stopifnot(is.im(x))
  stopifnot(is.function(colfun))
  colargnames <- names(formals(colfun))
  ## Convert it to a 'colourmap'
  colmap <- NULL
  xtype <- x$type
  if(xtype %in% c("real", "integer") && "range" %in% colargnames) {
    ## function(range) -> colourmap
    vrange <- range(range(x, finite=TRUE), zlim)
    cvals <- try(do.call.matched(colfun,
                                 append(list(range=vrange), colargs)),
                 silent=TRUE)
    if(!inherits(cvals, "try-error")) {
      colmap <- if(inherits(cvals, "colourmap")) cvals else
      if(is.character(cvals)) colourmap(cvals, range=vrange) else NULL
    }
  } else if(xtype != "real" && "inputs" %in% colargnames) {
    ## function(inputs) -> colourmap
    vpossible <- switch(xtype,
                        logical = c(FALSE, TRUE),
                        factor = levels(x),
                        unique(as.matrix(x)))
    if(!is.null(vpossible) && length(vpossible) < 256) {
      cvals <- try(do.call.matched(colfun,
                                   append(list(inputs=vpossible),
                                          colargs)),
                   silent=TRUE)
      if(!inherits(cvals, "try-error")) {
        colmap <- if(inherits(cvals, "colourmap")) cvals else
        if(is.character(cvals))
          colourmap(cvals, inputs=vpossible) else NULL
      }
    }
  }
  return(colmap)
}

########################################################################

image.im <- plot.im

######################################################################

contour.im <- function (x, ..., main, axes=FALSE, add=FALSE,
                        nlevels=10, levels=NULL, labels=NULL, log=FALSE, 
                        col=par("fg"), 
                        clipwin=NULL, show.all=!add, do.plot=TRUE)
{
  defaultmain <- if(missing(main)) short.deparse(substitute(x)) else NULL
  dotargs <- list(...)
  bb <- Frame(x)
  xtype <- x$type
  ## contour spacing
  do.log <- isTRUE(log)
  if(do.log && !(xtype %in% c("real", "integer")))
    stop(paste("Log transform is undefined for an image of type",
               sQuote(xtype)))
  ## return value
  result <- bb
  attr(result, "bbox") <- bb
  if(!do.plot) return(result)
  ## main title
  sop <- spatstat.options("par.contour")
  if(missing(main)) 
    main <- resolve.1.default(list(main=defaultmain), sop)
  pt <- prepareTitle(main)
  ## plotting parameters
  if(missing(add)) {
    force(add) ## use default in formal arguments, unless overridden
    add <- resolve.1.default(list(add=add), sop)
  }
  if(missing(axes)) {
    force(axes)
    axes <- resolve.1.default(list(axes=axes), sop)
  }
  axes <- axes && !add
  col0 <- if(inherits(col, "colourmap")) par("fg") else col
  ## clip to subset
  if(!is.null(clipwin))
    x <- x[clipwin, drop=FALSE]
  #' start plotting
  if(!add) {
    ## new plot - establish coordinate system
    if(axes && show.all) {
      #' standard plot initialisation in base graphics
      do.call.plotfun(plot.default,
                      resolve.defaults(
                                       list(x = range(x$xcol),
                                            y = range(x$yrow),
                                            type = "n"),
                                       list(...),
                                       list(asp = 1,
                                            xlab = "x",
                                            ylab = "y",
                                            col = col0,
                                            main = main)))
    } else {
      #' plot invisible bounding box
      do.call.plotfun(plot.owin,
                      resolve.defaults(list(x=quote(bb),
                                            type="n",
                                            main=pt$blank),
                                       dotargs),
                      extrargs=graphicsPars("owin"))
    }
  } 
  if(show.all && !axes) {
    ## plot title centred over contour region
    do.call.plotfun(plot.owin,
                    resolve.defaults(list(x=quote(bb),
                                          main=main,
                                          add=TRUE,
                                          show.all=TRUE),
                                     dotargs,
                                     list(col.main=col0)),
                    extrargs=graphicsPars("owin"))
  }
  ## Determine contour levels
  if(!do.log) {
    ## even spacing
    contourargs <- list(nlevels=nlevels, levels=levels, labels=labels)
  } else {
    ## logarithmic spacing
    rx <- range(x, finite=TRUE)
    ## take logarithm of pixel data
    if(all(rx > 0)) {
      x <- eval.im(log10(x))
    } else {
      if(do.plot && any(rx < 0)) 
        warning(paste("Negative pixel values",
                      "omitted from logarithmic colour map;",
                      "range of values =", prange(rx)),
                call.=FALSE)
      if(do.plot && !all(rx < 0))
        warning("Zero pixel values omitted from logarithmic colour map",
                call.=FALSE)
      x <- eval.im(log10orNA(x))
    }
    ## determine levels
    if(!is.null(levels)) {
      levels <- log10(levels)
      if(is.null(labels)) labels <- paste(levels)
    } else {
      logra <- range(x, finite=TRUE)
      ## default levels commensurate with logarithmic colour scale
      if(diff(logra) > 1.5 && missing(nlevels)) {
        wholepowers <- 10^(floor(logra[1]):ceiling(logra[2]))
        explevels <- sort(as.numeric(outer(wholepowers, c(1,2,5), "*")))
      } else {
        explevels <- pretty(10^logra, nlevels)
      }
      explevels <- explevels[inside.range(explevels, 10^logra)]
      levels <- log10(explevels)
      if(is.null(labels)) labels <- paste(explevels)
    }
    contourargs <- list(levels=levels, labels=labels)
  }
  #' plot contour lines
  xcol <- x$xcol
  yrow <- x$yrow
  zmat <- t(x$v)
  dont.complain.about(xcol, yrow, zmat)
  if(!inherits(col, "colourmap")) {
    do.call.plotfun(contour.default,
                    resolve.defaults(list(x=quote(xcol), 
					  y=quote(yrow), 
					  z=quote(zmat),
                                          add=TRUE,
                                          col=col),
                                     contourargs,
                                     list(...),
                                     .MatchNull=FALSE,
                                     .StripNULL=TRUE))
  } else {
    clin <- do.call.matched(contourLines,
                            resolve.defaults(list(x=quote(xcol), 
                                                  y=quote(yrow), 
                                                  z=quote(zmat)),
                                             contourargs,
                                             list(...),
                                             .MatchNull=FALSE,
                                             .StripNULL=TRUE))
    linpar <- graphicsPars("lines")
    for(i in seq_along(clin)) {
      lini <- clin[[i]]
      levi <- lini$level
      coli <- col(levi)
      argi <- resolve.defaults(lini[c("x", "y")],
                               list(...),
                               list(col=coli))
      do.call.matched(lines.default, argi, extrargs=linpar)
    }
  }
  return(invisible(result))
}


## not exported:

TenPower <- function(x) { 10^x }

log10orNA <- function(x) {
  y <- rep(NA_real_, length(x))
  ok <- !is.na(x) & (x > 0)
  y[ok] <- log10(x[ok])
  return(y)
}
