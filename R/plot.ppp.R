#
#	plot.ppp.R
#
#	$Revision: 1.124 $	$Date: 2025/04/24 03:33:11 $
#
#
#--------------------------------------------------------------------------

plot.ppp <- function(x, main, ..., clipwin=NULL,
                     chars=NULL, cols=NULL, use.marks=TRUE,
                     which.marks=NULL, add=FALSE, type=c("p", "n"), 
                     legend=TRUE, leg.side=c("left", "bottom", "top", "right"),
                     leg.args=list(),
                     symap=NULL, maxsize=NULL, meansize=NULL, markscale=NULL,
                     minsize=NULL, zerosize=NULL, zap=0.01, 
                     show.window=show.all, show.all=!add, do.plot=TRUE,
                     multiplot=TRUE)
{
  if(missing(main))
    main <- short.deparse(substitute(x))

  type <- match.arg(type)
  if(missing(legend)) legend <- (type == "p")

  if(clipped <- !is.null(clipwin)) {
    stopifnot(is.owin(clipwin))
    W <- Window(x)
    clippy <- if(is.mask(W)) intersect.owin(W, clipwin) else edges(W)[clipwin]
    x <- x[clipwin]
  } else clippy <- NULL
  
  ## sensible default position
  legend <- legend && show.all
  if(legend) {
    leg.side <- match.arg(leg.side)
    vertical <- (leg.side %in% c("left", "right"))
  }
  
  ## ................................................................
  ## Handle multiple columns of marks as separate plots
  ##  (unless add=TRUE or which.marks selects a single column
  ##   or multipage = FALSE)
  if(use.marks && is.data.frame(mx <- marks(x))) {
    implied.all <- is.null(which.marks)
    want.several <- implied.all || is.data.frame(mx <- mx[,which.marks])
    do.several <- want.several && !add && multiplot
    if(do.several) {
      ## generate one plot for each column of marks
      y <- solapply(mx, setmarks, x=x)
      out <- do.call(plot,
                     resolve.defaults(list(x=quote(y), main=main,
                                           show.window=show.window && !clipped,
                                           do.plot=do.plot,
                                           type=type,
                                           symap=symap),
                                      list(...),
                                      list(equal.scales=TRUE),
                                      list(panel.end=clippy),
                                      list(legend=legend,
                                           leg.side=leg.side,
                                           leg.args=leg.args),
                                      list(chars=chars, cols=cols,
                                           maxsize=maxsize,
                                           meansize=meansize,
                                           markscale=markscale,
                                           minsize=minsize,
                                           zerosize=zerosize,
                                           zap=zap)))
      return(invisible(out))
    } 
    if(is.null(which.marks)) {
      which.marks <- 1
      if(do.plot) message("Plotting the first column of marks")
    }
  }
  
  ## ............... unmarked, or single column of marks ....................

  ## Determine symbol map and mark values to be used
  y <- x
  if(!is.marked(x, na.action="ignore") || !use.marks) {
    ## Marks are not mapped.
    marx <- NULL
    if(is.null(symap))
      symap <- default.symbolmap(unmark(x), ..., chars=chars, cols=cols)
  } else {
    ## Marked point pattern
    marx <- marks(y, dfok=TRUE)
    if(is.data.frame(marx)) {
      ## select column or take first colum
      marx <- marx[, which.marks]
      y <- setmarks(y, marx)
    }
    if(npoints(y) > 0) {
      ok <- complete.cases(as.data.frame(y))
      if(!any(ok)) {
        warning("All mark values are NA; plotting locations only.")
        if(is.null(symap))
          symap <- default.symbolmap(unmark(x), ..., chars=chars, cols=cols)
      } else if(any(!ok)) {
        warning(paste("Some marks are NA;",
                      "corresponding points are omitted."))
        x <- x[ok]
        y <- y[ok]
        marx <- marks(y)
      }
    }
    ## apply default symbol map
    if(is.null(symap))
      symap <- default.symbolmap(y, chars=chars, cols=cols, 
                                 maxsize=maxsize, meansize=meansize,
                                 markscale=markscale,
                                 minsize=minsize, zerosize=zerosize,
                                 ...)
  }

  ## Determine bounding box for main plot
  BB <- as.rectangle(x)
  sick <- inherits(x, "ppp") && !is.null(rejects <- attr(x, "rejects"))
  if(sick) {
    ## Get relevant parameters
    par.direct <- list(main=main, use.marks=use.marks,
                       maxsize=maxsize, meansize=meansize, markscale=markscale,
                       minsize=minsize, zerosize=zerosize)
    par.rejects <- resolve.1.default(list(par.rejects=list(pch="+")),
                                     list(...))
    par.all <- resolve.defaults(par.rejects, par.direct)
    rw <- resolve.defaults(list(...), list(rejectwindow=NULL))$rejectwindow
    ## determine window for rejects
    rwin <-
      if(is.null(rw))
        rejects$window
      else if(is.logical(rw) && rw)
        rejects$window
      else if(inherits(rw, "owin"))
        rw
      else if(is.character(rw)) {
        switch(rw,
               box={boundingbox(rejects, x)},
               ripras={ripras(c(rejects$x, x$x), c(rejects$y, x$y))},
               stop(paste("Unrecognised option: rejectwindow=", rw)))
      } else stop("Unrecognised format for rejectwindow")
    if(is.null(rwin))
      stop("Selected window for rejects pattern is NULL")
    BB <- boundingbox(BB, as.rectangle(rwin))
  }

  ## Augment bounding box with space for legend, if appropriate
  legend <- legend && (symbolmaptype(symap) != "constant") 
  if(legend) {
    leg.args <- append(list(side=leg.side, vertical=vertical), leg.args)
    if(isTRUE(leg.args$colour.only)) {
      ## only the colour map will be plotted
      ## use layout similar to plot.im
      sizeguess <- NULL
      leg.args <- resolve.defaults(leg.args, list(sep.frac=0.15,
                                                  size.frac=0.05,
                                                  las=1))
    } else {
      ## symbols will be plotted
      ## guess maximum size of symbols
      maxsize <- invoke.symbolmap(symap, symbolmapdomain(symap),
                                  corners(as.rectangle(x)),
                                  add=add, do.plot=FALSE)
      sizeguess <- if(maxsize > 0) (1.5 * maxsize) else NULL
    }
    ## draw up layout
    layoutboxes <- do.call.matched(plan.legend.layout,
                                  append(list(B=quote(BB), size = sizeguess,
                                              started=FALSE, map=symap),
                                         leg.args))
    ## bounding box for everything
    BB <- layoutboxes[["A"]]
    ## bounding box for legend
    legbox <- layoutboxes[["b"]]
    attr(symap, "legbox") <- legbox
  }

  ## return now if not plotting
  attr(symap, "bbox") <- BB
  if(!do.plot)
    return(invisible(symap))
    
  ## ............. start plotting .......................
  pt <- prepareTitle(main)
  main <- pt$main
  nlines <- pt$nlines
  blankmain <- if(nlines == 0) "" else rep("  ", nlines)
  dflt <- list(cex.main=1, xlim=NULL, ylim=NULL,
               ann=FALSE, axes=FALSE, xlab="", ylab="")
  rez <- resolve.defaults(list(...), dflt)[names(dflt)]
  do.call(plot.owin,
          append(list(x=quote(BB), type="n", add=add,
                      main=blankmain, show.all=show.all),
                 rez))
  if(sick) {
    if(show.window) {
      ## plot windows
      if(!is.null(rw)) {
        ## plot window for rejects
        rwinpardefault <- list(lty=2,lwd=1,border=1)
        rwinpars <-
          resolve.defaults(par.rejects, rwinpardefault)[names(rwinpardefault)]
        dont.complain.about(rwin)
        do.call(plot.owin, append(list(quote(rwin), add=TRUE), rwinpars))
      }
      ## plot window of main pattern
      if(!clipped) {
        xwindow <- x$window
        dont.complain.about(xwindow)
        do.call(plot.owin,
                resolve.defaults(list(quote(xwindow), add=TRUE),
                                 list(...),
                                 list(invert=TRUE)))
      } else plot(clippy, add=TRUE, ...)
    }
    if(type != "n") {
      ## plot reject points
      do.call(plot.ppp, append(list(quote(rejects), add=TRUE), par.all))
      warning(paste(rejects$n, "illegal points also plotted"))
    }
    ## the rest is added
    add <- TRUE
  }

  ## Now convert to bona fide point pattern
  x <- as.ppp(x)
  xwindow <- x$window

  ## Plot observation window (or at least the main title)
  dont.complain.about(xwindow)
  do.call(plot.owin,
          resolve.defaults(list(x=quote(xwindow),
                                add=TRUE,
                                main=main,
                                type=if(show.window && !clipped) "w" else "n",
                                show.all=show.all),
                           list(...),
                           list(invert=TRUE)))
  ## If clipped, plot visible part of original window
  if(show.window && clipped)
    plot(clippy, add=TRUE, ...)
  # else if(show.all) fakemaintitle(as.rectangle(xwindow), main, ...)

  if(type != "n") {
    ## plot symbols ##
    invoke.symbolmap(symap, marx, x, add=TRUE)
  }
  
  ## add legend
  if(legend) {
    legendmap <- if(length(leg.args) == 0) symap else 
                 do.call(update, append(list(object=quote(symap)), leg.args))
    dont.complain.about(legendmap)
    do.call(plot.symbolmap,
            append(list(x=quote(legendmap), main="", add=TRUE,
                        xlim=legbox$xrange, ylim=legbox$yrange),
                   leg.args))
  }
  return(invisible(symap))
}

## determine symbol map for marks of points
default.symbolmap.ppp <- local({

  default.symbolmap.ppp <- function(x, ..., 
                                    chars=NULL, cols=NULL, 
                                    fixsize=FALSE,
                                    maxsize=NULL, meansize=NULL, markscale=NULL,
                                    minsize=NULL, zerosize=NULL,
                                    transform=NULL) {
    Y <- lapply(unstack(x),
                dsmEngine,
                ...,
                chars=chars,
                cols=cols,
                fixsize=fixsize,
                maxsize=maxsize,
                meansize=meansize,
                markscale=markscale,
                minsize=minsize,
                zerosize=zerosize,
                transform=transform)
    if(length(Y) == 1)
      Y <- Y[[1L]]
    return(Y)
  }

  ## full argument list
  dsmEngine <- function(x, ..., 
                        chars=NULL,
                        cols=NULL,
                        col=NULL,
                        fixsize=FALSE,
                        maxsize=NULL,
                        meansize=NULL,
                        markscale=NULL,
                        minsize=NULL,
                        zerosize=NULL,
                        markrange=NULL,
                        marklevels=NULL,
                        transform=NULL) {
    marx <- marks(x)
    if(is.null(marx) || npoints(x) == 0) {
      ## null or constant symbol map
      ## consider using transparent colours
      if(is.null(cols) && is.null(col) && 
         !any(c("fg", "bg") %in% names(list(...))) &&
         (nx <- npoints(x)) > 100 &&
         spatstat.options("transparent") &&
         isTRUE(safeDevCapabilities()$semiTransparency))
        cols <- rgb(0,0,0, default.transparency(nx))
      if(!is.null(cols) && !is.null(col)) col <- NULL
      symap <- symbolmap(..., chars=chars, cols=cols, col=col)
      pnames <- symbolmapparnames(symap)
      if("shape" %in% pnames && !("size" %in% pnames)) {
        ## symbols require a size parameter
        m <- symbol.sizes.default(rep(1, max(1, npoints(x))), Window(x),
                                  maxsize=maxsize, meansize=meansize,
                                  minsize=minsize,
                                  zerosize=zerosize)
        symap <- update(symap, size=m)
      }
      return(symap)
    }
    if(!is.null(dim(marx)))
      stop("Internal error: multivariate marks in dsmEngine")

    ## understand user's wishes
    argnames <- names(list(...))
    shapegiven <- "shape" %in% argnames
    chargiven <- (!is.null(chars)) || ("pch" %in% argnames)
    sizegiven <- ("size" %in% argnames) ||
                 (("cex" %in% argnames) && !shapegiven)
    assumecircles <- !(shapegiven || chargiven)
    sizeconstrained <- !all(sapply(list(maxsize, minsize, meansize, zerosize),
                                   is.null))

    ## set defaults
    shapedefault <- if(!assumecircles) list() else list(shape="circles")

    ## potential range of mark values
    if(is.factor(marx)) {
      if(is.null(marklevels)) marklevels <- levels(marx)
    } else {
      if(is.null(markrange)) markrange <- range(marx, na.rm=TRUE, finite=TRUE)
    }
    
    ## pre-transformation of mark values
    if(!is.null(transform)) stopifnot(is.function(transform))
    transforming <- is.function(transform) && !fixsize
    if(transforming) {
      Tmarx <- transform(marx)
    } else {
      Tmarx <- marx
      transform <- NULL
    }
    if(transforming && (is.factor(marx) || is.factor(Tmarx)))
      stop("Sorry, transformations are not yet implemented for factors")
    if(is.factor(Tmarx)) {
      Tmarklevels <- levels(Tmarx)
    } else {
      Tmarkrange <- if(is.null(markrange)) NULL else transform(markrange)
      Tmarkrange <- range(Tmarx, Tmarkrange, na.rm=TRUE, finite=TRUE)
    }

    if(inherits(Tmarx, c("Date", "POSIXt"))) {
      ## ......... transformed marks are dates or date/times ..............
      if(sizegiven) {
        g <- do.call(symbolmap,
                     resolve.defaults(list(range=markrange,
                                           transform=transform),
                                      list(...),
                                      shapedefault,
                                      list(chars=chars, cols=cols)))
        return(g)
      }
      ## attempt to determine a scale for the marks
      Timerange <- range(Tmarx, Tmarkrange, na.rm=TRUE, finite=TRUE)
      y <- scaletointerval(Tmarx, 0, 1, Timerange)
      y <- y[is.finite(y)]
      if(length(y) == 0) return(symbolmap(..., chars=chars, cols=cols))
      scal <- mark.scale.default(y, as.owin(x), markrange=c(0,1),
                                 markscale=markscale, maxsize=maxsize,
                                 meansize=meansize, 
                                 characters=chargiven)
      if(is.na(scal)) return(symbolmap(..., chars=chars, cols=cols))
      ## scale determined
      sizefun <- function(x, scal=1, Timerange=NULL) {
        (scal/2) * scaletointerval(x, 0, 1, Timerange)
      }
      formals(sizefun)[[2]] <- scal  ## ensures value of 'scal' is printed
      formals(sizefun)[[3]] <- Timerange
      ##
      g <- do.call(symbolmap,
                   resolve.defaults(list(range=markrange,
                                         transform=transform),
                                    list(...),
                                    shapedefault,
                                    list(size=sizefun)))
      return(g)
    }
    if(is.numeric(Tmarx)) {
      ## ............. marks are numeric values ...................
      Tmarx <- Tmarx[is.finite(Tmarx)]
      if(length(Tmarx) == 0)
        return(symbolmap(..., chars=chars, cols=cols))
      Tmarkrange <- range(Tmarx, Tmarkrange, na.rm=TRUE, finite=TRUE)
      ## 
      if(sizegiven) {
        ## size function is given
        g <- do.call(symbolmap,
                     resolve.defaults(list(range=markrange,
                                           transform=transform),
                                      list(...),
                                      shapedefault,
                                      list(chars=chars, cols=cols)))
        return(g)
      } else if(fixsize) {
        ## require symbols of equal size
        ## determine fixed physical size
        if(!is.null(meansize)) {
          size <- meansize
        } else if(!is.null(minsize) && !is.null(maxsize)) {
          size <- (minsize+maxsize)/2
        } else if(!is.null(minsize)) {
          size <- minsize
        } else if(!is.null(maxsize)) {
          size <- maxsize
        } else if(!is.null(zerosize) && zerosize > 0) {
          size <- zerosize
        } else {
          ## choose suitable size
          bb <- Frame(x)
          nn <- nndist(x)
          nn <- nn[nn > 0]
          size1 <- 1.4/sqrt(pi * length(marx)/area(bb))
          size2 <- 0.07 * diameter(bb)
          size3 <- if(length(nn)) median(nn) else Inf
          size <- min(size1, size2, size3)
        }
        g <- do.call(symbolmap,
                     resolve.defaults(list(range=markrange,
                                           transform=transform),
                                      list(...),
                                      list(size=size, chars=chars, cols=cols),
                                      shapedefault))
        return(g)
      }
      ## attempt to determine a scale for the (transformed) marks 
      ## degenerate?
      if(all(Tmarkrange == 0))
        return(symbolmap(..., chars=chars, cols=cols))
      ## try scaling
      scal <- mark.scale.default(Tmarx, as.owin(x),
                                 markrange=Tmarkrange,
                                 markscale=markscale, maxsize=maxsize,
                                 meansize=meansize,
                                 minsize=minsize, zerosize=zerosize,
                                 characters=chargiven)
      if(is.na(scal)) return(symbolmap(..., chars=chars, cols=cols))
      ## scale determined
      zerosize <- attr(scal, "zerosize") %orifnull% 0
      scal <- as.numeric(scal)
      if(Tmarkrange[1] >= 0) {
        ## all (transformed) marks are nonnegative
        cexfun <- function(x, scal=1, zerosize=0) { zerosize + scal * x }
        circfun <- function(x, scal=1, zerosize=0) { zerosize + scal * x }
        formals(cexfun)[[2]] <- formals(circfun)[[2]] <- scal
        formals(cexfun)[[3]] <- formals(circfun)[[3]] <- zerosize
        sizedefault <-
          if(sizegiven) list() else
          if(chargiven) list(cex=cexfun) else list(size=circfun)
      } else {
        ## some marks are negative
        shapedefault <-
          if(!assumecircles) list() else
          list(shape=function(x) { ifelse(x >= 0, "circles", "squares") })
        cexfun <- function(x, scal=1, zerosize=0) { zerosize + scal * abs(x) }
        circfun <- function(x, scal=1, zerosize=0) { zerosize + scal * abs(x) }
        formals(cexfun)[[2]] <- formals(circfun)[[2]] <- scal
        formals(cexfun)[[3]] <- formals(circfun)[[3]] <- zerosize
        sizedefault <-
          if(sizegiven) list() else
          if(chargiven) list(cex=cexfun) else list(size=circfun)
      }
      g <- do.call(symbolmap,
                   resolve.defaults(list(range=markrange,
                                         transform=transform),
                                    list(...),
                                    shapedefault,
                                    sizedefault,
                                    list(chars=chars, cols=cols)))
      return(g)
    }
    ##  ...........  non-numeric marks .........................
    if(transforming)
      stop(paste("Argument", sQuote("transform"),
                 "is not yet supported for non-numeric marks"),
           call.=FALSE)
    um <- marklevels %orifnull%
          if(is.factor(marx)) levels(marx) else sortunique(marx)
    ntypes <- length(um)
    if(!is.null(cols))
      cols <- rep.int(cols, ntypes)[1:ntypes]
    if(shapegiven && sizegiven) {
      #' values mapped to symbols (shape and size specified)
      g <- symbolmap(inputs=um, ..., cols=cols)
    } else if(!shapegiven) {
      #' values mapped to 'pch'
      chars <- default.charmap(ntypes, chars)
      g <- symbolmap(inputs=um, ..., chars=chars, cols=cols)
    } else {
      #' values mapped to symbols of equal size
      #' determine size
      scal <- symbol.sizes.default(rep(1, npoints(x)),
                                   Window(x), 
                                   maxsize=maxsize,
                                   meansize=meansize,
                                   minsize=minsize,
                                   zerosize=zerosize,
                                   characters=FALSE)
      g <- symbolmap(inputs=um, ..., size=scal, cols=cols)
    }
    return(g)
  }
                                  
  default.charmap <- function(n, ch=NULL) {
    if(!is.null(ch))
      return(rep.int(ch, n)[1:n])
    if(n <= 25)
      return(1:n)
    ltr <- c(letters, LETTERS)
    if(n <= 52)
      return(ltr[1:n])
    ## wrapped sequence of letters
    warning("Too many types to display every type as a different character")
    return(ltr[1 + (0:(n - 1) %% 52)])
  }

  default.transparency <- function(n) {
    if(n <= 100) 1 else (0.2 + 0.8 * exp(-(n-100)/1000))
  }
  
  default.symbolmap.ppp
})

## utility function to determine mark scale 
## (factor converting mark values to physical sizes on the plot)
## using a default rule

mark.scale.default <- function(marx, w, ...,
                               markrange=NULL,
                               markscale=NULL,
                               maxsize=NULL, meansize=NULL,
                               minsize=NULL, zerosize=NULL,
                               characters=FALSE) {
  ## establish values of parameters markscale, maxsize, meansize
  ngiven <- (!is.null(markscale)) +
            (!is.null(maxsize)) +
            (!is.null(meansize))
  if(ngiven > 1)
     stop("Only one of the arguments markscale, maxsize, meansize",
          " should be given", call.=FALSE)
  if(ngiven == 0) {
    ## if ALL are absent, enforce the spatstat defaults
    ## (which could also be null)
    pop <- spatstat.options("par.points")
    markscale <- pop$markscale
    maxsize   <- pop$maxsize
    meansize  <- pop$meansize
  }
  ng <- (!is.null(minsize)) + (!is.null(zerosize))
  if(ng > 1)
    stop("Arguments minsize and zerosize are incompatible", call.=FALSE)
  if(ng == 0) {
    pop <- spatstat.options("par.points")
    minsize <- pop$minsize
    zerosize <- pop$zerosize
    if(is.null(minsize) && is.null(zerosize))
      zerosize <- 0
  }
  if(!is.null(minsize)) stopifnot(minsize >= 0)
  ## determine range of absolute values of marks to be mapped
  absmarx <- abs(marx)
  ra <- range(absmarx)
  if(!is.null(markrange)) {
    check.range(markrange)
    ra <- range(ra, abs(markrange))
    if(inside.range(0, markrange))
      ra <- range(0, ra)
  }
  minabs <- ra[1L]
  maxabs <- ra[2L]
  ## determine linear map
  ## physical size = zerosize + scal * markvalue
  if(!is.null(markscale)) {
    ## mark scale is already given
    stopifnot(markscale > 0)
    scal <- markscale
    if(!is.null(minsize)) {
      ## required minimum physical size (of marks in range) is specified
      ## determine intercept 'zerosize'
      zerosize <- minsize - scal * ra[1L]
    } ## otherwise 'zerosize' is given or defaults to 0
  } else {
    ## mark scale is to be determined from desired maximum/mean physical size
    if(!is.null(maxsize)) {
      stopifnot(maxsize > 0)
    } else if(!is.null(meansize)) {
      stopifnot(meansize > 0)
    } else {
      ## No prescriptions specified.
      ## Compute default value of 'maxsize'
      ## First guess appropriate max physical size of symbols
      bb <- as.rectangle(w)
      maxradius <- 1.4/sqrt(pi * length(marx)/area(bb))
      maxsize <- 2 * min(maxradius, diameter(bb) * 0.07)
    }
    ## Examine mark values
    epsilon <- 4 * .Machine$double.eps
    if(maxabs < epsilon)
      return(NA)
    
    ## finally determine physical scale for symbols
    if(!is.null(maxsize)) {
      ## required maximum physical size (of marks in range) is specified
      if(!is.null(minsize)) {
        ## required minimum physical size (of marks in range) is specified
        ## map [minabs, maxabs] -> [minsize, maxsize]
        dv <- maxabs - minabs
        if(dv < epsilon) return(NA)
        scal <- (maxsize-minsize)/dv
        zerosize <- minsize - scal * minabs
      } else {
        ## minimum physical size not specified
        ## map [0, maxabs] to [zerosize, maxsize]
        ds <- maxsize - zerosize
        if(ds < epsilon) return(NA)
        scal <- ds/maxabs
        ## check minimum physical size is nonnegative
        if(zerosize + scal * minabs < 0)
          return(NA)
      }
    } else if(!is.null(meansize)) {
      ## required mean physical size (of marks in range) is specified
      meanabs <- mean(if(is.null(markrange)) absmarx else abs(markrange))
      if(!is.null(minsize)) {
        ## required minimum physical size (of marks in range) is specified
        ## map {minabs, meanabs} -> {minsize, meansize}
        dm <- meanabs - minabs
        if(dm < epsilon) return(NA)
        scal <- (meansize-minsize)/dm
        zerosize <- minsize - scal * minabs
      } else {
        ## minimum physical size not specified
        ## map {0, meanabs} -> {zerosize, meansize}
        ds <- meansize - zerosize
        if(ds < epsilon || meanabs < epsilon) return(NA)
        scal <- ds/meanabs
        ## check minimum physical size is nonnegative
        if(zerosize + scal * minabs < 0)
          return(NA)
      }
    } else stop("internal error - neither maxsize nor meansize determined")

    if(characters) {
      ## when using characters ('pch') we need to
      ## convert physical sizes to 'cex' values
      charsize <- max(sidelengths(as.rectangle(w)))/40
      scal <- scal/charsize
      zerosize <- zerosize/charsize
    }
  }

  attr(scal, "zerosize") <- zerosize
  
  return(scal)
}

## utility function to determine symbol sizes using default rule

symbol.sizes.default <- function(markvalues, ...) {
  scal <- mark.scale.default(markvalues, ...)
  if(is.na(scal)) return(NA)
  zerosize <- attr(scal, "zerosize") %orifnull% 0
  scal <- as.numeric(scal)
  sizes <- zerosize + scal * markvalues
  return(sizes)
}

fakemaintitle <- function(bb, main, ...) {
  ## Try to imitate effect of 'title(main=main)' above a specified box
  if(!any(nzchar(main))) return(invisible(NULL))
  bb <- as.rectangle(bb)
  x0 <- mean(bb$xrange)
  y0 <- bb$yrange[2] + length(main) * diff(bb$yrange)/12
  parnames <- c('cex.main', 'col.main', 'font.main')
  parlist <- par(parnames)
  parlist <- resolve.defaults(list(...), parlist)[parnames]
  names(parlist) <- c('cex', 'col', 'font')
  do.call.matched(text.default,
                  resolve.defaults(list(x=x0, y=y0, labels=main),
                                   parlist,    list(...)),
                  funargs=graphicsPars("text"))
  return(invisible(NULL))
}

text.ppp <- function(x, ...) {
  graphics::text.default(x=x$x, y=x$y, ...)
}
