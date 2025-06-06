##
## symbolmap.R
##
##   $Revision: 1.68 $  $Date: 2025/04/26 09:57:55 $
##

symbolmap <- local({

  known.unknowns <- c("shape", "pch", "chars",
                      "size", "cex",
                      "direction", "arrowtype", "headlength", "headangle", 
                      "col", "cols", "fg", "bg",
                      "lty", "lwd", "border", "fill",
                      "etch")

  trycolourmap <- function(...) {
    try(colourmap(...), silent=TRUE)
  }

  symbolmap <- function(..., range=NULL, inputs=NULL,
                        transform=NULL, compress=transform, decompress=NULL) {
    if(!is.null(range) && !is.null(inputs))
      stop("Arguments range and inputs are incompatible")
    ## graphics parameters
    parlist <- list(...)
    ## remove unrecognised parameters and NULL values 
    if(length(parlist) > 0) {
      ok <- names(parlist) %in% known.unknowns
      ok <- ok & !unlist(lapply(parlist, is.null))
      parlist <- parlist[ok]
    }
    got.pars <- (length(parlist) > 0)
    parnames <- names(parlist)
    type <- if(is.null(inputs) && is.null(range)) "constant" else
            if(!is.null(inputs)) "discrete" else "continuous"
    ## transformation applied to input data
    if(!is.null(transform)) {
      if(type != "continuous") {
        warning(paste("Argument 'transform' is ignored for",
                      type, "symbol maps"), call.=FALSE)
        transform <- NULL
      } else stopifnot(is.function(transform))
    }
    Transform <- transform %orifnull% identity
    ## 
    if(got.pars) {
      ## validate parameters
      if(is.null(parnames) || !all(nzchar(parnames)))
        stop("All graphics parameters must have names")
          atomic <- unlist(lapply(parlist, is.atomic))
      functions <- unlist(lapply(parlist, is.function))
      lenfs <- lengths(parlist)
      constants <- atomic & (lenfs == 1)
      if(any(bad <- !(constants | functions))) {
        if(type == "discrete" && any(repairable <- atomic & bad)) {
          ## recycle data to desired length
          parlist[repairable] <- lapply(parlist[repairable],
                                        reptolength,
                                        n=length(inputs))
          bad[repairable] <- FALSE
        }
        if(type == "continuous") {
          ## look for vectors of colour values
          iscol <- bad & sapply(parlist, is.colour) &
            (names(parlist) %in% c("cols", "col", "fg", "bg"))
          ## convert colour values to colour map
          if(any(iscol)) {
            Trange <- Transform(range)
            cmap <- lapply(parlist[iscol], trycolourmap, range=Trange)
            success <- sapply(cmap, inherits, what="colourmap")
            iscol[iscol] <- success
            if(any(iscol)) {
              parlist[iscol] <- cmap[success]
              bad[iscol] <- FALSE
              functions[iscol] <- TRUE
            }
          }
        }
        nbad <- sum(bad)
        if(nbad > 0) 
          stop(paste(ngettext(nbad, "Argument", "Arguments"),
                     commasep(sQuote(parnames[bad])),
                     ngettext(nbad, "is neither a function nor a constant",
                              "are neither functions nor constants")))
      }
      if(type == "constant" && any(functions))
        type <- "continuous"
    }
    switch(type,
           constant ={
             ## set of constant graphics parameters defining a single symbol
             stuff <- list(type=type, parlist=parlist)
             ConstantValue <- as.data.frame(parlist, stringsAsFactors=FALSE)
             f <- function(x) ConstantValue
           },
           discrete = {
             ## finite set of inputs mapped to symbols
             stuff <- list(type=type, inputs=inputs, parlist=parlist)
             f <- function(x) ApplyDiscreteSymbolMap(x, stuff)
           },
           continuous = {
             got.shape <- "shape" %in% parnames
             got.size <- "size" %in% parnames
             got.cha <- any(c("pch", "chars") %in% parnames)
             ## interval of real line (etc) mapped to symbols or characters
             if(!got.cha) {
               ## mapped to symbols
               if(!got.shape)
                 parlist$shape <- "circles"
               if(!got.size)
                 stop("Parameter 'size' is missing")
             }
             rangetype <- if(is.null(range)) "numeric" else
                          if(inherits(range, "POSIXt")) "datetime" else
                          if(inherits(range, "Date")) "date" else
                          if(is.numeric(range)) "numeric" else "unknown"
             stuff <- list(type=type, range=range, rangetype=rangetype,
                           parlist=parlist)
             f <- function(x) ApplyContinuousSymbolMap(x, stuff)
           })
    ## nonlinear transformation for display
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
    stuff <- append(stuff,
                    list(transform=transform,
                         compress=compress,
                         decompress=decompress))
    attr(f, "stuff") <- stuff
    class(f) <- c("symbolmap", class(f))
    f
  }

  reptolength <- function(z, n) { rep.int(z, n)[1:n] }
  
  MapDiscrete <- function(f, x, i) {
    if(is.function(f)) f(x) else if(length(f) == 1) rep.int(f, length(x)) else f[i]
  }
  
  MapContinuous <- function(f, x) {
    if(is.function(f)) f(x) else rep.int(f, length(x))
  }

  ApplyContinuousSymbolMap <- function(x, stuff) {
    with(stuff, {
      if(is.function(transform)) x <- transform(x)
      y <- as.data.frame(lapply(parlist, MapContinuous, x=x),
                         stringsAsFactors=FALSE)
      return(y)
    })
  }
  
  ApplyDiscreteSymbolMap <- function(x, stuff) {
    with(stuff, {
      ii <- match(x, inputs)
      if(anyNA(ii))
        stop("Some values do not belong to the domain of the symbol map")
      y <- as.data.frame(lapply(parlist, MapDiscrete, x=x, i=ii),
                         stringsAsFactors=FALSE)
      return(y)
    })
  }
  symbolmap
})

symbolmaptype <- function(x) { attr(x, "stuff")$type }

symbolmapdomain <- function(x) {
  stuff <- attr(x, "stuff")
  d <- switch(stuff$type,
              constant = { integer(0) },
              discrete = { stuff$inputs },
              continuous = { stuff$range })
  return(d)
}

symbolmapparnames <- function(x) { names(attr(x, "stuff")[["parlist"]]) }

update.symbolmap <- function(object, ...) {
  y <- attr(object, "stuff")
  gpars <- y[["parlist"]]
  oldargs <- y[names(y) %in% names(formals(symbolmap))]
  do.call(symbolmap, resolve.defaults(list(...), append(oldargs, gpars)))
}

print.symbolmap <- function(x, ...) {
  with(attr(x, "stuff"), {
    switch(type,
           constant = {
             if(length(parlist) == 0) {
               cat("Symbol map", "with no parameters", fill=TRUE)
             } else {
               cat("Symbol map", "with constant values", fill=TRUE)
             }
           },
           discrete = {
             cat("Symbol map", "for discrete inputs:", fill=TRUE)
             print(inputs)
           },
           continuous = {
             cat("Symbol map", "for",
                 switch(rangetype,
                        numeric="real numbers",
                        date = "dates",
                        datetime = "date/time values",
                        unknown = "unrecognised data"),
                 if(!is.null(range)) paste("in", prange(range)) else NULL,
                 fill=TRUE)
           })
    if(!is.null(transform)) {
      parbreak()
      PrintTransformation("Transformation of input values:", transform)
      parbreak()
    }

    ## print graphics parameters
    if(length(parlist) > 0) {
      for(i in seq_along(parlist)) {
        cat(paste0(names(parlist)[i], ": "))
        pari <- parlist[[i]]
        if(!is.function(pari) && length(pari) == 1)
          cat(pari, fill=TRUE) else print(pari)
      }
    }

    if(!is.null(compress)) {
      parbreak()
      PrintTransformation("Scale compression for legend:", compress)
    }
  })
  return(invisible(NULL))
}

PrintTransformation <- function(header, f) {
  cdfclasses <- c("interpolatedCDF", "ewcdf", "ecdf")
  if(samefunction(f, log10)) {
    splat(header, "log10")
  } else if(inherits(f, cdfclasses)) {
    hit <- (inherits(f, cdfclasses, which=TRUE) != 0)
    classname <- cdfclasses[which(hit)][1L]
    splat(header, paren(classname, "["))
  } else {
    splat(header)
    print(f)
  }
  invisible(NULL)
}
  

## Function which actually plots the symbols.
## Called by plot.ppp and plot.symbolmap
## Returns maximum size of symbols

invoke.symbolmap <- local({

  ## plot points, handling various arguments
  do.points <- function(x, y, ...,
                        cex=size, size=NULL, 
                        col=cols, pch=chars, cols=NULL, chars=NULL,
                        lwd=1, etch=FALSE, 
                        do.plot=TRUE) {
    if(do.plot) {
      if(length(cex) == 0) cex <- 1
      if(length(col) == 0) col <- par("col")
      if(length(pch) == 0) pch <- 1
      if(length(lwd) == 0) lwd <- 1
      n <- length(x)
      if(length(cex) == 1) cex <- rep(cex, n)
      if(length(col) == 1) col <- rep(col, n)
      if(length(pch) == 1) pch <- rep(pch, 1)
      if(length(lwd) == 1) lwd <- rep(lwd, n)
      if(length(etch) == 1) etch <- rep(etch, n)
      ## infer which arguments are parallelised
      other <- append(list(...), list(cex=cex, pch=pch))
      isvec <- (lengths(other) == n)
      other.fixed <- other[!isvec]
      other.vec   <- other[isvec]
      ##
      if(any(i <- as.logical(etch))) {
        anti.col <- complementarycolour(col)
        anti.lwd <- if(is.numeric(etch)) etch else 2 * lwd
        do.call.matched(points.default,
                        resolve.defaults(list(x=x[i], y=y[i]),
                                         other.fixed,
                                         lapply(other.vec, "[", i=i),
                                         list(col=anti.col[i],
                                              lwd=anti.lwd[i])),
                        extrargs=c("col", "pch", "type", "bg",
                                   "cex", "lwd", "lty"))
      }
      do.call.matched(points.default,
                    resolve.defaults(list(x=x, y=y),
                                     other,
                                     list(col=col, lwd=lwd)),
                    extrargs=c("col", "pch", "type", "bg", "cex", "lwd", "lty"))
    }
    return(max(cex %orifnull% 1))
  }
  ## plot symbols likewise
  known.shapes <-  c("circles", "squares", "arrows", "crossticks")
  do.symbols <- function(x, y, ..., 
                         shape,
                         size=cex, cex=NULL,
                         fg=col, col=cols, cols=NULL,
                         lwd=1, etch=FALSE,
                         angleref=0,
                         do.plot=TRUE) {
    if(do.plot) {
      if(is.null(size))
        stop("parameter 'size' is required", call.=FALSE)
      ## zap tiny sizes
      tiny <- (size < (max(size)/1000))
      size[tiny] <- 0
      ## collect arguments
      n <- length(x)
      if(length(lwd) == 1) lwd <- rep(lwd, n)
      if(length(etch) == 1) etch <- rep(etch, n)
      if(length(fg) == 0) fg <- rep(par("col"), n) else if(length(fg) == 1) fg <- rep(fg, n)
      if(length(angleref) == 1) angleref <- rep(angleref, n)
      other <- resolve.defaults(list(...),
                                list(add=TRUE, inches=FALSE))
      ## infer which arguments are parallelised
      isvec <- (lengths(other) == n)
      other.fixed <- other[!isvec]
      other.vec   <- other[isvec]
      ##
      if(any(as.logical(etch))) {
        anti.fg <- complementarycolour(fg)
        anti.lwd <- if(is.numeric(etch)) etch else 2 * lwd
      }
      ## plot
      if(any(i <- (shape == "circles") & as.logical(etch))) 
        do.call.matched(symbols,
                        c(list(x=x[i], y=y[i], circles=size[i]/2),
                          other.fixed,
                          lapply(other.vec, "[", i=i),
                          list(lwd=anti.lwd[i], fg=anti.fg[i])),
                        extrargs=c("lwd", "lty"))
      if(any(i <- (shape == "circles")))
        do.call.matched(symbols,
                        c(list(x=x[i], y=y[i], circles=size[i]/2),
                          other.fixed,
                          lapply(other.vec, "[", i=i),
                          list(lwd=lwd[i], fg=fg[i])),
                        extrargs=c("lwd", "lty"))
      if(any(i <- (shape == "squares") & as.logical(etch)))
        do.call.matched(symbols,
                        c(list(x=x[i], y=y[i], squares=size[i]),
                          other.fixed,
                          lapply(other.vec, "[", i=i),
                          list(lwd=anti.lwd[i], fg=anti.fg[i])),
                        extrargs=c("lwd", "lty"))
      if(any(i <- (shape == "squares"))) 
        do.call.matched(symbols,
                        c(list(x=x[i], y=y[i], squares=size[i]),
                          other.fixed,
                          lapply(other.vec, "[", i=i),
                          list(lwd=lwd[i], fg=fg[i])),
                        extrargs=c("lwd", "lty"))
      if(any(i <- (shape == "arrows") & as.logical(etch)))
        do.call.matched(do.arrows,
                        c(list(x=x[i], y=y[i], len=size[i]),
                          other.fixed,
                          lapply(other.vec, "[", i=i),
                          list(lwd=anti.lwd[i], cols=anti.fg[i])),
                        extrargs=c("cols", "col", "lwd", "lty"))
      if(any(i <- (shape == "arrows"))) 
        do.call.matched(do.arrows,
                        c(list(x=x[i], y=y[i], len=size[i]),
                          other.fixed,
                          lapply(other.vec, "[", i=i),
                          list(lwd=lwd[i], cols=fg[i])),
                        extrargs=c("cols", "col", "lwd", "lty"))
      if(any(i <- (shape == "crossticks") & as.logical(etch)))
        do.call.matched(do.crossticks,
                        c(list(x=x[i], y=y[i], len=size[i], angleref=angleref[i]),
                          other.fixed,
                          lapply(other.vec, "[", i=i),
                          list(lwd=anti.lwd[i], cols=anti.fg[i])),
                        extrargs=c("cols", "col", "lwd", "lty"))
      if(any(i <- (shape == "crossticks"))) 
        do.call.matched(do.crossticks,
                        c(list(x=x[i], y=y[i], len=size[i], angleref=angleref[i]),
                          other.fixed,
                          lapply(other.vec, "[", i=i),
                          list(lwd=lwd[i], cols=fg[i])),
                        extrargs=c("cols", "col", "lwd", "lty"))
      if(any(nomatch <- is.na(match(shape, known.shapes)))) {
        unknown <- unique(shape[nomatch])
        nun <- length(unknown)
        warning(paste("Unrecognised",
                      ngettext(nun, "shape", "shapes"),
                      paste0(commasep(sQuote(unknown)), ";"),
                      "recognised values are",
                      commasep(sQuote(known.shapes))),
                call.=FALSE)
      }
    }
    return(max(size))
  }

  do.arrows <- function(x, y, len, direction=0, arrowtype=2, ...,
                        headlength=len * 0.4, 
                        headangle=40,
                        cols=col, col=par('fg'),
                        lwd=1, lty=1) {
    #' vectorise all arguments
    df <- data.frame(x=x, y=y, len=len, direction=direction,
                     arrowtype=arrowtype, headlength=headlength,
                     headangle=headangle, cols=cols, lwd=lwd, lty=lty)
    with(df, {
      alpha <- direction * pi/180
      dx <- len * cos(alpha)/2
      dy <- len * sin(alpha)/2
      x0 <- x - dx
      x1 <- x + dx
      y0 <- y - dy
      y1 <- y + dy
      segments(x0, y0, x1, y1, ..., col=cols, lty=lty, lwd=lwd)
      if(any(arrowtype != 0)) {
        halfangle <- (headangle/2) * pi/180
        beta1 <- alpha + halfangle
        beta2 <- alpha - halfangle
        hx1 <- headlength * cos(beta1)
        hy1 <- headlength * sin(beta1)
        hx2 <- headlength * cos(beta2)
        hy2 <- headlength * sin(beta2)
        if(any(left <- (arrowtype %in% c(1,3)))) {
          segments(x0[left], y0[left], (x0 + hx1)[left], (y0 + hy1)[left],
                   ..., col=cols[left], lwd=lwd[left], lty=lty[left])
          segments(x0[left], y0[left], (x0 + hx2)[left], (y0 + hy2)[left],
                   ..., col=cols[left], lwd=lwd[left], lty=lty[left])
        }
        if(any(right <- (arrowtype %in% c(2,3)))) {
          segments(x1[right], y1[right], (x1 - hx1)[right], (y1 - hy1)[right],
                   ..., col=cols[right], lwd=lwd[right], lty=lty[right])
          segments(x1[right], y1[right], (x1 - hx2)[right], (y1 - hy2)[right],
                   ..., col=cols[right], lwd=lwd[right], lty=lty[right])
        }
      }
    })
    return(invisible(NULL))
  }

  do.crossticks <- function(x, y, len, angleref=0, arrowtype=0, ...,
                            direction = angleref + 90) {
    do.arrows(x=x, y=y, len=len,
              direction=direction,
              arrowtype=arrowtype, ...)
  }

  sanitycheck <- function(df, forbidden, kind) {
    ## 'df' is the result of a symbol map
    ##  (applies to numeric, factor, character)
    bad <- sapply(lapply(df, forbidden), any)
    if(any(bad)) {
      stop(paste("Symbol map produced", kind, "values for",
                 ngettext(sum(bad), "parameter", "parameters"),
                 commasep(sQuote(colnames(df)[bad]))),
           call.=FALSE)
    }
    return(NULL)
  }
    

  ## main function

  invoke.symbolmap <- function(map, values, x=NULL, y=NULL, ...,
                               angleref=NULL, 
                               add=FALSE, do.plot=TRUE,
                               started = add && do.plot) {
    if(!inherits(map, "symbolmap"))
      stop("Argument 'map' should be an object of class 'symbolmap'")
    if(hasxy <- (!is.null(x) || !is.null(y))) {
      xy <- xy.coords(x, y)
      x <- xy$x
      y <- xy$y
    }
    if(is.null(angleref)) 
      angleref <- numeric(length(x)) ## zeroes, or numeric(0)
    ## function will return maximum size of symbols plotted.
    maxsize <- 0
    if(do.plot && !add) plot(x, y, type="n", ...)
    ## map numerical/factor values to graphical parameters
    g <- map(values)
    parnames <- colnames(g)
    ## trap user coding errors etc
    sanitycheck(g, is.na, "NA")
    sanitycheck(g, is.nan, "NaN")
    sanitycheck(g, is.infinite, "infinite")
    ##
    if(do.plot) {
      ## add spatial coordinates
      xydf <- data.frame(x=x, y=y, angleref=angleref)
      if(nrow(xydf) == 0)
        return(invisible(maxsize))
      g <- if(prod(dim(g)) == 0) xydf else 
           do.call(data.frame,
                   c(as.list(g), as.list(xydf), list(stringsAsFactors=FALSE)))
    }
    n <- nrow(g)
    ## figure out which function does the graphics job
    need.points <- any(c("pch", "chars") %in% parnames)
    need.symbols <- "shape" %in% parnames
    if(need.symbols && need.points) {
      worker <- with(g, ifelse(!is.na(shape), "symbols", "points"))
    } else if(need.symbols) {
      worker <- rep.int("symbols", n)
    } else {
      worker <- rep.int("points", n)
    } 
    ## split data according to graphics function involved
    z <- split(g, factor(worker))
    ## display using 'pch'
    zpoints <- z[["points"]]
    if(!is.null(zpoints) && nrow(zpoints) > 0) {
      ms <- do.call(do.points,
                    resolve.defaults(as.list(zpoints),
                                     list(...),
                                     list(do.plot=do.plot)))
      ## value is max(cex)
      ## guess size of one character
      charsize <- if(started || do.plot) max(par('cxy')) else
                  if(hasxy) max(sidelengths(boundingbox(x,y))/40) else 1/40
      maxsize <- max(maxsize, charsize * ms)
    }
    ## display using 'symbols'
    zsymbols <- z[["symbols"]]
    if(!is.null(zsymbols) && nrow(zsymbols) > 0) {
      ms <- do.call(do.symbols,
                    resolve.defaults(as.list(zsymbols),
                                     list(...),
                                     list(do.plot=do.plot)))
      ## ms value is max physical size.
      maxsize <- max(maxsize, ms)
    }
    return(invisible(maxsize))
  }

  invoke.symbolmap
})


## Display the symbol map itself (`legend' style)

plot.symbolmap <- function(x, ..., main,
                           xlim=NULL, ylim=NULL,
                           vertical=FALSE,
                           side=c("bottom", "left", "top", "right"),
                           annotate=TRUE, labelmap=NULL, add=FALSE,
                           nsymbols=NULL, warn=TRUE,
                           colour.only=FALSE, representatives=NULL) {
  if(missing(main))
    main <- short.deparse(substitute(x))
  miss.side <- missing(side)
  if(is.numeric(side)) {
    check.1.integer(side)
    side <- c("bottom", "left", "top", "right")[side]
  }
  side <- match.arg(side) # this may be overwritten later if miss.side=TRUE

  if(colour.only) {
    ## extract only the colour map and plot it
    cmap <- as.colourmap(x)
    if(inherits(cmap, "colourmap")) {
      if(miss.side)
        side <- if(vertical) "right" else "bottom"
      if(!is.numeric(side))
        side <- match(side, c("bottom", "left", "top", "right"))
      result <- plot.colourmap(cmap, ..., main=main,
                               xlim=xlim, ylim=ylim, vertical=vertical,
                               side=side, labelmap=labelmap, add=add)
      return(result)
    }
  }
    
  type <- symbolmaptype(x)
  map <- x
  stuff <- attr(map, "stuff")

  if(type == "constant" && length(stuff$parlist) == 0)
    return(invisible(NULL))

  if(is.null(labelmap)) {
    labelmap <- function(x) x
  } else if(type == "continuous" &&
            is.numeric(labelmap) && length(labelmap) == 1) {
    labscal <- labelmap
    labelmap <- function(x) { x * labscal }
  } else stopifnot(is.function(labelmap))

  ## determine the 'example' input values and their graphical representations
  switch(type,
         constant = {
           vv <- NULL
         },
         continuous = {
           ra <- stuff$range
           vv <- NULL
           if(!is.null(representatives)) {
             vv <- representatives
             if(!all(ok <- inside.range(vv, ra))) {
               nbad <- sum(!ok)
               vv <- vv[ok]
               warning(paste(nbad, "out of", length(vv),
                             ngettext(nbad, "value", "values"),
                             "in the argument",
                             sQuote("representatives"),
                             ngettext(nbad, "was", "were"),
                             "outside the range of the symbol map, and",
                             ngettext(nbad, "was", "were"),
                             "removed"),
                       call.=FALSE)
             }
           } else {
             compress <- stuff$compress
             decompress <- stuff$decompress
             ## check for colour map information 
             if(!is.null(cmap <- as.colourmap(x, warn=FALSE))) {
               st <- attr(cmap, "stuff")
               if(st$discrete) {
                 ## use discrete inputs as representative values
                 vv <- st$inputs
               } else {
                 ## use colour map to determine default range and transformation
                 if(is.null(ra)) ra <- range(st$breaks)
                 if(is.null(compress) && !is.null(st$compress)) {
                   ## nonlinear colour map
                   compress <- st$compress 
                   decompress <- st$decompress
                 }
               }
             }
             if(is.null(vv)) {
               ## representative values not yet determined
               if(is.null(ra))
                 stop("Cannot plot symbolmap with an infinite range")
               if(is.null(compress)) {
                 ## representative values evenly spaced 
                 vv <- prettyinside(ra, n = nsymbols)
               } else {
                 ## representative values evenly spaced on compressed scale
                 vv <- decompress(prettyinside(compress(ra), n = nsymbols))
               }
             }
             if(is.numeric(vv))
               vv <- signif(vv, 4)
           }
         },
         discrete = {
           dd <- stuff$inputs
           if(!is.null(representatives)) {
             vv <- representatives
             if(!all(ok <- vv %in% dd)) {
               nbad <- sum(!ok)
               vv <- vv[ok]
               warning(paste(nbad, "out of", length(vv),
                             ngettext(nbad, "value", "values"),
                             "in the argument",
                             sQuote("representatives"),
                             ngettext(nbad, "was", "were"),
                             "outside the domain of the symbol map, and",
                             ngettext(nbad, "was", "were"),
                             "removed"),
                       call.=FALSE)
             }
           } else {
             vv <- if(is.null(nsymbols)) prettydiscrete(dd) else
                   prettydiscrete(dd, n = nsymbols)
           }
           if(warn && (length(vv) < length(dd))) {
             warning(paste("Only", length(vv), "out of", length(dd),
                           "symbols are shown in the symbol map"),
                     call.=FALSE)
           }
           if(vertical) vv <- rev(vv)
         })
  nn <- length(vv)
  ##    gg <- map(vv)
  ll <- paste(labelmap(vv))
    
  ## determine position of plot and symbols
  if(add) {
    ## x and y limits must respect existing plot space
    usr <- par('usr')
    if(is.null(xlim)) xlim <- usr[1:2]
    if(is.null(ylim)) ylim <- usr[3:4]
  } else {
    ## create new plot
    maxdiam <- invoke.symbolmap(map, vv, do.plot=FALSE, started=FALSE, ...)
    zz <- c(0, max(1, maxdiam))
    if(is.null(xlim) && is.null(ylim)) {
      if(vertical) {
        xlim <- zz
        ylim <- length(vv) * zz
      } else {
        xlim <- length(vv) * zz
        ylim <- zz
      }
    } else if(is.null(ylim)) {
      ylim <- zz
    } else if(is.null(xlim)) {
      xlim <- zz
    }
  }

  ## .......... initialise plot ...............................
  if(!add)
    do.call.matched(plot.default,
                    resolve.defaults(list(x=xlim, y=ylim,
                                          type="n", main=main,
                                          axes=FALSE, xlab="", ylab="",
                                          asp=1.0),
                                     list(...)))
  ## maximum symbol diameter
  maxdiam <- invoke.symbolmap(map, vv, do.plot=FALSE, started=TRUE, ...)

  ## .......... plot symbols ....................
  if(type == "constant") {
    xp <- mean(xlim)
    yp <- mean(ylim)
  } else if(vertical) {
    ## vertical arrangement
    xp <- rep(mean(xlim), nn)
    vskip <- 1.1 * max(maxdiam, 3 * max(strheight(labelmap(vv))))
    if(diff(ylim) > nn * vskip) {
      yp <- (1:nn) * vskip
      yp <- yp - mean(yp) + mean(ylim)
    } else {
      z <- seq(ylim[1], ylim[2], length=nn+1)
      yp <- z[-1] - diff(z)/2
    }
  } else {
    ## horizontal arrangement
    yp <- rep(mean(ylim), nn)
    hskip <- 1.1 * max(maxdiam, max(strwidth(labelmap(vv))))
    if(diff(xlim) > nn * hskip) {
      xp <- (1:nn) * hskip
      xp <- xp - mean(xp) + mean(xlim)
    } else {
      z <- seq(xlim[1], xlim[2], length=nn+1)
      xp <- z[-1] - diff(z)/2
    }
  }
  invoke.symbolmap(map, vv, xp, yp, ..., add=TRUE,
                   angleref=if(vertical) 90 else 0)

  ## ................. draw annotation ..................
  dotargs <- list(...)
  axiscol <- dotargs$col
  nac <- length(axiscol)
  if(nac > 0 && (!is.colour(axiscol) || nac > 1)) {
    ## only a single colour is permitted for 'axis'
    dotargs$col <- NULL
  }
  if(annotate && length(ll) > 0) {
    if(vertical) {
      ## default axis position is to the right 
      if(miss.side) side <- "right"
      sidecode <- match(side, c("bottom", "left", "top", "right"))
      if(!(sidecode %in% c(2,4)))
        warning(paste("side =", sQuote(side),
                      "is not consistent with vertical orientation"))
      pos <- c(ylim[1], xlim[1], ylim[2], xlim[2])[sidecode]
      ## draw axis
      do.call.matched(graphics::axis,
                      resolve.defaults(dotargs,
                                       list(side=sidecode, pos=pos, at=yp,
                                            labels=ll, tick=FALSE, las=1)),
                      extrargs=graphicsPars("axis"))
    } else {
      ## default axis position is below 
      if(miss.side) side <- "bottom"
      sidecode <- match(side, c("bottom", "left", "top", "right"))
      if(!(sidecode %in% c(1,3)))
        warning(paste("side =", sQuote(side),
                      "is not consistent with horizontal orientation"))
      pos <- c(ylim[1], xlim[1], ylim[2], xlim[2])[sidecode]
      ## draw axis
      do.call.matched(graphics::axis,
                      resolve.defaults(dotargs,
                                       list(side = sidecode, pos = pos,
                                            at = xp, labels=ll, tick=FALSE)),
                      extrargs=graphicsPars("axis"))
    } 
  }
  return(invisible(NULL))
}

plan.legend.layout <- function(B, 
                               ..., 
                               side=c("bottom", "left", "top", "right"),
                               sep=NULL,
                               leg.size=NULL,
                               sep.frac=0.05,
                               size.frac=0.05,
                               started=FALSE,
                               map=NULL) {
  ## Determine size and position of a box containing legend or symbolmap
  ## attached to a plot in region 'B'.
  ##   sep, leg.size are absolute distances;
  ##   sep.frac, size.frac are fractions of the maximum sidelength of B.
  side <- match.arg(side)
  B <- as.rectangle(B)
  Bsize <- max(sidelengths(B))
  if(is.null(leg.size)) {
    leg.size <- size.frac * Bsize
  } else {
    check.1.real(leg.size)
    stopifnot(leg.size > 0)
  }
  if(is.null(sep)) {
    sep <- sep.frac * Bsize
  } else {
    check.1.real(sep)
    stopifnot(sep > 0)
  }
  if(is.null(map) || !inherits(map, "symbolmap")) {
    vv <- NULL
    textlength <- 8
  } else {
    vv <- with(attr(map, "stuff"),
               if(type == "discrete") inputs else prettyinside(range))
    textlength <- max(nchar(paste(vv)))
  }
  if(started && !is.null(vv)) {
    textwidth <- max(strwidth(vv))
    textheight <- max(strheight(vv))
  } else {
    ## the plot has not been initialised: guess character size
    charsize <- diff(if(side %in% c("left", "right")) B$yrange else B$xrange)/40
    textwidth <- charsize * textlength
    textheight <- charsize
  }
  switch(side,
         right={
           ## symbols to right of image
           b <- owinInternalRect(B$xrange[2] + sep + c(0, leg.size),
                     B$yrange)
           ## text to right of symbols
           tt <- owinInternalRect(b$xrange[2] + sep + c(0, textwidth),
                      b$yrange)
           iside <- 4
         },
         left={
           ## symbols to left of image
           b <- owinInternalRect(B$xrange[1] - sep - c(leg.size, 0),
                     B$yrange)
           ## text to left of symbols
           tt <- owinInternalRect(b$xrange[1] - sep - c(textwidth, 0),
                      b$yrange)
           iside <- 2
         },
         top={
           ## symbols above image
           b <- owinInternalRect(B$xrange,
                     B$yrange[2] + sep + c(0, leg.size))
           ## text above symbols
           tt <- owinInternalRect(b$xrange,
                      b$yrange[2] + 3* charsize + c(0, textheight))
           iside <- 3
         },
         bottom={
           ## symbols below image
           b <- owinInternalRect(B$xrange,
                     B$yrange[1] - sep - c(leg.size, 0))
           ## text below symbols
           tt <- owinInternalRect(b$xrange,
                      b$yrange[1] - 3 * charsize - c(textheight, 0))
           iside <- 1
         })
  A <- boundingbox(B, b, tt)
  return(list(A=A, B=B, b=b, tt=tt,
              iside=iside, side=side, size=leg.size, charsize=charsize, sep=sep))
}


as.colourmap.symbolmap <- function(x, ..., warn=TRUE) {
  ## extract only the colour map and plot it
  parlist <- attr(x, "stuff")$parlist
  iscol <- sapply(parlist, inherits, what="colourmap")
  nc <- sum(iscol)
  if(nc == 0) {
    if(warn) warning("No colour map information was detected", call.=FALSE)
    return(NULL)
  }
  used <- which(iscol)[[1L]]
  cmap <- parlist[[used]]
  if(nc > 1 && warn && length(unique(parlist[iscol])) > 1)
    warning(paste("More than one colour map was detected;",
                  "using the colour map for",
                  sQuote(names(parlist)[used])),
            call.=FALSE)
    return(cmap)
}

summary.symbolmap <- function(object, ...) {
  st <- attr(object, "stuff")
  typ <- st[["type"]]
  dom <- switch(typ,
                constant = { integer(0) },
                discrete = { st$inputs },
                continuous = { st$range })
  parlist <- st[["parlist"]]
  parnames <- names(parlist)
  iscolmap <- sapply(parlist, inherits, what="colourmap")
  isatom <- sapply(parlist, is.atomic)
  lenfs <- lengths(parlist)
  isconstant <- isatom & (lenfs == 1)
  if(any(iscolmap)) 
    isconstant[iscolmap] <- (lengths(lapply(parlist[iscolmap], colouroutputs)) == 1)
  colournames <- c("col", "cols", "fg", "bg", "border", "fill")
  shapenames <- c("shape", "pch", "chars", "direction", "arrowtype", "headlength", "headangle", "etch")
  sizenames <- c("size", "cex", "headlength")
  physicalsizenames <- c("size", "headlength")
  iscolour <- iscolmap | (parnames %in% colournames)
  isshape <- parnames %in% shapenames
  issize <- parnames %in% sizenames
  isphysical <- parnames %in% physicalsizenames
  fixedcolour <- all(isconstant[iscolour])
  fixedshape <- all(isconstant[isshape])
  fixedsize <- all(isconstant[issize])
  z <- list(type        = typ,
            domain      = dom, 
            pars        = parnames,
            colmaps     = parnames[iscolmap],
            rangetype   = if(typ == "continuous") st[["rangetype"]] else NULL,
            range       = if(typ == "continuous") st[["range"]] else NULL,
            isconstant  = isconstant, # vector
            iscolour    = iscolour, # vector
            issize      = issize, # vector
            isshape     = isshape, # vector
            isphysical  = isphysical, # vector
            fixedshape  = fixedshape, # logical(1)
            fixedsize   = fixedsize,  # logical(1)
            fixedcolour = fixedcolour # logical(1)
            )
  class(z) <- c("summary.symbolmap", class(z))
  return(z)
}

print.summary.symbolmap <- function(x, ...) {
  with(x, {
    switch(type,
           constant = {
             if(length(pars) == 0) {
               cat("Symbol map", "with no parameters", fill=TRUE)
             } else {
               cat("Symbol map", "with constant values", fill=TRUE)
             }
           },
           discrete = {
             cat("Symbol map", "for discrete inputs:", fill=TRUE)
             print(domain)
           },
           continuous = {
             cat("Symbol map", "for",
                 switch(rangetype,
                        numeric="real numbers",
                        date = "dates",
                        datetime = "date/time values",
                        unknown = "unrecognised data"),
                 if(!is.null(range)) paste("in", prange(range)) else NULL,
                 fill=TRUE)
           })
    if(length(pars) > 0) {
      splat("Graphics parameters defined:")
      splat(paste("\t",
                  sQuote(pars),
                  "\t",
                  paren(ifelse(isconstant, "constant", "variable")),
                  "\t",
                  ifelse(pars %in% colmaps, "(colour map)", ""),
                  "\n"))
      if(all(isconstant)) {
        splat("All graphics parameters are constant")
      } else {
        att <- c("size", "shape", "colour")
        fux <- c(fixedsize, fixedshape, fixedcolour)
        attfux <- att[fux]
        attnon <- att[!fux]
        blurb <- "Symbols have"
        if(length(attfux))
          blurb <- paste(blurb, "fixed", commasep(attfux, " and "))
        if(length(attnon)) {
          if(length(attfux)) blurb <- paste0(blurb, ", but")
          blurb <- paste(blurb, "variable", commasep(attnon, " and "))
        }
        splat(blurb)
      }
    }
    return(invisible(NULL))
  })
}

default.symbolmap <- function(x, ...) {
  UseMethod("default.symbolmap")
}
