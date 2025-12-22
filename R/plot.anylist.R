##
##  plot.anylist.R
##
##  Plotting functions for 'solist', 'anylist', 'imlist'
##       and legacy class 'listof'
##
##  $Revision: 1.49 $ $Date: 2025/12/22 08:00:18 $
##

plot.anylist <- plot.solist <- plot.listof <-
  local({

  ## auxiliary functions
  classes.with.do.plot <- c("colourmap",
                            "diagramobj", "diagppm",
                            "fv", 
                            "im", "layered",
                            "linnet", "lintess", "linim", "linfun", "lpp",
                            "linearquadratcount", 
                            "msr",
                            "owin",
                            "ppp", "psp",
                            "ssf", "symbolmap",
                            "tess", "texturemap")
  
  classes.with.multiplot <- c("ppp", "lpp", "msr", "tess",
                              "leverage.ppm", "influence.ppm")

  has.multiplot <- function(x) {
    inherits(x, classes.with.multiplot) ||
      (is.function(x) && "multiplot" %in% names(formals(x)))
  }
  
  has.do.plot <- function(x) {
    inherits(x, classes.with.do.plot) ||
      (is.function(x) && "do.plot" %in% names(formals(x)))
  }
  
  extraplot <- function(nnn, x, ..., add=FALSE, do.plot=TRUE, extrargs=list(),
                        panel.args=NULL, plotcommand="plot", pctype) {
    argh <- list(...)
    if(missing(pctype))
      pctype <- recognise.plotcommand.type(plotcommand)
    if(has.multiplot(x) && pctype == "plot")
      argh <- c(argh, list(multiplot=FALSE))
    if(!do.plot) {
      if(has.do.plot(x) && (pctype %in% c("plot", "image", "contour"))) {
        argh <- c(argh, list(do.plot=FALSE))
      } else {
        return(NULL)
      }
    }
    if(!is.null(panel.args)) {
      xtra <- if(is.function(panel.args)) panel.args(nnn) else panel.args
      if(!is.list(xtra))
        stop(paste0("panel.args",
                    if(is.function(panel.args)) "(i)" else "",
                    " should be a list"))
      argh <- resolve.defaults(xtra, argh)
    }
    if(length(extrargs) > 0)
      argh <- resolve.defaults(argh, extrargs)
    ## some plot commands don't recognise 'add'
    if(add)
      argh <- append(argh, list(add=TRUE))
    do.call(plotcommand, append(list(x=x), argh))
  }

  exec.or.plot <- function(cmd, i, xi, ..., extrargs=list(), add=FALSE) {
    if(is.null(cmd)) return(NULL)
    argh <-
      resolve.defaults(list(...),
                       extrargs,
                       ## some plot commands don't recognise 'add' 
                       if(add) list(add=TRUE) else NULL,
                       if(has.multiplot(cmd)) list(multiplot=FALSE) else NULL)
    if(is.function(cmd)) {
      force(xi)
      do.call(cmd, resolve.defaults(list(i, quote(xi)), argh))
    } else {
      do.call(plot, resolve.defaults(list(cmd), argh))
    }
  }

  exec.or.plotshift <- function(cmd, i, xi, ..., vec=vec,
                                extrargs=list(), add=FALSE) {
    if(is.null(cmd)) return(NULL)
    argh <-
      resolve.defaults(list(...),
                       extrargs,
                       ## some plot commands don't recognise 'add' 
                       if(add) list(add=TRUE) else NULL,
                       if(has.multiplot(cmd)) list(multiplot=FALSE) else NULL)
    if(is.function(cmd)) {
      force(xi)
      do.call(cmd, resolve.defaults(list(i, quote(xi)), argh))
    } else {
      cmd <- shift(cmd, vec)
      do.call(plot, resolve.defaults(list(quote(cmd)), argh))
    }
  }

  ## bounding box, including ribbon for images, legend for point patterns
  getplotbox <- function(x, ..., do.plot,
                         plotcommand="plot", pctype, multiplot) {
    if(missing(pctype))
      pctype <- recognise.plotcommand.type(plotcommand)
    if(inherits(x, classes.with.do.plot)) {
      switch(pctype,
             plot = {
               y <- if(has.multiplot(x)) {
                      plot(x, ..., multiplot=FALSE, do.plot=FALSE)
                    } else {
                      plot(x, ..., do.plot=FALSE)
                    }
               return(as.owin(y))
             },
             contour = {
               y <- contour(x, ..., do.plot=FALSE)      
               return(as.owin(y))
             },
             image = {
               y <- image(x, ..., do.plot=FALSE)
               return(as.owin(y))
             },
             unknown = {
               plc <- plotcommand
               if(is.character(plc)) plc <- get(plc)
               if(!is.function(plc)) stop("Unrecognised plot function")
               if("do.plot" %in% names(args(plc))) {
                 if(has.multiplot(plc)) {
                   y <- do.call(plc, list(x=x, ..., do.plot=FALSE,
                                          multiplot=FALSE))
                 } else {
                   y <- do.call(plc, list(x=x, ..., do.plot=FALSE))
                 }
               }
               return(as.owin(y))
             })
    } 
    return(try(as.rectangle(x), silent=TRUE))
  }

  # calculate bounding boxes for each panel using intended arguments!
  getPlotBoxes <- function(xlist, ..., panel.args=NULL, extrargs=list()) {
    userargs <- list(...)
    n <- length(xlist)
    result <- vector(length=n, mode="list")
    for(i in seq_len(n)) {
      pai <- if(is.function(panel.args)) panel.args(i) else list()
      argh <- resolve.defaults(pai, userargs, extrargs)
      xxi <- xlist[[i]]
      result[[i]] <- do.call(getplotbox, append(list(x=quote(xxi)), argh))
    }
    return(result)
  }
    
  is.shiftable <- function(x) {
    if(is.null(x)) return(TRUE)
    if(is.function(x)) return(FALSE)
    y <- try(as.rectangle(x), silent=TRUE)
    return(!inherits(y, "try-error"))
  }

  maxassigned <- function(i, values) max(-1, values[i[i > 0]])

  plotadornment <- function(adorn, adorn.args, ...) {
    aname <- deparse(substitute(adorn))
    if(is.null(adorn)) {
      z <- NULL
    } else if(is.function(adorn)) {
      z <- do.call(adorn, resolve.defaults(adorn.args, list(...)))
    } else if(inherits(adorn, c("colourmap", "symbolmap"))) {
      z <- do.call(plot, resolve.defaults(list(x=adorn),
                                          adorn.args,
                                          list(...)))
    } else warning("Unrecognised format for", sQuote(aname))
    return(z)
  }
  
  plot.anylist <- function(x, ..., main, arrange=TRUE,
                            nrows=NULL, ncols=NULL,
                            main.panel=NULL,
                            mar.panel=c(2,1,1,2),
                            hsep = 0,
                            vsep = 0,
                            panel.begin=NULL,
                            panel.end=NULL,
                            panel.args=NULL,
                            panel.begin.args=NULL,
                            panel.end.args=NULL,
                            panel.vpad = 0.2,
                            plotcommand="plot",
                            do.plot=TRUE,
                            adorn.left=NULL,
                            adorn.right=NULL,
                            adorn.top=NULL,
                            adorn.bottom=NULL,
                            adorn.size=0.2,
                            adorn.args=list(),
                            equal.scales=FALSE,
                            halign=FALSE, valign=FALSE
                           ) {
    xname <- short.deparse(substitute(x))

    ## recursively expand entries which are 'anylist' etc
    while(any(sapply(x, inherits, what="anylist"))) 
      x <- as.solist(expandSpecialLists(x, "anylist"), demote=TRUE)
    
    isSo <- inherits(x, "solist")
    isIm <- inherits(x, "imlist") || (isSo && all(unlist(lapply(x, is.im))))
    
    cl <- match.call()
    if(isIm && missing(plotcommand)) {
      ## dispatch to image.imlist instead
      cl[[1]] <- as.name("image.imlist")
      parenv <- sys.parent()
      return(invisible(eval(cl, envir=parenv)))
    }

    ## recognise type of plot
    pctype <- recognise.plotcommand.type(plotcommand)
    
    ## determine whether 'fv' objects are present
    if(isSo) {
      allfv <- somefv <- FALSE
    } else {
      isfv <- unlist(lapply(x, is.fv))
      allfv <- all(isfv)
      somefv <- any(isfv)
      if(somefv && !requireNamespace("spatstat.explore"))
        stop(paste("Package 'spatstat.explore' is required",
                   "for plotting objects of class 'fv'"),
             call.=FALSE)
    }
    
    ## panel margins
    if(!missing(mar.panel)) {
      nm <- length(mar.panel)
      if(nm == 1) mar.panel <- rep(mar.panel, 4) else
      if(nm == 2) mar.panel <- rep(mar.panel, 2) else
      if(nm != 4) stop("mar.panel should have length 1, 2 or 4")
    } else if(somefv) {
      ## change default
      mar.panel <- 0.25+c(4,4,2,2)
    }
    
    n <- length(x)
    names(x) <- good.names(names(x), "Component_", 1:n)
    if(is.null(main.panel))
      main.panel <- names(x)
    else {
      if(!is.expression(main.panel))
        main.panel <- as.character(main.panel)
      nmp <- length(main.panel)
      if(nmp == 1)
        main.panel <- rep.int(main.panel, n)
      else if(nmp != n)
        stop("Incorrect length for main.panel")
    }

    if(allfv && equal.scales) {
      ## all entries are 'fv' objects: determine their plot limits
      fvlims <- lapply(x, plot, ..., limitsonly=TRUE)
      ## establish common x,y limits for all panels
      xlim <- range(unlist(lapply(fvlims, getElement, name="xlim")))
      ylim <- range(unlist(lapply(fvlims, getElement, name="ylim")))
      extrargs <- list(xlim=xlim, ylim=ylim)
    } else extrargs <- list()

    extrargs.begin <- resolve.defaults(panel.begin.args, extrargs)
    extrargs.end <- resolve.defaults(panel.end.args, extrargs)

    ## adornments
    adornments <- list(adorn.left   = adorn.left,
                       adorn.right  = adorn.right,
                       adorn.top    = adorn.top,
                       adorn.bottom = adorn.bottom)
    adornments <- adornments[!sapply(adornments, is.null)]
    nadorn <- length(adornments)
    adorable <- all(sapply(adornments, inherits,
                           what=c("symbolmap", "colourmap")))
    ## "texturemap" not yet supported by plan.legend.layout
    
    if(!arrange) {
      ## ------------  sequence of plots -----------------------
      if(!do.plot) {
        result <- vector(mode="list", length=n)
        if(allfv) return(result) # physical size is undefined
        switch(pctype,
               persp = {
                 return(result) # physical size undefined
               },
               image = ,
               contour = ,
               plot = {
                 ## compute colour map (etc) and layout boxes
                 for(i in 1:n) {
                   xi <- x[[i]]
                   result[[i]] <-
                     extraplot(i, xi, ...,
                               do.plot=FALSE,
                               add=!is.null(panel.begin),
                               main=main.panel[i],
                               panel.args=panel.args, extrargs=extrargs,
                               plotcommand=plotcommand,
                               pctype=pctype) %orifnull% list()
                 }
               },
               unknown = {
                 ## attempt to determine physical size of each plot
                 result <- getPlotBoxes(x, ...,
                                        plotcommand=plotcommand,
                                        pctype=pctype,
                                        panel.args=panel.args,
                                        extrargs=extrargs)
                 if(any(bad <- sapply(result, inherits, what="try-error")))
                   result[bad] <- list(NULL)
               })
        ## no plotting performed - exit
        return(result)
      }
      ## ----------- start actual plotting ------------------
      result <- vector(mode="list", length=n)
      ## generate each plot
      for(i in 1:n) {
        xi <- x[[i]]
        exec.or.plot(panel.begin, i, xi, main=main.panel[i],
                     extrargs=extrargs.begin)
        result[[i]] <-
          extraplot(i, xi, ...,
                    add=!is.null(panel.begin),
                    main=main.panel[i],
                    panel.args=panel.args, extrargs=extrargs,
                    plotcommand=plotcommand, pctype=pctype) %orifnull% list()
        exec.or.plot(panel.end, i, xi, add=TRUE, extrargs=extrargs.end)
      }
      ## 
      if(nadorn > 0)
        warning(paste(ngettext(nadorn, "Argument", "Arguments"),
                      commasep(sQuote(names(adornments))),
                      ngettext(nadorn, "was", "were"),
                      "ignored because arrange=FALSE"),
                call.=FALSE)

      return(invisible(result))
    }

    ## ----------------------------------------------------
    ## >>>>>>>>>>>>>>>> ARRAY of plots <<<<<<<<<<<<<<<<<<<<
    ## ----------------------------------------------------

    ## decide whether to plot a main header
    main <- if(!missing(main) && !is.null(main)) main else xname
    if(!is.character(main)) {
      ## main title could be an expression
      nlines <- 1
      banner <- TRUE
    } else {
      ## main title is character string/vector, possibly ""
      banner <- any(nzchar(main))
      if(length(main) > 1)
        main <- paste(main, collapse="\n")
      nlines <- length(unlist(strsplit(main, "\n")))
    }
    ## determine arrangement of plots
    ## arrange like mfrow(nrows, ncols) plus a banner at the top
    if(is.null(nrows) && is.null(ncols)) {
      nrows <- as.integer(floor(sqrt(n)))
      ncols <- as.integer(ceiling(n/nrows))
    } else if(!is.null(nrows) && is.null(ncols))
      ncols <- as.integer(ceiling(n/nrows))
    else if(is.null(nrows) && !is.null(ncols))
      nrows <- as.integer(ceiling(n/ncols))
    else stopifnot(nrows * ncols >= length(x))
    nblank <- ncols * nrows - n
    if(allfv || pctype == "persp") {
      ## Function plots do not have physical 'size'
      sizes.known <- FALSE
    } else {
      ## Determine dimensions of objects
      ##     (including space for colour ribbons, if they are images)
      boxes <- getPlotBoxes(x, ..., plotcommand=plotcommand, pctype=pctype,
                            panel.args=panel.args, extrargs=extrargs)
      sizes.known <- !any(sapply(boxes, inherits, what="try-error"))
      sizes.known <- sizes.known && (nadorn == 0 || adorable)
      if(sizes.known) {
        extrargs <- resolve.defaults(extrargs, list(claim.title.space=TRUE))
        boxes <- getPlotBoxes(x, ..., plotcommand=plotcommand, pctype=pctype,
                              panel.args=panel.args, extrargs=extrargs)
      }
      if(equal.scales && !sizes.known) {
        warning("Ignored equal.scales=TRUE; scales could not be determined")
        equal.scales <- FALSE
      }
    }
    if(sizes.known) {
      ## determine size of each panel
      if(equal.scales) {
        ## do not rescale panels
        scaledboxes <- boxes
      } else {
        ## rescale panels
        sides <- lapply(boxes, sidelengths)
        bwidths <- unlist(lapply(sides, "[", 1))
        bheights <- unlist(lapply(sides, "[", 2))
        ## Force equal heights, unless there is only one column
        scales <- if(ncols > 1) 1/bheights else 1/bwidths
        if(all(is.finite(scales))) {
          scaledboxes <- vector(mode="list", length=n)
          for(i in 1:n)
            scaledboxes[[i]] <- scalardilate(boxes[[i]], scales[i])
        } else {
          #' uh-oh
          equal.scales <- sizes.known <- FALSE
          scaledboxes <- boxes
        }
      }
    }
    ## determine whether to display all objects in one enormous plot
    ## Precondition is that everything has a spatial bounding box
    single.plot <- equal.scales && sizes.known
    if(equal.scales && !single.plot && !allfv)
      warning("equal.scales=TRUE ignored ", "because bounding boxes ",
              "could not be determined", call.=FALSE)
    ## enforce alignment by expanding boxes
    if(halign) {
      if(!equal.scales)
        warning("halign=TRUE ignored because equal.scales=FALSE")
      ## x coordinates align in each column
      xr <- range(sapply(scaledboxes, getElement, name="xrange"))
      scaledboxes <- lapply(scaledboxes, "[[<-", i="xrange", value=xr)
    }
    if(valign) {
      if(!equal.scales)
        warning("valign=TRUE ignored because equal.scales=FALSE")
      ## y coordinates align in each column
      yr <- range(sapply(scaledboxes, getElement, name="yrange"))
      scaledboxes <- lapply(scaledboxes, "[[<-", i="yrange", value=yr)
    }
    ## set up layout
    mat <- matrix(c(seq_len(n), integer(nblank)),
                  byrow=TRUE, ncol=ncols, nrow=nrows)
    if(sizes.known) {
      boxsides <- lapply(scaledboxes, sidelengths)
      xwidths <- sapply(boxsides, "[", i=1)
      xheights <- sapply(boxsides, "[", i=2)
      heights <- apply(mat, 1, maxassigned, values=xheights)
      widths <- apply(mat, 2, maxassigned, values=xwidths)
    } else {
      heights <- rep.int(1, nrows)
      widths <- rep.int(1, ncols)
    }
    #' negative heights/widths arise if a row/column is not used.
    meanheight <- mean(heights[heights > 0])
    meanwidth  <- mean(widths[heights > 0])
    heights[heights <= 0] <- meanheight
    widths[widths <= 0] <- meanwidth
    nall <- n
    ##
    if(single.plot) {
      ## .........  create a single plot ..................
      ## determine sizes
      ht <- max(heights)
      wd <- max(widths)
      marpar <- mar.panel * c(ht, wd, ht, wd)/6
      vsep <- vsep * ht/6
      hsep <- hsep * wd/6
      mainheight <- any(nzchar(main.panel)) * ht/5
      ewidths <- marpar[2] + widths + marpar[4]
      eheights <- marpar[1] + heights + marpar[3] + mainheight
      ## create box delimiting all panels
      Width <- sum(ewidths) + hsep * (length(ewidths) - 1)
      Height <- sum(eheights) + vsep * (length(eheights) - 1)
      bigbox <- owinInternalRect(c(0, Width), c(0, Height))
      ## bottom left corner of each panel 
      ox <- marpar[2] + cumsum(c(0, ewidths + hsep))[1:ncols]
      oy <- marpar[1] + cumsum(c(0, rev(eheights) + vsep))[nrows:1]
      panelorigin <- as.matrix(expand.grid(x=ox, y=oy))
      ## add space for adornments (colour maps or symbol maps)
      if(nadorn > 0) {
        ## box containing spatial objects but excluding annotations
        rx <- c(min(ox), max(ox) + widths[length(widths)])
        ry <- c(min(oy), max(oy) + heights[length(heights)])
        actionbox <-  owinInternalRect(rx, ry)
        ## calculate extensions 
        sidestrings <- sub("adorn.", "", names(adornments))
        sideplans <- mapply(plan.legend.layout,
                            MoreArgs=list(B=actionbox),
                            side=sidestrings,
                            map=unname(adornments),
                            SIMPLIFY=FALSE)
        ## extract box for each adornment
        sideboxes <- lapply(sideplans, getElement, name="b")
        ## update bigbox to contain them all
        coveringboxes <- lapply(sideplans, getElement, name="A")
        bigbox <- do.call(boundingbox, append(list(bigbox), coveringboxes))
      }
      ## initialise, with banner
      if(do.plot) {
        cex <- resolve.1.default(list(cex.title=1.5), list(...))/par('cex.main')
        plot(bigbox, type="n", main=main, cex.main=cex)
      }
      ## plot individual objects
      result <- vector(mode="list", length=n)
      for(i in 1:n) {
        ## determine shift vector that moves bottom left corner of spatial box
        ## to bottom left corner of target area on plot device
        vec <- panelorigin[i,] - with(scaledboxes[[i]], c(xrange[1], yrange[1]))
        ## shift panel contents
        xi <- x[[i]]
        xishift <- shift(xi, vec)
        ## let rip
        if(do.plot && !is.null(panel.begin))
          exec.or.plotshift(panel.begin, i, xishift,
                            add=TRUE,
                            main=main.panel[i], show.all=TRUE,
                            extrargs=extrargs.begin,
                            vec=vec)
        result[[i]] <-
          extraplot(i, xishift, ..., do.plot=do.plot,
                    add=TRUE, show.all=is.null(panel.begin),
                    main=main.panel[i],
                    extrargs=extrargs,
                    panel.args=panel.args,
                    plotcommand=plotcommand, pctype=pctype) %orifnull% list()
        if(do.plot)
          exec.or.plotshift(panel.end, i, xishift, add=TRUE,
                            extrargs=extrargs.end,
                            vec=vec)
      }
      ## add adornments if any
      if(do.plot) {
        for(i in seq_len(nadorn)) {
          bi <- sideboxes[[i]]
          do.call(plot,
                  resolve.defaults(list(x=adornments[[i]]),
                                   list(...),
                                   list(add=TRUE,
                                        xlim=bi$xrange, ylim=bi$yrange,
                                        side=sidestrings[i]),
                                   adorn.args))
        }
      }
      return(invisible(result))
    }
    ## ................. multiple logical plots using 'layout' ..............
    ## adjust panel margins to accommodate desired extra separation
    mar.panel <- pmax(0, mar.panel + c(vsep, hsep, vsep, hsep)/2)
    ## increase heights to accommodate panel titles
    if(sizes.known && any(nzchar(main.panel))) 
      heights <- heights * (1 + panel.vpad)
    ## check for adornment
    if(!is.null(adorn.left)) {
      ## add margin at left, of width adorn.size * meanwidth
      nall <- i.left <- n+1
      mat <- cbind(i.left, mat)
      widths <- c(adorn.size * meanwidth, widths)
    } 
    if(!is.null(adorn.right)) {
      ## add margin at right, of width adorn.size * meanwidth
      nall <- i.right <- nall+1
      mat <- cbind(mat, i.right)
      widths <- c(widths, adorn.size * meanwidth)
    } 
    if(!is.null(adorn.bottom)) {
      ## add margin at bottom, of height adorn.size * meanheight
      nall <- i.bottom <- nall+1
      mat <- rbind(mat, i.bottom)
      heights <- c(heights, adorn.size * meanheight)
    } 
    if(!is.null(adorn.top)) {
      ## add margin at top, of height adorn.size * meanheight
      nall <- i.top <- nall + 1
      mat <- rbind(i.top, mat)
      heights <- c(adorn.size * meanheight, heights)
    } 
    if(banner) {
      ## Increment existing panel numbers
      ## New panel 1 is the banner
      panels <- (mat > 0)
      mat[panels] <- mat[panels] + 1
      mat <- rbind(1, mat)
      heights <- c(0.1 * meanheight * (1 + nlines), heights)
    }
    ## -------- start actual plotting -----------------
    if(do.plot) {
      ## declare layout
      layout(mat, heights=heights, widths=widths, respect=sizes.known)
      ## .... plot banner
      if(banner) {
        opa <- par(mar=rep.int(0,4), xpd=TRUE)
        on.exit(par(opa))
        plot(numeric(0),numeric(0),type="n",ann=FALSE,axes=FALSE,
             xlim=c(-1,1),ylim=c(-1,1))
        cex <- resolve.1.default(list(cex.title=1.5), list(...))/par('cex')
        text(0,0,main, cex=cex)
      }
      ## plot panels
      npa <- par(mar=mar.panel)
      if(!banner) on.exit(par(npa))
    }
    result <- vector(mode="list", length=n)
    for(i in 1:n) {
      xi <- x[[i]]
      if(do.plot)
        exec.or.plot(panel.begin, i, xi, main=main.panel[i],
                     extrargs=extrargs.begin)
      result[[i]] <-
        extraplot(i, xi, ..., do.plot=do.plot,
                  add = !is.null(panel.begin), 
                  main = main.panel[i],
                  extrargs=extrargs,
                  panel.args=panel.args,
                  plotcommand=plotcommand, pctype=pctype) %orifnull% list()
      if(do.plot)
        exec.or.plot(panel.end, i, xi, add=TRUE, extrargs=extrargs.end)
    }
    ## adornments
    if(do.plot && nall > n) {
      par(mar=rep.int(0,4), xpd=TRUE)
      plotadornment(adorn.left, adorn.args)
      plotadornment(adorn.right, adorn.args)
      plotadornment(adorn.bottom, adorn.args)
      plotadornment(adorn.top, adorn.args)
    }
    ## revert
    if(do.plot)
      layout(1)
    return(invisible(result))
  }

  
  plot.anylist
})


contour.imlist <- contour.listof <- function(x, ...) {
  xname <- short.deparse(substitute(x))
  force(x)
  do.call(plot.solist,
          resolve.defaults(list(x=quote(x), plotcommand="contour"),
                           list(...),
                           list(main=xname)))
}

plot.imlist <- local({

  plot.imlist <- function(x, ..., plotcommand="image",
                          equal.ribbon = FALSE,
                          equal.scales = FALSE,
                          ribmar=NULL) {
    xname <- short.deparse(substitute(x))
    force(x)
    if(missing(plotcommand) &&
       any(sapply(x, inherits, what=c("linim", "linfun"))))
      plotcommand <- "plot"
    pctype <- recognise.plotcommand.type(plotcommand)
    if(equal.ribbon && pctype %in% c("image", "plot")) {
      out <- imagecommon(x, ...,
                         xname=xname,
                         ribmar=ribmar,
                         equal.scales=equal.scales)
    } else {
      out <- do.call(plot.solist,
                     resolve.defaults(list(x=quote(x),
                                           plotcommand=plotcommand), 
                                      list(...),
                                      list(main=xname,
                                           equal.scales=equal.scales)))
    }
    return(invisible(out))
  }

  imagecommon <- function(x, ...,
                          xname,
                          zlim=NULL, log=FALSE,
                          equal.scales=FALSE,
                          do.plot=TRUE, 
                          ribbon=TRUE,
                          ribside=c("right", "left", "bottom", "top"),
                          ribsep=NULL, ribwid=0.5, ribn=1024,
                          ribscale=NULL, ribargs=list(),
                          ribmar = NULL, mar.panel = c(2,1,1,2)) {
    if(missing(xname))
      xname <- short.deparse(substitute(x))
    force(x)
    ribside <- match.arg(ribside)
    stopifnot(is.list(ribargs))
    if(!is.null(ribsep))
      warning("Argument ribsep is not yet implemented for image arrays")
    ## ascertain types of pixel values
    xtypes <- sapply(x, getElement, name="type")
    ischar <- (xtypes == "character")
    if(any(ischar)) {
      ## convert character-valued images to factor-valued
      strings <- unique(unlist(lapply(x[ischar], "[")))
      x[ischar] <- lapply(x[ischar], factorimage, levels=strings)
      xtypes[ischar] <- "factor"
    }
    isfactor <- xtypes == "factor"
    isnumeric <- xtypes %in% c("real", "integer", "logical")
    if(all(isnumeric)) {
      ## determine range of values for colour map
      if(is.null(zlim))
        zlim <- range(unlist(lapply(x, range)))
      ## determine common colour map based on zlim
      imcolmap <- plot.im(x[[1L]], do.plot=FALSE, zlim=zlim, ...,
                          log=log, ribn=ribn)
    } else if(all(isfactor)) {
      x <- harmoniseLevels(x)
      ## determine common colour map based on factor levels
      imcolmap <- plot.im(x[[1L]], do.plot=FALSE, ..., ribn=ribn)
    } else warning("Could not determine a common colour map for these images",
                   call.=FALSE)
    ## plot ribbon?
    if(!ribbon) {
      ribadorn <- list()
    } else if(equal.scales) {
      ## colour ribbon will be aligned with objects in plot
      ribadorn <- list(adorn=imcolmap,
                       adorn.args=append(ribargs, list(labelmap=ribscale)))
      names(ribadorn)[1] <- paste("adorn", ribside, sep=".")
    } else {
      ## colour ribbon will be "free-floating"
      ## Determine plot arguments for ribbon
      vertical <- (ribside %in% c("right", "left"))
      scaleinfo <- if(!is.null(ribscale)) list(labelmap=ribscale) else list()
      sidecode <- sideCode(ribside)
      ribstuff <- c(list(x=imcolmap, main="", vertical=vertical),
                    ribargs,
                    scaleinfo,
                    list(side=sidecode))
      if (is.null(mar.panel)) 
        mar.panel <- c(2, 1, 1, 2)
      if (length(mar.panel) != 4) 
        mar.panel <- rep(mar.panel, 4)[1:4]
      if (is.null(ribmar)) {
        ribmar <- mar.panel/2
        newmar <- c(2, 0)
        switch(ribside,
               left   = { ribmar[c(2, 4)] <- newmar },
               right  = { ribmar[c(4, 2)] <- newmar },
               bottom = { ribmar[c(1, 3)] <- newmar },
               top    = { ribmar[c(3, 1)] <- newmar }
               )
      }
      ## bespoke function executed to plot colour ribbon
      do.ribbon <- function() {
        opa <- par(mar=ribmar)
        on.exit(par(opa))
        do.call(plot, ribstuff)
      }
      ## ribbon plot function encoded as 'adorn' argument
      ribadorn <- list(adorn=do.ribbon, adorn.size=ribwid)
      names(ribadorn)[1] <- paste("adorn", ribside, sep=".")
    }
    ##
    result <- do.call(plot.solist,
                      resolve.defaults(list(x=quote(x), plotcommand="image"),
                                       list(...),
                                       list(equal.scales=equal.scales,
                                            mar.panel=mar.panel,
                                            main=xname,
                                            do.plot=do.plot,
                                            col=imcolmap,
                                            zlim=zlim,
                                            ribbon=FALSE),
                                       ribadorn))
    return(invisible(result))
  }

  factorimage <- function(X, levels=NULL) {
    eval.im(factor(X, levels=levels))
  }

  plot.imlist
})

image.imlist <- image.listof <-
  function(x, ..., equal.ribbon = FALSE, equal.scales=FALSE, ribmar=NULL) {
    plc <- resolve.1.default(list(plotcommand="image"), list(...))
    if(list(plc) %in% list("image", "plot", image, plot)) {
      out <- do.call(plot.imlist,
                     resolve.defaults(
                       list(x=quote(x), plotcommand="image"),
                       list(...),
                       list(equal.ribbon=equal.ribbon,
                            equal.scales=equal.scales,
                            ribmar=ribmar)))
    } else {
      out <- plot.solist(x, ..., equal.scales=equal.scales, ribmar=ribmar)
    }
    return(invisible(out))
  }

recognise.plotcommand.type <- function(s) {
  if(is.name(s)) s <- as.character(s)
  if(is.character(s) && length(s) == 1) {
    for(a in c("plot", "contour", "image", "persp")) {
      ## recognise name of generic
      if(s == a) return(a)
      ## recognise name of a method
      if(nchar(s) > nchar(a) && substr(s, 1, nchar(a)+1) == paste0(a, "."))
        return(a)
    }
  }
  if(is.function(s)) {
    ## recognise the generic
    m <- match(list(s), list(plot, contour, image, persp))
    if(!is.na(m))
      return(c("plot", "contour", "image", "persp")[m])
  }
  return("unknown")
}

