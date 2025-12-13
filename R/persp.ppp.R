#'  persp.ppp.R
#'
#'  Perspective plot for marked point pattern
#'
#'  Copyright (C) Adrian Baddeley 2024
#'  GPL Public Licence >= 2.0
#'
#' $Revision: 1.12 $ $Date: 2025/12/13 07:52:15 $

persp.ppp <- local({
  
  persp.ppp <- function(x, ..., main, type=c("l", "b"), grid=TRUE, ngrid=10,
                        col.grid="grey", col.base="white",
                        win.args=list(),
                        spike.args=list(), 
                        neg.args=list(), point.args=list(), which.marks=1,
                        zlab=NULL, zlim=NULL,
                        zadjust=1, show.window=TRUE,
                        legend=TRUE, legendpos="bottomleft",
                        leg.args=list(lwd=4),
                        leg.col=c("black", "orange")) {
    if(missing(main)) main <- short.deparse(substitute(x))
    type <- match.arg(type)
    W <- Window(x)
    R <- Frame(x)
    marx <- marks(x)
    dotargs <- list(...)
    #' ensure numeric marks
    mf <- markformat(marx)
    switch(mf,
           none = {
             stop("point pattern must have marks", call.=FALSE)
           },
           vector = {
             if(is.null(zlab)) zlab <- "mark"
           },
           dataframe = {
             if(missing(which.marks)) {
               marx <- numeric.columns(marx)
             }
             cn <- colnames(marx)
             stopifnot(length(which.marks) == 1)
             if(is.character(which.marks)) {
               k <- match(which.marks, cn)
               if(is.na(k))
                 stop(paste("unrecognised selection of mark column",
                            sQuote(which.marks)),
                      call.=FALSE)
               which.marks <- k
             }
             if(missing(zlab)) zlab <- colnames(marx)[which.marks]
             marx <- marx[,which.marks]
           },
           stop("marks should be a vector or a data frame", call.=FALSE)
           )
    marx <- as.numeric(marx)
    if(is.null(zlim)) zlim <- range(marx, 0)
    check.range(zlim)
    #' rescale marks to a scale commensurate with window
    #'  (to achieve appropriate default scale in persp.default)
    maxmark <- max(abs(marx))
    if(maxmark > .Machine$double.eps) {
      scal <- max(sidelengths(R))/maxmark
      scaled.marx <- scal * marx
      scaled.zlim <- scal * zlim
    } else {
      scaled.marx <- marx
      scaled.zlim <- zlim
    }
    #' Set up objects to be plotted in perspective
    Rplus <- grow.rectangle(R, fraction=1/(2*ngrid))
    #' spikes
    S <- xyzsegmentdata(x$x, x$y, 0,
                        x$x, x$y, scaled.marx)
    #' bubbles
    if(type == "b")
      P <- data.frame(x=x$x, y=x$y, z=scaled.marx)

    #' Assemble arguments for persp.default
    col.grid.used <- if(grid && (zlim[1] >= 0)) col.grid else NA
    if(!is.na(k <- match("adj.main", names(dotargs))))
      names(dotargs)[k] <- "adj"
    if(is.im(col.base)) {
      #' Horizontal plane will be painted by a colour image
      if(!is.subset.owin(Window(x), Window(col.base)))
        Window(x) <- boundingbox(Window(x), Window(col.base))
      #' base plane height
      Z <- as.im(0, W=col.base)
      BaseInfo <- list(colin=col.base)
    } else if(length(col.base) == 1 && is.colour(col.base)) {
      #' Usual case: horizontal plane will be painted a single colour
      #' Base plane image
      Z <- as.im(0, W=Rplus, dimyx=rev(ngrid)+1)
      #' Draw grid lines by setting 'border' argument of persp.default
      #' provided the function has no negative values
      border <- if(grid && (zlim[1] >= 0)) col.grid else NA
      BaseInfo <- list(col=col.base, border=border)
    } else {
      stop("Argument col.base should be a single colour or a pixel image",
           call.=FALSE)
    }
    argh <- resolve.defaults(list(x=quote(Z), main=main),
                             BaseInfo,
                             dotargs,
                             list(axes=FALSE, box=FALSE,
                                  zlim=scaled.zlim, zlab=zlab,
                                  #' do not independently rescale x & y
                                  scale=FALSE,
                                  #' expand=0.1 is default in persp.default
                                  expand=zadjust * 0.1))
    #' Start perspective plot; plot horizontal plane
    M <- do.call.matched(persp.im, argh,
                         extrargs=graphicsPars("persp"))
    
    #' Start drawing objects
    if(grid) {
      if(scaled.zlim[1] < 0) {
        #' first draw downward spikes
        downward <- (scaled.marx < 0)
        if(any(downward)) {
          SD <- S[downward, , drop=FALSE]
          spectiveSegments(SD, neg.args, spike.args, dotargs, M=M)
          S <- S[!downward, , drop=FALSE]
          if(type == "b") {
            PD <- P[downward, , drop=FALSE]
            spectivePoints(PD, point.args, dotargs, M=M)
            P <- P[!downward, , drop=FALSE]
          }
        }
        #' plot baseline grid 
        spectiveFlatGrid(R, ngrid, M, col=col.grid)
      }
    }
    if(isTRUE(show.window) && !is.rectangle(W)) {
      #' plot window
      spectiveFlatPolygons(W, M, win.args, dotargs)
    }
    #' draw upward spikes
    if(nrow(S) > 0) {
      spectiveSegments(S, spike.args, dotargs, M=M)
      if(type == "b")
        spectivePoints(P, point.args, dotargs, M=M)
    }
    #'
    if(legend) {
      #' draw a reference scale as another spike
      #' determine spike position
      if(is.character(legendpos)) {
        legendpos <- match.arg(legendpos, c("bottomleft", "bottomright",
                                            "topleft", "topright",
                                            "bottom", "left", "top", "right"))
        B <- Frame(x)
        xr <- B$xrange
        yr <- B$yrange
        legxy <- switch(legendpos,
                        bottomleft = c(xr[1], yr[1]),
                        bottomright = c(xr[2], yr[1]),
                        topleft = c(xr[1], yr[2]),
                        topright = c(xr[2], yr[2]),
                        bottom = c(mean(xr), yr[1]),
                        left = c(xr[1], mean(yr)),
                        top = c(mean(xr), yr[2]),
                        right = c(xr[2], mean(yr)))
      } else legxy <- ensure2vector(unlist(legendpos))
      #' determine tickmarks
      tix <- unique(sort(c(zlim, prettyinside(zlim))))
      ntix <- length(tix)
      scaled.tix <- scal * tix
      tixseg <- xyzsegmentdata(legxy[1], legxy[2], scaled.tix[-ntix],
                               legxy[1], legxy[2], scaled.tix[-1])
      tixcol <- rep(leg.col, ntix)[1:ntix]
      #' inclination of vertical spike relative to projection plane
      a <- trans3dz(legxy[1], legxy[2], scal*zlim, M)
      spikeangle <- with(a, atan2(diff(y), diff(x))) * 180/pi - 90
      #' draw zebra segments
      spectiveSegments(tixseg, list(col=tixcol), leg.args, M=M)
      #' draw tickmark values
      spectiveText(legxy[1], legxy[2], scaled.tix[-c(1,ntix)],
                   labels=tix[-c(1,ntix)], pos=4, M=M, srt=spikeangle)
      #' add text for z axis?
      if(!isTRUE(argh$box) && nchar(zlab) > 0) {
        #' persp.default(box=FALSE) suppresses axis labels,
        #' so draw a vertical axis label now
        spectiveText(legxy[1], legxy[2], scal * zlim[2],
                     labels=zlab, M=M,
                     pos=2, offset=0.75, srt=spikeangle + 90)
      }
    }
    invisible(M)
  }

  xyzsegmentdata <- function(x0, y0, z0, x1, y1, z1) {
    data.frame(x0=x0, y0=y0, z0=z0, x1=x1, y1=y1, z1=z1)
  }
  
  trans3dz <- function(x,y,z,pmat) {
    tr <- cbind(x, y, z, 1) %*% pmat
    list(x = tr[, 1]/tr[, 4],
         y = tr[, 2]/tr[, 4],
         z = tr[, 3]/tr[, 4])
  }

  spectiveFlatGrid <- function(B, ngrid, M, ...) {
    ## arguments ... should be lists of parameters
    B <- Frame(B)
    xr <- B$xrange
    yr <- B$yrange
    ngrid <- ensure2vector(ngrid)
    xx <- seq(xr[1], xr[2], length.out=ngrid[1]+1)
    yy <- seq(yr[1], yr[2], length.out=ngrid[2]+1)
    horiz <- xyzsegmentdata(xr[1], yy,    0, xr[2], yy,    0)
    vert  <- xyzsegmentdata(xx,    yr[1], 0, xx,    yr[2], 0)
    spectiveSegments(horiz, ..., M=M)
    spectiveSegments(vert,  ..., M=M)
    invisible(NULL)
  }

  spectiveFlatPolygons <- function(W, M, ...) {
    ## arguments ... should be lists of parameters
    Wbdry <- as.polygonal(W)$bdry
    Pbdry <- lapply(Wbdry,
                    function(p, M) {
                      as.list(trans3dz(p$x, p$y, 0, M))[c("x","y")]
                    }, M=M)
    P <- owin(poly=Pbdry, check=FALSE, fix=FALSE)
    do.call(plot.owin,
            resolve.defaults(list(quote(P)),
                             ..., 
                             list(add=TRUE)))
  }

  spectiveSegments <- function(df, ..., M) {
    ## arguments ... should be lists of parameters
    a0 <- with(df, trans3dz(x0, y0, z0, M))
    a1 <- with(df, trans3dz(x1, y1, z1, M))
    do.call.matched(segments,
                    resolve.defaults(
                      list(x0=a0$x,
                           y0=a0$y,
                           x1=a1$x,
                           y1=a1$y),
                      ...))
    invisible(NULL)
  }

  spectivePoints <- function(df, ..., M) {
    ## arguments ... should be lists of parameters
    p <- with(df, trans3dz(x, y, z, M))
    do.call.matched(points.default,
                    resolve.defaults(
                      list(x=p$x, y=p$y),
                      ...),
                    extrargs=graphicsPars("points"))
  }
  
  spectiveText <- function(x,y,z, ..., M) {
    p <- trans3dz(x, y, z, M)
    text(p$x, p$y, ...)
  }

  persp.ppp
})


