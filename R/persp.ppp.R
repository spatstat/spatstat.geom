#'  persp.ppp.R
#'
#'  Perspective plot for marked point pattern
#'
#'  Copyright (C) Adrian Baddeley 2024
#'  GPL Public Licence >= 2.0
#'
#' $Revision: 1.3 $ $Date: 2024/09/07 01:34:04 $

persp.ppp <- local({
  
  persp.ppp <- function(x, ..., main, grid=TRUE, ngrid=10,
                        col.grid="grey", col.base="white",
                        win.args=list(),
                        spike.args=list(),
                        neg.args=list(), which.marks=1,
                        zlab=NULL, zlim=NULL,
                        zadjust=1) {
    if(missing(main)) main <- short.deparse(substitute(x))
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
    #' rescale marks to a scale commensurate with window
    #'  (to achieve appropriate default scale in persp.default)
    maxmark <- max(abs(marx))
    if(maxmark > .Machine$double.eps) {
      marx <- marx * (max(sidelengths(R))/maxmark)
    }
    #' set up perspective transformation and plot horizontal plane
    Rplus <- grow.rectangle(R, fraction=1/(2*ngrid))
    Z <- as.im(0, W=Rplus, dimyx=rev(ngrid)+1)
    if(is.null(zlim)) zlim <- range(marx, 0)
    check.range(zlim)
    col.grid.used <- if(grid && (zlim[1] >= 0)) col.grid else NA
    argh <- resolve.defaults(list(x=Z, main=main,
                                  border=col.grid.used,
                                  col=col.base),
                             dotargs,
                             list(axes=FALSE, box=FALSE,
                                  zlim=zlim, zlab=zlab,
                                  #' do not independently rescale x & y
                                  scale=FALSE,
                                  #' expand=0.1 is default in persp.default
                                  expand=zadjust * 0.1))
    M <- do.call.matched(persp.im, argh,
                         funargs=graphicsPars("persp"))
    #' create spikes
    S <- xyzsegmentdata(x$x, x$y, 0,
                        x$x, x$y, marx)
    if(grid) {
      if(zlim[1] < 0) {
        #' first draw downward spikes
        downward <- (marx < 0)
        if(any(downward)) {
          SD <- S[downward, , drop=FALSE]
          spectiveSegments(SD, neg.args, spike.args, ..., M=M)
          S <- S[!downward, , drop=FALSE]
        }
        #' plot baseline grid 
        spectiveFlatGrid(R, ngrid, M, col=col.grid)
      }
    }
    if(!is.rectangle(W)) {
      #' plot window
      spectiveFlatPolygons(W, M, win.args, dotargs)
    }
    #' draw upward spikes
    if(nrow(S) > 0) {
      spectiveSegments(S, spike.args, dotargs, M=M)
    }
    #'   
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

  persp.ppp
})


