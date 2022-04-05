#'
#'   bufftess.R
#'
#'   Buffer (Distance) Tessellation
#'
#'   $Revision: 1.3 $ $Date: 2022/03/31 09:05:10 $
#'
#'   Copyright (c) 2021 Adrian Baddeley, Ege Rubak and Rolf Turner
#'   GNU Public Licence >= 2.0

bufftess <- function(X, breaks, W=Window(X), ..., polygonal=TRUE) {
  breaks <- as.numeric(breaks)
  Wgiven <- !missing(W)
  if(!polygonal || length(breaks) == 1L) {
    ## Determine break points from distance values in distmap
    D <- distmap(X, ...)
    if(Wgiven) D <- D[W, drop=FALSE]
    drange <- c(0, range(D))
    breaks <- exactCutBreaks(drange, breaks)
    ## ensure break points are nonzero
    breaks <- unique(pmax(0, breaks))
  }
  if(!polygonal) {
    ## pixel image tessellation 
    G <- cut(x=D, breaks=breaks, ...)
    Y <- tess(image=G)
    attr(Y, "breaks") <- breaks
    return(Y)
  } else {
    ## polygonal tiles tessellation 
    W <- as.polygonal(W)
    dbig <- diameter(Frame(W))
    nbreaks <- length(breaks)
    nbands <- nbreaks - 1L
    Ytiles <- vector(mode="list", length=nbands)
    for(ibreak in seq_len(nbreaks)) {
      d <- breaks[ibreak]
      #' dilation
      if(d > dbig) {
        B <- W
      } else if(d > 0) {
        B <- dilation(X, d, polygonal=TRUE)
        B <- intersect.owin(B, W)
      } else {
        B <- NULL
      }
      #' set difference
      if(ibreak == 1L) {
        Bprev <- Bmin <- B
      } else {
        iband <- ibreak - 1L
        Ytiles[[iband]] <- setminus.owin(B, Bprev)
        Bprev <- B
      }
    }
    names(Ytiles) <- levels(cut(breaks, breaks, ...))
    Wfinal <- rescue.rectangle(setminus.owin(B, Bmin))
    Y <- tess(tiles=Ytiles, window=Wfinal)
    attr(Y, "breaks") <- breaks
    return(Y)
  }
}

