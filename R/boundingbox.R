##
## boundingbox.R
##
## $Revision: 1.12 $ $Date: 2024/02/04 08:04:51 $

# bounding.box <- function(...) {
#  .Deprecated("boundingbox", "spatstat")
#  boundingbox(...)
# }

boundingbox <- function(...) {
  ## remove any NULL arguments
  arglist <- list(...)
  if(any(isnull <- sapply(arglist, is.null))) {
    if(length(arglist[!isnull]))
       return(do.call(boundingbox, arglist[!isnull]))
    stop("No non-null arguments given.\n")
  }
  UseMethod("boundingbox")
}

boundingbox.solist <- function(...) {
  argh <- list(...)
  issl <- sapply(argh, inherits, what="solist")
  yarg <- c(do.call(c, argh[issl]), argh[!issl])
  do.call(bbEngine, yarg)
}

boundingbox.ppp  <-
boundingbox.psp  <-
boundingbox.owin <-
boundingbox.list <-
boundingbox.linnet <-
boundingbox.lpp <-
boundingbox.im   <- function(...) {
   bbEngine(...)
}

recognise.spatstat.type <- local({

  knowntypes <- c("ppp","psp","owin","im", "lpp", "linnet")

  function(x) {
    for(kt in knowntypes)
      if(inherits(x, kt)) return(kt)
    if(is.list(x) && checkfields(x, c("x", "y"))
       && is.numeric(x$x) && is.numeric(x$y) &&
       is.vector(x$x) && is.vector(x$y) && length(x$x) == length(x$y))
        return("listxy")
    aso <- try(as.owin(x), silent=TRUE)
    if(!inherits(aso, "try-error")) return("as.owin")
    return("unknown")
  }
})

bbEngine <- local({

  bb.listxy <- function(X) owinInternalRect(range(X$x), range(X$y))

  bb.linnet <- function(X) boundingbox(vertices(X))

  bb.lpp <- function(X) boundingbox(as.ppp(X))
  
  bbEngine <- function(...) {
    wins <- list(...)
    ## first detect any numeric vector arguments
    if(any(isnumvec <- unlist(lapply(wins, is.vector)) &
           unlist(lapply(wins, is.numeric)))) {
      ## invoke default method on these arguments
      bb <- do.call(boundingbox, wins[isnumvec])
      ## repack
      wins <- append(wins[!isnumvec], list(bb))
    }
    if(length(wins) > 1) {
      ## multiple arguments -- compute bounding box for each argument.
      objtype <- unlist(lapply(wins, recognise.spatstat.type))
      nbad <- sum(objtype == "unknown")
      if(nbad > 0) {
        whinge <- paste("Function boundingbox called with",
                        nbad,"unrecognised",
                        ngettext(nbad,"argument","arguments"))
        stop(whinge, call.=FALSE)
      }
      if(any(isppp <- (objtype == "ppp"))) 
        wins[isppp] <- lapply(wins[isppp], boundingbox)
      if(any(islistxy <- (objtype == "listxy")))
        wins[islistxy] <- lapply(wins[islistxy], bb.listxy)
      if(any(isnet <- (objtype == "linnet")))
        wins[isnet] <- lapply(wins[isnet], bb.linnet)
      if(any(islpp <- (objtype == "lpp")))
        wins[islpp] <- lapply(wins[islpp], bb.lpp)
      ## then convert all windows to owin
      wins <- lapply(wins, as.owin)
      ## then take bounding box of each window
      boxes <- lapply(wins, boundingbox)
      ## discard NULL values
      isnull <- unlist(lapply(boxes, is.null))
      boxes <- boxes[!isnull]
      ## take bounding box of these boxes
      xrange <- range(unlist(lapply(boxes, getElement, name="xrange")))
      yrange <- range(unlist(lapply(boxes, getElement, name="yrange")))
      W <- owinInternalRect(xrange, yrange)
      ## If all of the windows have a common unit name, give
      ## that unit name to the bounding box.
      youse <- unique(t(sapply(boxes,unitname)))
      if(nrow(youse)==1) {
        ute <- unlist(youse[1L,])
        unitname(W) <- ute
      }
      return(W)
    }

    ## single argument
    w <- wins[[1L]]
    if(is.null(w))
      return(NULL)
    
    wtype <- recognise.spatstat.type(w)
    ## point pattern?
    if(wtype == "ppp")
      return(boundingbox(coords(w)))
    
    ## line segment pattern?
    if(wtype == "psp")
      return(boundingbox(endpoints.psp(w)))
    
    ## list(x,y)
    if(wtype == "listxy")
      return(bb.listxy(w))

    if(wtype == "linnet")
      w <- return(bb.linnet(w))

    if(wtype == "lpp")
      w <- return(bb.lpp(w))
    
    ## convert to window
    w <- as.owin(w)

    ## determine a tight bounding box for the window w
    switch(w$type,
           rectangle = {
             return(w)
           },
           polygonal = {
             bdry <- w$bdry
             if(length(bdry) == 0)
               return(NULL)
             xr <- range(unlist(lapply(bdry, rangeofx)))
             yr <- range(unlist(lapply(bdry, rangeofy)))
             return(owinInternalRect(xr, yr, unitname=unitname(w)))
           },
           mask = {
             m <- w$m
             x <- rasterx.mask(w)
             y <- rastery.mask(w)
             xr <- range(x[m]) + c(-1,1) * w$xstep/2
             yr <- range(y[m]) + c(-1,1) * w$ystep/2
             return(owinInternalRect(xr, yr, unitname=unitname(w)))
           },
           stop("unrecognised window type", w$type)
           )
  }

  rangeofx <- function(a) range(a$x)
  rangeofy <- function(a) range(a$y)
  
  bbEngine
})


boundingbox.default <- local({

  bb.listxy <- function(X) owinInternalRect(range(X$x), range(X$y))

  boundingbox.default <- function(...) {
    arglist <- list(...)
    bb <- NULL
    if(length(arglist) == 0)
      return(bb)
    ## handle numeric vector arguments
    if(any(isnumvec <- unlist(lapply(arglist, is.vector)) &
           unlist(lapply(arglist, is.numeric)))) {
      nvec <- sum(isnumvec)
      if(nvec != 2)
        stop(paste("boundingbox.default expects 2 numeric vectors:",
                   nvec, "were supplied"),
             call.=FALSE)
      vecs <- arglist[isnumvec]
      x <- vecs[[1L]]
      y <- vecs[[2L]]
      bb <- if(length(x) == length(y)) owinInternalRect(range(x), range(y)) else NULL
      arglist <- arglist[!isnumvec]
    }
    if(length(arglist) == 0)
      return(bb)
    ## other objects are present
    objtype <- unlist(lapply(arglist, recognise.spatstat.type))
    ## Unrecognised?
    nbad <- sum(objtype == "unknown")
    if(nbad > 0) {
      whinge <- paste("Function boundingbox called with",
                      nbad,"unrecognised",
                      ngettext(nbad,"argument","arguments"))
      stop(whinge, call.=FALSE)
    }
    if(any(aso <- (objtype == "as.owin"))) {
      ## promote objects to owin (to avoid infinite recursion!)
      arglist[aso] <- lapply(arglist[aso], as.owin)
    }
    if(any(lxy <- (objtype == "listxy"))) {
      ## handle list(x,y) objects 
      arglist[lxy] <- lapply(arglist[lxy], bb.listxy)
    }
    result <- do.call(boundingbox,
                      if(is.null(bb)) arglist else append(list(bb), arglist))
    return(result)
  }

  boundingbox.default
})


