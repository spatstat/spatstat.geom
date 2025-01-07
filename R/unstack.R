#'
#'   unstack.R
#'
#'   Methods for generic 'unstack'
#' 
#'   $Revision: 1.7 $  $Date: 2025/01/07 07:12:17 $

unstack.ppp <- unstack.psp <-
    unstack.tess <- function(x, ...) {
  trap.extra.arguments(...)
  marx <- marks(x, drop=FALSE)
  d <- dim(marx)
  if(is.null(d)) return(solist(x))
  y <- rep(list(unmark(x)), d[2])
  for(j in seq_along(y))
    marks(y[[j]]) <- marx[,j,drop=FALSE]
  names(y) <- colnames(marx)
  return(as.solist(y))
}


unstackFilter <- function(x) {
  ## deal with a whole swag of classes that do not need to be unstacked
  nonvectorclasses <- c("im", "owin", "quad", 
                        "quadratcount", "quadrattest", 
                        "funxy", "distfun", "nnfun", 
                        "linnet", "linfun", 
                        "influence.ppm", "leverage.ppm")
  y <- if(inherits(x, nonvectorclasses)) solist(x) else unstack(x)
  return(y)
}

unstack.solist <- function(x, ...) {
  trap.extra.arguments(...)
  y <- lapply(x, unstackFilter)
  z <- as.solist(unlist(y, recursive=FALSE))
  return(z)
}

unstack.layered <- function(x, ...) {
  trap.extra.arguments(...)
  y <- lapply(x, unstackFilter)
  ny <- lengths(y)
  nx <- length(ny)
  if(all(ny == 1) || nx == 0) return(solist(x))
  pax <- layerplotargs(x)
  pay <- rep(pax, times=ny)
  z <- unlist(y, recursive=FALSE)
  z <- layered(LayerList=z, plotargs=pay)
  return(z)
}


  




  


  
