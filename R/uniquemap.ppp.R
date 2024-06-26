#'
#'   uniquemap.ppp.R
#'
#'   Methods for 'uniquemap' for classes ppp, lpp, ppx
#' 
#'   Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2019
#'   Licence: GNU Public Licence >= 2
#'
#'   $Revision: 1.1 $  $Date: 2024/06/26 08:26:30 $

uniquemap.ppp <- function(x) {
  n <- npoints(x)
  seqn <- seq_len(n)
  if(n <= 1) return(seqn)
  marx <- marks(x)
  switch(markformat(marx),
         none = {
           useC <- TRUE
         },
         vector = {
           #' convert to integers if possible
           if(is.integer(marx) || is.factor(marx)) {
             marx <- as.integer(marx)
             useC <- TRUE
           } else {
             um <- unique(marx)
             if(length(um) <= 2^30) {
               marx <- match(marx, um)
               useC <- TRUE
             } else {
               useC <- FALSE
             }
           }
         },
         {
           useC <- FALSE
         })

  if(!useC) {
    #' first find duplicated spatial coordinates
    u <- uniquemap(unmark(x))
    #' add marks
    df <- cbind(data.frame(ind=seqn, uni=u), as.data.frame(marx))
    bb <- split(df, factor(u))
    #' consider each set of duplicated locations
    for(b in bb) {
      #' find unique rows of marks, as a list
      mrows <- lapply(seq_len(nrow(b)), function(i) b[i, -(1:2)])
      um <- unique(mrows)
      #' match other rows to them
      ma <- match(mrows, um)
      #' map to original index
      u[b$ind] <- b$ind[ma]
    }
    return(u)
  }

  #' unmarked or integer/factor marked
  xx <- x$x
  yy <- x$y
  o <- order(xx, seqn)

  if(is.null(marx)) {
    umap <- .C(SG_uniqmapxy,
               n=as.integer(n),
               x=as.double(xx[o]),
               y=as.double(yy[o]),
               uniqmap=as.integer(integer(n)),
               PACKAGE="spatstat.geom")$uniqmap
  } else {
    #' marks are (converted to) integers
    umap <- .C(SG_uniqmap2M,
               n=as.integer(n),
               x=as.double(xx[o]),
               y=as.double(yy[o]),
               marks=as.integer(marx[o]),
               uniqmap=as.integer(integer(n)),
               PACKAGE="spatstat.geom")$uniqmap
  }
  nodup <- (umap == 0)
  umap[nodup] <- which(nodup)
  result <- integer(n)
  result[o] <- o[umap]
  return(result)
}

uniquemap.lpp <- function(x) {
  n <- npoints(x)
  if(n <= 1 || !anyDuplicated(as.ppp(x))) return(seq_len(n))
  result <- uniquemap(as.data.frame(x))
  return(result)
}

uniquemap.ppx <- function(x) {
  uniquemap(as.data.frame(x))
}

