##
## solist.R
##
## Methods for class `solist' (spatial object list)
##
##      and related classes 'anylist', 'ppplist', 'imlist', 'linimlist'
##
## plot.solist is defined in plot.solist.R
##
## $Revision: 1.36 $ $Date: 2025/09/28 03:57:38 $

anylist <- function(...) {
  x <- list(...)
  class(x) <- unique(c("anylist", "listof", class(x)))
  return(x)
}

print.anylist <- function (x, ...) {
  ll <- length(x)
  if(ll == 0) {
    splat("(Zero length list)")
    return(invisible(NULL))
  }
  nn <- names(x)
  if (length(nn) != ll) 
    nn <- paste("Component", seq.int(ll))
  spaceok <- waxlyrical('space')
  for (i in seq_len(ll)) {
    splat(paste0(nn[i], ":"))
    print(x[[i]], ...)
    if(spaceok && i < ll) cat("\n")
  }
  return(invisible(NULL))
}

as.anylist <- function(x) {
  if(inherits(x, "anylist")) return(x)
  if(!is.list(x))
    x <- list(x)
  class(x) <- unique(c("anylist", "listof", class(x)))
  return(x)
}
  
"[.anylist" <- function(x, i, ...) {
  cl <- oldClass(x)
  ## invoke list method
  y <- NextMethod("[")
  if(length(y) == 0) return(list())
  class(y) <- cl
  return(y)
}

"[<-.anylist" <- function(x, i, value) {
  as.anylist(NextMethod("[<-"))
}

summary.anylist <- function(object, ...) {
  as.anylist(lapply(object, summary, ...))
}

## .................... solist .............................

is.sob <- local({
  ## test whether x is a spatial object suitable for solist
  sobjectclasses <- c("ppp", "psp", "im", "owin", 
                      "quad", "tess", "msr",
                      "quadratcount", "quadrattest", 
                      "layered",
                      "funxy", "distfun", "nnfun", 
                      "lpp", "linnet", "linfun", "lintess",  
                      "influence.ppm", "leverage.ppm")
  ## Note 'linim' inherits 'im'
  ##      'dfbetas.ppm' inherits 'msr'
  ##      diagram objects typically inherit 'ppp'
  is.sob <- function(x) { inherits(x, what=sobjectclasses) }
  is.sob
})
  
solist <- function(..., check=TRUE, promote=TRUE, demote=FALSE, .NameBase) {
  stuff <- list(...)
  if(length(stuff)) {
    if(!missing(.NameBase) && !any(nzchar(names(stuff))))
      names(stuff) <- paste(.NameBase, seq_along(stuff))
    if(check || promote)
      stuff <- coerceNAtoObject(stuff)
    if((check || demote) && !all(sapply(stuff, is.sob))) {
      if(demote) return(as.anylist(stuff)) else
      stop("Some arguments of solist() are not 2D spatial objects")
    }
  }
  class(stuff) <- unique(c("solist", "anylist", "listof", class(stuff)))
  if(promote && length(stuff)) {
    if(all(sapply(stuff, is.ppp))) {
      class(stuff) <- unique(c("ppplist", class(stuff)))
    } else if(all(sapply(stuff, is.im))) {
      class(stuff) <- unique(c("imlist", class(stuff)))
      if(all(sapply(stuff, is.linim)))
        class(stuff) <- unique(c("linimlist", class(stuff)))
    }
  }
  return(stuff)
}

coerceNAtoObject <- function(x, cls=NULL) {
  ## coerce vanilla NA entries x[[i]]
  ## to a spatial NA object (or atomic NA value), if class is unambiguous
  if(any(isvna <- sapply(x, identical, y=NA))) {
    if(is.null(cls)) 
      cls <- unique(sapply(x[!isvna], classIgnoringNA, first=TRUE))
    if(length(cls) == 1) {
      #' coerce NA entries to NA object of this class
      nao <- NAobject(cls)
      #' unless ...
      if(all(sapply(x, is.atomic)) && all(lengths(x) == 1)) {
        ## each entry of x is a single atomic value
        nao <- switch(cls,
                      logical = NA,
                      integer = NA_integer_,
                      numeric = NA_real_,
                      character = NA_character_,
                      complex = NA_complex_,
                      nao)
      }
      x[isvna] <- rep(list(nao), sum(isvna))
    }
  }
  return(x)
}

as.solist <- function(x, ...) {
  if(inherits(x, "solist") && length(list(...)) == 0) {
    #' wipe superfluous info
    if(inherits(x, "ppplist")) 
      attributes(x)[c("fsplit", "fgroup")] <- NULL
    class(x) <- c("solist", "anylist", "listof")
    return(x)
  }
  #' needs to be enclosed in list() ?
  if(!is.list(x) || (is.sob(x) && !inherits(x, "layered")))
    x <- list(x)
  return(do.call(solist, append(x, list(...))))
}

is.solist <- function(x) inherits(x, "solist")

print.solist <- function (x, ...) {
  what <- if(inherits(x, "ppplist")) "point patterns" else
          if(inherits(x, "linimlist")) "pixel images on a network" else
          if(inherits(x, "imlist")) "pixel images" else "spatial objects"
  splat(paste("List of", what))
  parbreak()
  NextMethod("print")
}


"[.solist" <- function(x, i, ...) {
  cl <- oldClass(x)
  if(!missing(i) && is.owin(i)) {
    ## spatial subset
    y <- lapply(unclass(x), "[", i=i, ...)
  } else {
    ## invoke list method
    y <- NextMethod("[")
  }
  if(length(y) == 0) return(list())
  class(y) <- cl
  return(y)
}
  
"[<-.solist" <- function(x, i, value) {
  ## invoke list method
  y <- NextMethod("[<-")
  ## check again
  return(do.call(solist, y))
}
  
"[[<-.solist" <- function(x, i, value) {
  ## invoke list method
  y <- NextMethod("[[<-")
  ## check again
  return(do.call(solist, y))
}

is.na.solist <- function(x) { sapply(x, is.NAobject) }

summary.solist <- function(object, ...) {
  x <- lapply(object, summary, ...)
  attr(x, "otype") <-
    if(inherits(object, "ppplist")) "ppp" else
    if(inherits(object, "linimlist")) "im" else ""
    if(inherits(object, "imlist")) "im" else ""
  class(x) <- c("summary.solist", "anylist")
  x
}

print.summary.solist <- function(x, ...) {
  what <- switch(attr(x, "otype"),
                 ppp="point patterns",
                 im="pixel images",
                 "spatial objects")
  splat("Summary of", length(x), what)
  parbreak()
  NextMethod("print")
}

as.layered.solist <- function(X) {
  layered(LayerList=X)
}

#'  ----- ppplist and imlist methods ----------------------------

as.data.frame.ppplist <- local({

  rnf <- function(x) {
    n <- nrow(x)
    if(n > 0) 
      row.names(x) <- paste("Point", seq_len(n))
    return(x)
  }
  
  as.data.frame.ppplist <- function(x, row.names = NULL, ...) {
    y <- lapply(lapply(x, as.data.frame.ppp), rnf)
    if(is.null(row.names)) {
      #' work around a quirk of 'rbind'
      singleton <- (sapply(y, nrow) == 1)
      if(any(singleton))
        names(y)[singleton] <- paste0(names(x)[singleton], ".Point 1")
    }
    z <- do.call(rbind, y)
    if(!is.null(row.names)) 
      row.names(z) <- row.names
    return(z)
  }

  as.data.frame.ppplist
  
})

#'  ----- ppplist and imlist ----------------------------
#'  for efficiency only

as.ppplist <- function(x, check=TRUE) {
  if(check) {
    x <- as.solist(x, promote=TRUE, check=TRUE)
    if(!inherits(x, "ppplist"))
      stop("some entries are not point patterns")
  }
  class(x) <- unique(c("ppplist", "solist", "anylist", "listof", class(x)))
  return(x)
}

is.ppplist <- function(x) inherits(x, "ppplist")

as.imlist <- function(x, check=TRUE) {
  if(check) {
    x <- as.solist(x, promote=TRUE, check=TRUE)
    if(!inherits(x, "imlist"))
      stop("some entries are not images")
  } else {
    #' just apply required classes, in required order
    reqd <- c("imlist", "solist", "anylist", "listof")
    class(x) <- unique(c(reqd, class(x)))
  }
  return(x)
}

is.imlist <- function(x) inherits(x, "imlist")

as.linimlist <- function(x, check=TRUE) {
  x <- as.imlist(x, check=check)
  if(check) {
    if(!all(sapply(x, is.linim)))
      stop("All entries must be pixel images on a network")
  } 
  class(x) <- unique(c("linimlist", class(x)))
  return(x)
}


# --------------- counterparts of 'lapply' --------------------

anylapply <- function(X, FUN, ...) {
  ok <- !sapply(X, is.NAobject)
  if(all(ok)) {
    v <- lapply(X, FUN, ...)
  } else {
    v <- rep(list(NA), length(X))
    v[ok] <- lapply(X[ok], FUN, ...)
    v <- coerceNAtoObject(v)
  }
  return(as.anylist(v))
}

solapply <- function(X, FUN, ..., check=TRUE, promote=TRUE, demote=FALSE) {
  ok <- !sapply(X, is.NAobject)
  if(all(ok)) {
    v <- lapply(X, FUN, ...)
  } else {
    v <- rep(list(NA), length(X))
    v[ok] <- lapply(X[ok], FUN, ...)
    v <- coerceNAtoObject(v)
  }
  u <- as.solist(v, check=check, promote=promote, demote=demote)
  return(u)
}

expandSpecialLists <- function(x, special="solist") {
  ## x is a list which may include entries which are lists, of class 'special'
  ## unlist these entries only
  hit <- sapply(x, inherits, what=special) 
  if(!any(hit)) return(x)
  # wrap each *non*-special entry in list()
  x[!hit] <- lapply(x[!hit], list)
  # now strip one layer of list() from all entries
  return(unlist(x, recursive=FALSE))
}
