#'
#'   factors.R
#'
#'  Tools for manipulating factors and factor-valued things
#'
#'  $Revision: 1.7 $  $Date: 2023/01/05 23:33:47 $

relevel.im <- function(x, ref, ...) {
  if(x$type != "factor")
    stop("Only valid for factor-valued images")
  x[] <- relevel(x[], ref, ...)
  return(x)
}

relevel.ppp <- relevel.ppx <- function(x, ref, ...) {
  stopifnot(is.multitype(x))
  marks(x) <- relevel(marks(x), ref, ...)
  return(x)
}

mergeLevels <- function(.f, ...) {
  if(is.im(.f)) {
    aa <- mergeLevels(.f[], ...)
    .f[] <- aa
    return(.f)
  }
  if(is.multitype(.f)) {
    marks(.f) <- mergeLevels(marks(.f), ...)
    return(.f)
  }
  stopifnot(is.factor(.f))
  map <- list(...)
  n <- length(map)
  if(n == 0) return(.f)
  # mapping for 'other'
  if(any(isnul <- (lengths(map) == 0))) {
    if(sum(isnul) > 1)
      stop("At most one argument should be NULL or character(0)")
    otherlevels <- setdiff(levels(.f), unlist(map))
    map[[which(isnul)]] <- otherlevels
  }
  newlevels <- names(map)
  oldlevels <- levels(.f)
  mappedlevels <- unlist(map)
  if(sum(nzchar(newlevels)) != n)
    stop("Arguments must be in the form name=value")
  if(!all(mappedlevels %in% oldlevels))
    stop("Argument values must be levels of .f")
  ## construct mapping
  fullmap <- oldlevels
  for(i in seq_len(n)) {
    relevant <- oldlevels %in% map[[i]]
    fullmap[relevant] <- newlevels[i]
  }
  ## apply mapping
  newf <- factor(fullmap[.f], levels=unique(fullmap))
  return(newf)
}

levelsAsFactor <- function(x) {
  lev <- levels(x)
  if(is.null(lev)) return(NULL)
  return(factor(lev, levels=lev))
}

harmoniseLevels <- function(...) {
  x <- list(...)
  if(length(x) == 1) {
    x <- x[[1L]]
    if(!is.null(levels(x))) return(x) ## single factor or object
  }
  ## extract factor levels for each factor
  levlist <- lapply(x, levels)
  if(any(sapply(levlist, is.null)))
    stop("Some of the arguments do not have factor levels")
  if(length(unique(levlist)) == 1)
    return(x) # levels are already identical
  ## pool factor levels of all factors
  pooledlevels <- unique(unlist(levlist))
  matchlist <- lapply(levlist, match, table=pooledlevels)
  if(anyNA(unlist(matchlist)))
    stop("Unable to harmonise levels")
  ## map each factor to the pooled levels
  xentries <- lapply(x, "[")
  oldcodelist <- lapply(xentries, as.integer)
  newcodelist <- mapply("[", matchlist, oldcodelist, SIMPLIFY=FALSE)
  newfactors <- lapply(newcodelist, factor,
                       levels=seq_along(pooledlevels),
                       labels=pooledlevels)
  isim <- sapply(x, is.im)
  xnew <- newfactors
  xnew[isim] <- mapply("[<-", x=x[isim], value=newfactors, SIMPLIFY=FALSE)
  names(xnew) <- names(x)
  if(is.solist(x) || all(isim))
    xnew <- as.solist(xnew)
  return(xnew)
}
