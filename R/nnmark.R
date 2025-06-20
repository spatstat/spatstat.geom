#
# nnmark.R
#
# $Revision: 1.11 $ $Date: 2025/06/20 05:26:24 $

nnmark <- local({

  nnmark <- function(X, ..., k=1, at=c("pixels", "points"),
                     ties=c("first", "mean", "min", "max"),
                     proper=FALSE) {
    stopifnot(is.ppp(X))
    stopifnot(is.marked(X))
    at <- match.arg(at)
    if(!missing(ties) && is.function(ties)) {
      #' undocumented option: 'ties' is a function
      tiesfun <- ties
      ties <- "functie"
    } else {
      #' usual case
      ties <- match.arg(ties)
    }
    mX <- marks(X)
    if(ties != "first" && anyDuplicated(P <- unmark(X))) {
      ## pool marks of coincident points
      um <- uniquemap(P)
      uniq <- (um == seq_along(um))
      poolid <- cumsum(uniq)
      NX <- coerce.marks.numeric(X)
      mX <- marks(NX)
      smX <- split(as.data.frame(mX), factor(um))
      pmX <- switch(ties,
                    mean = { lapply(smX, colMeans) },
                    max =  {
                      lapply(smX, function(z) { apply(z, 2, max) })
                    },
                    min =  {
                      lapply(smX, function(z) { apply(z, 2, min) })
                    },
                    functie = {
                      lapply(smX, function(z, f=tiesfun) { apply(z, 2, f) })
                    })
      pmX <- do.call(rbind, unname(pmX))
      ## reassign to original pattern
      mX <- marksubset(pmX, poolid[um])
    }
    switch(at,
           pixels = {
             Y <- nnmap(X, k=k, what="which", ...)
             switch(markformat(mX),
                    vector={
                      result <- eval.im(mX[Y])
                    },
                    dataframe = {
                      mX <- as.list(as.data.frame(mX))
                      result <- solapply(mX, lookedup, indeximage=Y)
                      if(length(result) == 1)
                        result <- result[[1]]
                    },
                    stop("Marks must be a vector or dataframe"))
           },
           points = {
             if(!proper) {
               Y <- nnwhich(X, k=k)
             } else {
               ## find distinct points
               P <- unmark(X)
               um <- uniquemap(P)
               uid <- which(um == seq_along(um))
               U <- P[uid]
               ## find neighbours 
               Z <- nnwhich(U, k=k)
               ## map back
               Y <- uid[Z[um]]
             }
             switch(markformat(X),
                    vector={
                      result <- mX[Y]
                    },
                    dataframe = {
                      result <- mX[Y,, drop=FALSE]
                      row.names(result) <- NULL
                    },
                    stop("Marks must be a vector or dataframe"))
           })
    return(result)
  }

  lookedup <- function(xvals, indeximage) eval.im(xvals[indeximage])

  nnmark
})




  
