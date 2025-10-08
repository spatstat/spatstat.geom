##
## hypersub.R
##
##
##  subset operations for hyperframes
##
##  $Revision: 1.37 $    $Date: 2025/07/07 06:14:38 $
##

"[.hyperframe" <- function(x, i, j, drop=FALSE, strip=drop, ...) {
  x <- unclass(x)
  if(!missing(i)) {
    if(length(dim(i)) > 1)
      stop("Matrix index i is not supported in '[.hyperframe'", call.=FALSE)
    y <- x
    y$df     <- x$df[i, , drop=FALSE]
    y$ncases <- nrow(y$df)
    y$hypercolumns <- lapply(x$hypercolumns, "[", i=i)
    x <- y
  }
  if(!missing(j)) {
    if(length(dim(j)) > 1)
      stop("Matrix index j is not supported in '[.hyperframe'", call.=FALSE)
    y <- x
    patsy <- seq_len(y$nvars)
    names(patsy) <- y$vname
    jj <- patsy[j]
    names(jj) <- NULL
    y$nvars <- length(jj)
    y$vname <- vname <- x$vname[jj]
    y$vtype <- vtype <- x$vtype[jj]
    y$vclass <- x$vclass[jj]
    if(ncol(x$df) != 0) 
      y$df    <- x$df[ , vname[vtype == "dfcolumn"], drop=FALSE]
    y$hyperatoms <- x$hyperatoms[ vname[ vtype == "hyperatom" ]]
    y$hypercolumns <- x$hypercolumns[ vname [ vtype == "hypercolumn" ] ]
    x <- y
  }
  if(drop) {
    nrows <- x$ncases
    ncols <- x$nvars
    if(nrows == 1 && ncols == 1 && strip) {
      ## return a single object 
      y <- switch(as.character(x$vtype),
                  dfcolumn    = x$df[, , drop=TRUE],
                  hypercolumn = (x$hypercolumns[[1L]])[[1L]],
                  hyperatom   = x$hyperatoms[[1L]])
      return(y)
    } else if(nrows == 1) {
      ## return the row as a vector or a list
      if(strip && all(x$vtype == "dfcolumn"))
        return(x$df[ , , drop=TRUE])
      n <- x$nvars
      y <- vector(mode="list", length=n)
      names(y) <- nama <- x$vname
      for(k in seq_len(n)) {
        namk <- nama[k]
        y[[k]] <- switch(as.character(x$vtype[k]),
                         dfcolumn = x$df[ , namk, drop=TRUE],
                         hyperatom = x$hyperatoms[[namk]],
                         hypercolumn = (x$hypercolumns[[namk]])[[1L]]
                         )
      }
      return(as.solist(y, demote=TRUE))
    } else if(ncols == 1) {
      ## return a column as an 'anylist'/'solist' or a vector
      switch(as.character(x$vtype),
             dfcolumn = {
               return(x$df[, , drop=TRUE])
             },
             hypercolumn = {
               y <- as.solist(x$hypercolumns[[1L]], demote=TRUE)
               names(y) <- row.names(x$df)
               return(y)
             },
             hyperatom = {
               ## replicate it to make a hypercolumn
               ha <- x$hyperatoms[1L]
               names(ha) <- NULL
               hc <- rep.int(ha, x$ncases)
               hc <- as.solist(hc, demote=TRUE)
               names(hc) <- row.names(x$df)
               return(hc)
             }
           )
    }
  }
  class(x) <- unique(c("hyperframe", class(x)))
  return(x)
}

"$.hyperframe" <- function(x,name) {
  m <- match(name, unclass(x)$vname)
  if(is.na(m))
    return(NULL)
  return(x[, name, drop=TRUE, strip=FALSE])
}

"$<-.hyperframe" <- function(x, name, value) {
  y <- as.list(x, expandatoms=FALSE)
  if(is.hyperframe(value)) {
    if(ncol(value) == 1) {
      y[name] <- as.list(value)
    } else {
      y <- insertinlist(y, name, as.list(value))
    }
  } else {
    if(!is.null(value) && !can.be.dfcolumn(value)) {
      value <- as.list(value)
      if(inherits(value, "anylist")) {
        ## all entries must have common class
        ## catch vanilla NA entries and coerce to appropriate type
        value <- coerceNAtoObject(value)
        ## check for conflicting types
        cls <- unique(sapply(value, classIgnoringNA, first=TRUE))
        if(length(cls) > 1)
          stop(paste("Column entries have conflicting types",
                     commasep(sQuote(cls))),
               call.=FALSE)
      }
    }
    y[[name]] <- value
  }
  z <- do.call(hyperframe, append(y, list(row.names=row.names(x),
                                            stringsAsFactors=FALSE)))
  return(z)
}



"[<-.hyperframe" <- 
function (x, i, j, value)
{
  sumry <- summary(x)
  colnam <- sumry$col.names
  dimx <- sumry$dim
  coltype <- sumry$storage
  colclass <- sumry$classes
  igiven <- !missing(i)
  jgiven <- !missing(j)
  if(igiven) {
    if(length(dim(i)) > 1)
      stop("Matrix index i is not supported in '[<-.hyperframe'", call.=FALSE)
    singlerow    <- ((check.1.integer(i, fatal=FALSE, warn=FALSE) && i > 0)
                    || (is.character(i) && length(i) == 1)
                    || (is.logical(i) && sum(i) == 1))
  } else {
    i <- seq_len(dimx[1L])
    singlerow <- FALSE
  }
  if(jgiven) {
    if(length(dim(j)) > 1)
      stop("Matrix index j is not supported in '[<-.hyperframe'", call.=FALSE)
    singlecolumn <- ((check.1.integer(j, fatal=FALSE, warn=FALSE) && j > 0)
                    || (is.character(j) && length(j) == 1)
                    || (is.logical(j) && sum(j) == 1))
  } else {
    j <- seq_len(dimx[2L])
    singlecolumn <- FALSE
  }
    
  if(!igiven && jgiven) {
    # x[, j] <- value
    if(singlecolumn) {
      # expecting single hypercolumn
      if(is.logical(j)) j <- names(x)[j]
      y <- get("$<-.hyperframe")(x, j, value)
    } else {
      # expecting hyperframe 
      xlist <- as.list(x)
      xlist[j] <- as.list(as.hyperframe(value))
      # the above construction accepts all indices including extra entries
      y <- do.call(hyperframe, append(xlist,
                                        list(row.names=row.names(x))))
    }
  } else {
    ## x[, ] <- value or x[i, ] <- value or x[i,j] <- value 
    ## convert indices to positive integers
    rowseq <- seq_len(dimx[1L])
    colseq <- seq_len(dimx[2L])
    names(rowseq) <- row.names(x)
    names(colseq) <- colnam
    I <- as.integer(unname(rowseq[i]))
    J <- as.integer(unname(colseq[j]))
    ## convert to lists (automatically replicates hyperatoms to make columns)
    xlist <- as.list(x)
    ## modify xlist
    if(singlerow && singlecolumn) {
      ## modify single entry x[I,J]
      if(coltype[J] != "dfcolumn") {
        classJ <- colclass[J]
        ## ensure replacement entry belongs to same class as other entries
        if(identical(value, NA)) {
          ## coerce 'NA' to an NA object of the appropriate class
          value <- NAobject(classJ)
        } else {
          ## check replacement value has correct class
          if(!inherits(value, classJ)) {
            ## perhaps a list of 1 item of the required class
            if(is.list(value) && length(value) == 1 &&
               inherits(value[[1L]], classJ)) {
              value <- value[[1L]]
            } else {
              stop(paste("Replacement value does not have class",
                         sQuote(classJ)),
                   call.=FALSE)
            }
          }
        }
      }
      xlist[[J]][[I]] <- value
    } else {
      ## modify multiple entries
      hv <- if(is.hyperframe(value)) value else
            if(is.atomic(value)) as.hyperframe(as.data.frame(value)) else
            as.hyperframe(as.solist(value, demote=TRUE))
      vlist <- as.list(hv)
      nrowV <- dim(hv)[1L]
      ncolV <- dim(hv)[2L]
      if(nrowV != length(I)) {
        if(nrowV == 1) {
          ## replicate
          vlist <- lapply(vlist, rep, times=length(I))
        } else stop(paste("Replacement value has wrong number of rows:",
                          nrowV, "should be", length(I)),
                    call.=FALSE)
      }
      if(ncolV != length(J)) {
        if(ncolV == 1) {
          ## replicate
          vlist <- rep(vlist, times=length(J))
        } else stop(paste("Replacement value has wrong number of columns:",
                          ncolV, "should be", length(J)),
                    call.=FALSE)
      }
      ## 
      fullcolumn <- identical(I, rowseq)
      ## replace entries
      for(k in seq_along(J)) {
        jj <- J[k]
        valuesjjI <- vlist[[k]]
        ## replace entries in subset I of column jj by entries of 'valuesjjI'   
        if(!can.be.dfcolumn(valuesjjI)) {
          ## replacement values are objects
          if(fullcolumn) {
            ## Replace entire column jj by another list of objects
            ## (possibly of a different class)
            ## Coerce vanilla NA entries to appropriate class, whatever that is
            xlist[[jj]] <- coerceNAtoObject(valuesjjI)
          } else {
            ## Replace entries in proper subset of column jj by objects
            if(coltype[jj] == "dfcolumn")
              stop("Cannot replace entries in an atomic vector by objects",
                   call.=FALSE)
            ## Replacing objects by objects
            classjj <- colclass[jj]
            ## catch vanilla NA entries and coerce to spatial objects
            valuesjjI <- coerceNAtoObject(valuesjjI, classjj)
            if(!all(sapply(valuesjjI, inherits, what=classjj)))
            stop(paste("Replacement value in column", jj,
                       "is not of the required class", classjj),
                 call.=FALSE)
          }
        }
        xlist[[jj]][I] <- valuesjjI
      }
    }
    ## re-compress any hyperatoms that were not affected
    isatom <- (coltype == "hyperatom")
    changed <- seq_len(dimx[2L]) %in% J
    if(any(ha <- isatom & !changed)) {
      for(k in which(ha))
        xlist[[k]] <- xlist[[k]][[1L]]
    }
    ## put back together
    y <- do.call(hyperframe, append(xlist,
                                      list(row.names=row.names(x))))
  }
  return(y)
}

    
"[[.hyperframe" <-
function(x, ...)
{
  rr <- as.data.frame(row(x))
  cc <- as.data.frame(col(x))
  dimnames(rr) <- dimnames(cc) <- dimnames(x)
  chosen.rows <- unique(as.integer(rr[[...]]))
  chosen.cols <- unique(as.integer(cc[[...]]))
  nr <- length(chosen.rows)
  nc <- length(chosen.cols)
  if(nc == 0 || nr == 0) {
    ## should never be reached
    stop("No data selected", call.=FALSE)
  } else if(nc > 1) {
    ## should never be reached
    stop("More than one item (or column) of data selected", call.=FALSE)
  }
  if(nr == 1) {
    ## single item
    result <- x[chosen.rows, chosen.cols, drop=TRUE, strip=TRUE]
  } else if(length(chosen.rows) == nrow(rr)) {
    ## column
    result <- x[,chosen.cols, drop=TRUE, strip=FALSE]
  } else {
    ## subset of a column
    stop("Cannot select part of a column in '[['", call.=FALSE)
  }
  return(result)
}
                            
"[[<-.hyperframe" <-
function(x, i, j, value)
{
  ## detect 'blank' arguments like second argument in x[i, ] 
  ngiven <- length(sys.call())
  nmatched <- length(match.call())
  nblank <- ngiven - nmatched
  itype <- if(missing(i)) "absent" else "given"
  jtype <- if(missing(j)) "absent" else "given"
  if(nblank == 1) {
    if(!missing(i)) jtype <- "blank"
    if(!missing(j)) itype <- "blank"
  } else if(nblank == 2) {
    itype <- jtype <- "blank"
  }
  ## detect idiom x[[ ]] or x[[ , ]]
  if(itype != "given" && jtype != "given" && prod(dim(x)) > 1)
    stop("More than one cell or column of cells selected", call.=FALSE)
  ## find selected rows and columns
  rr <- as.data.frame(row(x))
  cc <- as.data.frame(col(x))
  dimnames(rr) <- dimnames(cc) <- dimnames(x)
  switch(paste0(itype, jtype),
         givengiven = {
           chosen.rows <- rr[[i, j]]
           chosen.cols <- cc[[i, j]]
         },
         givenabsent = {
           chosen.rows <- rr[[i]]
           chosen.cols <- cc[[i]]
         },
         givenblank = {
           chosen.rows <- rr[[i, ]]
           chosen.cols <- cc[[i, ]]
         },
         absentgiven = {
           ## cannot occur
           chosen.rows <- rr[[, j]]
           chosen.cols <- cc[[, j]]
         },
         absentabsent = {
           chosen.rows <- rr[[ ]]
           chosen.cols <- cc[[ ]]
         },
         absentblank = {
           ## cannot occur
           chosen.rows <- rr[[ , ]]
           chosen.cols <- cc[[ , ]]
         },
         blankgiven = {
           chosen.rows <- rr[[, j]]
           chosen.cols <- cc[[, j]]
         },
         blankabsent = {
           ## cannot occur
           chosen.rows <- rr[[]]
           chosen.cols <- cc[[]]
         },
         blankblank = {
           chosen.rows <- rr[[ , ]]
           chosen.cols <- cc[[ , ]]
         })
  chosen.rows <- unique(as.integer(chosen.rows))
  chosen.cols <- unique(as.integer(chosen.cols))
  nr <- length(chosen.rows)
  nc <- length(chosen.cols)
  if(nc == 0 || nr == 0) {
    ## should never be reached
    stop("No cells selected", call.=FALSE)
  } else if(nc > 1) {
    ## should never be reached
    stop("More than one cell or column of cells selected", call.=FALSE)
  }
  if(nr == 1) {
    ## single item
    xj <- x[[chosen.cols]]
    if(!is.atomic(xj)) {
      ## ensure replacement value has same class
      vcj <- unclass(x)$vclass[[chosen.cols]]
      if(identical(value, NA)) {
        ## coerce to NA object of required class
        value <- NAobject(vcj)
      } else if(!inherits(value, vcj)) {
        stop(paste("Replacement value does not have required class",
                   sQuote(vcj)),
             call.=FALSE)
      }
    }
    xj[[chosen.rows]] <- value
    x[,chosen.cols] <- xj
  } else if(length(chosen.rows) == nrow(rr)) {
    ## column
    x[,chosen.cols] <- value
  } else {
    ## subset of a column
    stop("Cannot assign part of a column in '[[<-'", call.=FALSE)
  }
  return(x)
}

split.hyperframe <- local({

  split.hyperframe <- function(x, f, drop=FALSE, ...) {
    y <- data.frame(id=seq_len(nrow(x)))
    z <- split(y, f, drop=drop)
    z <- lapply(z, getElement, name="id")
    out <- lapply(z, indexi, x=x)
    return(out)
  }

  indexi <- function(i, x) x[i,]
  
  split.hyperframe
})


"split<-.hyperframe" <- function(x, f, drop=FALSE, ..., value) {
  ix <- split(seq_len(nrow(x)), f, drop = drop, ...)
  n <- length(value)
  j <- 0
  for (i in ix) {
    j <- j%%n + 1L
    x[i, ] <- value[[j]]
  }
  x
}
  
