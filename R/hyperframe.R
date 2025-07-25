#
#  hyperframe.R
#
# $Revision: 1.83 $  $Date: 2025/07/26 01:54:10 $
#

## ------------------  utilities -------------------------

## recognise data that would be acceptable to data.frame

can.be.dfcolumn <- function(x) {
  is.atomic(x) &&
    (is.vector(x) ||
     is.factor(x) ||
     inherits(x, c("POSIXlt", "POSIXct", "Date", "Surv")))
}

## recognise data that would require a hypercolumn

can.be.hypercolumn <- function(x) {
  if(!is.list(x))
    return(FALSE)
  if(is.NAobject(x))
    return(FALSE)
  if(inherits(x, c("listof", "anylist")))
    return(TRUE)
  if(length(x) <= 1)
    return(TRUE)
  cla <- sapply(x, classIgnoringNA, first=TRUE)
  return(length(unique(cla)) == 1)
}


## ---------------  main definition -------------------------

hyperframe <- function(...,
                       row.names=NULL, check.rows=FALSE, check.names=TRUE,
                       stringsAsFactors=NULL) {
  aarg <- list(...)
  nama <- names(aarg)

  stringsAsFactors <- resolve.stringsAsFactors(stringsAsFactors)
    
  ## number of columns (= variables)
  nvars <- length(aarg)
  
  if(nvars == 0) {
    ## zero columns - return
    result <- list(nvars=0,
                   ncases=0,
                   vname=character(0),
                   vtype=factor(,
                                levels=c("dfcolumn","hypercolumn","hyperatom")),
                   vclass=character(0),
                   df=data.frame(),
                   hyperatoms=list(),
                   hypercolumns=list())
    class(result) <- c("hyperframe", class(result))
    return(result)
  }

  ## check column names
  if(is.null(nama)) {
    nama <- paste("V", 1:nvars, sep="")
  } else if(any(unnamed <- (nama == ""))) {
    nama[unnamed] <- paste("V", seq_len(sum(unnamed)), sep="")
  }
  nama <- make.names(nama, unique=TRUE)
  names(aarg) <- nama
  
  ## Each argument must be either
  ##    - a vector suitable as a column in a data frame
  ##    - a list of objects of the same class
  ##    - a single object of some class
  
  dfcolumns    <- sapply(aarg, can.be.dfcolumn)
  hypercolumns <- sapply(aarg, can.be.hypercolumn)
  hyperatoms   <- !(dfcolumns | hypercolumns)
  
  ## Determine number of rows (= cases) 
  columns <- dfcolumns | hypercolumns
  if(!any(columns)) {
    ncases <- 1
  } else {
    heights <- rep.int(1, nvars)
    heights[columns] <-  lengths(aarg[columns])
    u <- unique(heights)
    if(length(u) > 1) {
      u <- u[u != 1]
      if(length(u) > 1)
        stop(paste("Column lengths are inconsistent:",
                   paste(u, collapse=",")))
    }
    ncases <- u
    if(ncases > 1 && all(heights[dfcolumns] == 1)) {
      ## force the data frame to have 'ncases' rows
      aarg[dfcolumns] <- lapply(aarg[dfcolumns], rep, ncases)
      heights[dfcolumns] <- ncases
    }
    if(any(stubs <- hypercolumns & (heights != ncases))) {
      ## hypercolumns of height 1 should be hyperatoms
      aarg[stubs] <- lapply(aarg[stubs], "[[", i=1L)
      hypercolumns[stubs] <- FALSE
      hyperatoms[stubs] <- TRUE
    }
  }
  
  ## Collect the data frame columns into a data frame
  if(!any(dfcolumns)) {
    df <- as.data.frame(matrix(, ncases, 0))
  } else {
    df <- do.call(data.frame,
                  append(aarg[dfcolumns],
                         list(check.rows=check.rows,
                              check.names=check.names,
                              stringsAsFactors=stringsAsFactors)))
    names(df) <- nama[dfcolumns]
  }
  if(length(row.names)) row.names(df) <- row.names

  ## Storage type of each variable
  vtype <- character(nvars)
  vtype[dfcolumns] <- "dfcolumn"
  vtype[hypercolumns] <- "hypercolumn"
  vtype[hyperatoms] <- "hyperatom"
  vtype=factor(vtype, levels=c("dfcolumn","hypercolumn","hyperatom"))
  
  ## Class of each variable
  vclass <- character(nvars)
  if(any(dfcolumns))
    vclass[dfcolumns] <- sapply(as.list(df), classIgnoringNA, first=TRUE)
  if(any(hyperatoms))
    vclass[hyperatoms] <- sapply(aarg[hyperatoms], classIgnoringNA, first=TRUE)
  if(any(hypercolumns))
    vclass[hypercolumns] <- sapply(lapply(aarg[hypercolumns], "[[", i=1L),
                                   classIgnoringNA, first=TRUE)
  ## Put the result together
  result <- list(nvars=nvars,
                 ncases=ncases,
                 vname=nama,
                 vtype=vtype,
                 vclass=vclass,
                 df=df,
                 hyperatoms=aarg[hyperatoms],
                 hypercolumns=aarg[hypercolumns])
    
  class(result) <- c("hyperframe", class(result))
  return(result)
}


## ...................   methods -------------------------

is.hyperframe <- function(x) inherits(x, "hyperframe")

print.hyperframe <- function(x, ...) {
  ux <- unclass(x)
  nvars <- ux$nvars
  ncases <- ux$ncases
  if(nvars * ncases == 0) {
    splat("NULL hyperframe with", ncases,
          ngettext(ncases, "row (=case)", "rows (=cases)"),
          "and", nvars,
          ngettext(nvars, "column (=variable)", "columns (=variables)"))
  } else {
    if(waxlyrical('gory')) cat("Hyperframe:\n")
    print(as.data.frame(x, discard=FALSE), ...)
  }
  return(invisible(NULL))
}

dim.hyperframe <- function(x) {
  with(unclass(x), c(ncases, nvars))
}

summary.hyperframe <- function(object, ..., brief=FALSE) {
  x <- unclass(object)
  y <- list(
            nvars = x$nvars,
            ncases = x$ncases,
            dim = c(x$ncases, x$nvars),
            typeframe = data.frame(VariableName=x$vname, Class=x$vclass),
            storage = x$vtype,
            col.names = x$vname)
  classes <- x$vclass
  names(classes) <- x$vname
  y$classes <- classes
  # Ordinary data frame columns
  df <- x$df
  y$dfnames <- colnames(df)
  y$df <- if(length(df) > 0 && !brief) summary(df) else NULL
  y$row.names <- row.names(df)
  # insert into full array
  if(!brief && x$nvars > 0) {
    isobject <- (x$vtype != "dfcolumn")
    nobj <- sum(isobject)
    if(nobj == 0) {
      allcols <- y$df
    } else {
      nas <- rep(list(NA_character_), nobj)
      names(nas) <- x$vname[isobject]
      allcols <- do.call(cbind, append(list(y$df), nas))
      acnames <- c(colnames(df), names(nas))
      allcols <- allcols[ , match(x$vname, acnames), drop=FALSE]
    }
    pclass <- padtowidth(paren(classes), colnames(allcols), justify="right")
    allcols <- as.table(rbind(class=pclass, as.table(allcols)))
    row.names(allcols) <- rep("", nrow(allcols))
    y$allcols <- allcols
  }
  class(y) <- c("summary.hyperframe", class(y))
  return(y)
}

print.summary.hyperframe <- function(x, ...) {
  nvars <- x$nvars
  ncases <- x$ncases
  splat(if(nvars * ncases == 0) "NULL hyperframe" else "hyperframe",
        "with", ncases,
        ngettext(ncases, "row", "rows"),
        "and", nvars,
        ngettext(nvars, "column", "columns"))
  if(nvars == 0)
    return(invisible(NULL))
  print(if(any(x$storage == "dfcolumn")) x$allcols else noquote(x$classes))
  return(invisible(NULL))
}

names.hyperframe <- function(x) { unclass(x)$vname }

"names<-.hyperframe" <- function(x, value) {
  x <- unclass(x)
  stopifnot(is.character(value))
  value <- make.names(value)
  if(length(value) != x$nvars)
    stop("Incorrect length for vector of names")
  vtype <- x$vtype
  names(x$df)           <- value[vtype == "dfcolumn"]
  names(x$hyperatoms)   <- value[vtype == "hyperatom"]
  names(x$hypercolumns) <- value[vtype == "hypercolumn"]
  x$vname <- value
  class(x) <- c("hyperframe", class(x))
  return(x)
}

row.names.hyperframe <- function(x) {
  return(row.names(unclass(x)$df))
}

"row.names<-.hyperframe" <- function(x, value) {
  y <- unclass(x)
  row.names(y$df) <- value
  class(y) <- c("hyperframe", class(y))
  return(y)
}

dimnames.hyperframe <- function(x) {
  ux <- unclass(x)
  return(list(row.names(ux$df), ux$vname))
}

"dimnames<-.hyperframe" <- function(x, value) {
  if(!is.list(value) || length(value) != 2 || !all(sapply(value, is.character)))
    stop("Invalid 'dimnames' for a hyperframe", call.=FALSE)
  rn <- value[[1L]]
  cn <- value[[2L]]
  d <- dim(x)
  if(length(rn) != d[1L])
    stop(paste("Row names have wrong length:",
               length(rn), "should be", d[1L]),
         call.=FALSE)
  if(length(cn) != d[2L])
    stop(paste("Column names have wrong length:",
               length(cn), "should be", d[2L]),
         call.=FALSE)
  y <- unclass(x)
  row.names(y$df) <- value[[1L]]
  y$vname <- value[[2]]
  class(y) <- c("hyperframe", class(y))
  return(y)
}

## conversion to hyperframe

as.hyperframe <- function(x, ...) {
  UseMethod("as.hyperframe")
}

as.hyperframe.hyperframe <- function(x, ...) {
  return(x)
}

as.hyperframe.data.frame <- function(x, ..., stringsAsFactors=FALSE) {
  if(missing(x) || is.null(x)) {
    xlist <- rona <- NULL
  } else {
    rona <- row.names(x)
    xlist <- as.list(x)
  }
  do.call(hyperframe,
          resolve.defaults(xlist,
                           list(...),
                           list(row.names=rona,
                                stringsAsFactors=stringsAsFactors),
                           .StripNull=TRUE))
}

as.hyperframe.anylist <- 
as.hyperframe.listof <- function(x, ...) {
  if(!missing(x)) {
    xname <- sensiblevarname(short.deparse(substitute(x)), "x")
    xlist <- list(x)
    names(xlist) <- xname
  } else xlist <- NULL
  do.call(hyperframe,
          resolve.defaults(
                           xlist,
                           list(...),
                           .StripNull=TRUE))
}

as.hyperframe.default <- function(x, ...) {
  as.hyperframe(as.data.frame(x, ...))
}

#### conversion to other types

as.data.frame.hyperframe <- function(x, row.names = NULL,
                                     optional = FALSE, ...,
                                     discard=TRUE, warn=TRUE) {
  ux <- unclass(x)
  if(is.null(row.names))
    row.names <- row.names(ux$df)
  vtype <- ux$vtype
  vclass <- ux$vclass
  ishyper <- (vtype != "dfcolumn")
  nhyper <- sum(ishyper)
  if(discard) { 
    if(nhyper > 0 && warn)
      warning(paste(nhyper, 
                    ngettext(nhyper, "variable", "variables"),
                    "discarded in conversion to data frame"))
    df <- as.data.frame(ux$df, row.names=row.names, optional=optional, ...)
  } else {
    ## convert to a list to be passed to 'data.frame'
    lx <- as.list(x)
    ## handle non-data-frame columns
    if(nhyper > 0) {
      ## objects of class 'foo' are represented by "(foo)"
      vclassstring <- paren(vclass)
      nrows <- ux$ncases
      lx[ishyper] <- lapply(as.list(vclassstring[ishyper]),
                            rep.int, times=nrows)
    }
    df <- do.call(data.frame, append(lx, list(row.names=row.names)))
    colnames(df) <- ux$vname
    ## detect NA objects using is.na.hyperframe
    if(nhyper > 0 && any(isna <- is.na(x))) {
      if(TRUE) {
        ## print <NA> for each missing object
        df[isna] <- NA
      } else {
        ## print <NA foo> for missing object of class foo
        J <- col(df)
        isnaobj <- isna & (ishyper[J])
        if(any(isnaobj))
          df[isnaobj] <- paste0("<NA_", vclass[J[isnaobj]], ">")
      }
    }
  }
  return(df)
}

is.na.hyperframe <- function(x) {
  y <- matrix(FALSE, nrow(x), ncol(x), dimnames=dimnames(x))
  ux <- unclass(x)
  vname <- ux$vname
  vtype <- ux$vtype
  vtype <- levels(vtype)[as.integer(vtype)]
  hyperatoms <- ux$hyperatoms
  hypercolumns <- ux$hypercolumns
  for(j in seq_len(ncol(x))) {
    switch(vtype[j],
           dfcolumn = {
             y[,j] <- is.na(x[,j,drop=TRUE])
           },
           hyperatom = {
             vn <- vname[[j]]
             if(is.NAobject(hyperatoms[[vn]]))
               y[,j] <- TRUE
           },
           hypercolumn = {
             vn <- vname[[j]]
             y[,j] <- sapply(hypercolumns[[vn]], is.NAobject)
           })
  }
  return(y)
}
             
as.list.hyperframe <- function(x, ..., expandatoms=TRUE) {
  ux <- unclass(x)
  out <- vector(mode="list", length=ux$nvars)
  vtype <- ux$vtype
  df <- ux$df
  if(any(dfcol <- (vtype == "dfcolumn")))
    out[dfcol] <- as.list(df)
  if(any(hypcol <- (vtype == "hypercolumn"))) {
    hc <- lapply(ux$hypercolumns, as.solist, demote=TRUE)
    out[hypcol] <- hc
  }
  if(any(hatom <- (vtype == "hyperatom"))) {
    ha <- ux$hyperatoms
    names(ha) <- NULL
    if(expandatoms) {
      hacol <- lapply(ha, list)
      hacol <- lapply(hacol, rep.int, times=ux$ncases)
      hacol <- lapply(hacol, as.solist, demote=TRUE)
      out[hatom] <- hacol
    } else {
      out[hatom] <- ha
    }
  }
  fullrows <- !hatom | expandatoms
  out[fullrows] <- lapply(out[fullrows], "names<-", value=row.names(df))
  names(out) <- names(x)
  return(out)
}

# evaluation

# eval.hyper <- function(e, h, simplify=TRUE, ee=NULL) {
#   .Deprecated("with.hyperframe", package="spatstat")
#   if(is.null(ee))
#     ee <- as.expression(substitute(e))
#   with.hyperframe(h, simplify=simplify, ee=ee)
# }

with.hyperframe <- function(data, expr, ..., simplify=TRUE, ee=NULL,
                            enclos=NULL) {
  if(!inherits(data, "hyperframe"))
    stop("data must be a hyperframe")
  if(is.null(ee))
    ee <- as.expression(substitute(expr))
  if(is.null(enclos))
    enclos <- parent.frame()
  n <- nrow(data)
  out <- vector(mode="list", length=n)
  #' consider all variables or functions provided in 'data'
  nama <- intersect(all.names(ee), colnames(data))
  if(length(nama)) {
    #' check for NA value or NAobject among these variables or functions
    bad <- apply(is.na(data)[, nama, drop=FALSE], 1, any)
    goodrows <- which(!bad)
    out[bad] <- NA
  } else {
    goodrows <- seq_len(n)
  }
  datalist <- as.list(data)
  for(i in goodrows) {
    rowi <- lapply(datalist, "[[", i=i)  # ensures the result is always a list
    outi <- eval(ee, rowi, enclos)
    if(!is.null(outi))
      out[[i]] <- outi
  }
  names(out) <- row.names(data)
  if(simplify && all(unlist(lapply(out, is.vector)))) {
    # if all results are atomic vectors of equal length,
    # return a matrix or vector.
    lenfs <- lengths(out)
    if(all(unlist(lapply(out, is.atomic))) &&
            length(unique(lenfs)) == 1) {
      out <- t(as.matrix(as.data.frame(out)))
      row.names(out) <- row.names(data)
      out <- out[,,drop=TRUE]
      return(out)
    }
  }
  ## apply rules for a column in a hyperframe 
  out <- hyperframe(result=as.anylist(out), row.names=row.names(data))$result
  return(out)
}

cbind.hyperframe <- function(...) {
  aarg <- list(...)
  narg <- length(aarg)
  if(narg == 0)
    return(hyperframe())
  namarg <- names(aarg)
  if(is.null(namarg))
    namarg <- rep.int("", narg)
  ishyper <- unlist(lapply(aarg, inherits, what="hyperframe"))
  isdf <- unlist(lapply(aarg, inherits, what="data.frame"))
  columns <- list()
  for(i in 1:narg) {
    if(ishyper[i] || isdf[i]){
      if(ncol(aarg[[i]]) > 0) {
        newcolumns <- as.list(aarg[[i]])
        if(namarg[i] != "")
          names(newcolumns) <- paste(namarg[i], ".", names(newcolumns), sep="")
        columns <- append(columns, newcolumns)
      }
    } else {
      nextcolumn <- list(aarg[[i]])
      names(nextcolumn) <- namarg[i]
      columns <- append(columns, nextcolumn)
    }
  }
  result <- do.call(hyperframe, columns)
  ## tack on row names
  rona <- lapply(aarg, row.names)
  good <- (lengths(rona) == nrow(result))
  if(any(good)) 
    row.names(result) <- rona[[min(which(good))]]
  return(result)
}

rbind.hyperframe <- function(...) {
  argh <- list(...)
  if(length(argh) == 0)
    return(NULL)
  # convert them all to hyperframes
  argh <- lapply(argh, as.hyperframe)
  #
  nargh <- length(argh)
  if(nargh == 1)
    return(argh[[1L]])
  # check for compatibility of dimensions & names
  dfs <- lapply(argh, as.data.frame, discard=FALSE)
  dfall <- do.call(rbind, dfs)
  # check that data frame columns also match
  dfs0 <- lapply(argh, as.data.frame, discard=TRUE, warn=FALSE)
  df0all <- do.call(rbind, dfs0)
  # assemble data
  rslt <- list()
  nam <- names(dfall) 
  nam0 <- names(df0all)
  for(k in seq_along(nam)) {
    nama <- nam[k]
    if(nama %in% nam0) {
      # data frame column: already made
      rslt[[k]] <- dfall[,k]
    } else {
      ## hypercolumns or hyperatoms: extract them
      hdata <- lapply(argh, "[", j=nama, drop=FALSE)
      hdata <- lapply(lapply(hdata, as.list), getElement, name=nama)
      ## bind them
      hh <- Reduce(append, hdata)
      rslt[[k]] <- hh
    }
  }
  ## make hyperframe
  names(rslt) <- nam
  rona <- row.names(dfall)
  out <- do.call(hyperframe, append(rslt,
                                    list(stringsAsFactors=FALSE,
                                         row.names=rona)))
  return(out)
}

plot.hyperframe <-
  function(x, e, ..., main, arrange=TRUE,
           nrows=NULL, ncols=NULL,
           parargs=list(mar=mar * marsize),
           marsize=1, mar=c(1,1,3,1)) {
  xname <- short.deparse(substitute(x))
  main <- if(!missing(main)) main else xname
  mar <- rep(mar, 4)[1:4]
  
  if(missing(e)) {
    # default: plot first column that contains objects
    ok <- (summary(x)$storage %in% c("hypercolumn", "hyperatom"))
    if(any(ok)) {
      j <- min(which(ok))
      x <- x[,j, drop=TRUE, strip=FALSE]
      x <- as.solist(x, demote=TRUE)
      plot(x, ..., main=main, arrange=arrange, nrows=nrows, ncols=ncols)
      return(invisible(NULL))
    } else {
      # hyperframe does not contain any objects
      # invoke plot.data.frame
      x <- as.data.frame(x)
      plot(x, ..., main=main)
      return(invisible(NULL))
    }
  }

  if(!is.language(e))
    stop(paste("Argument e should be a call or an expression;",
               "use quote(...) or expression(...)"))
  ee <- as.expression(e)

  if(!arrange) {
    # No arrangement specified: just evaluate the plot expression 'nr' times
    with(x, ee=ee)
    return(invisible(NULL))
  }

  # Arrangement
  # Decide whether to plot a main header
  banner <- (sum(nchar(as.character(main))) > 0)
  if(length(main) > 1)
    main <- paste(main, collapse="\n")
  nlines <- if(!is.character(main)) 1 else length(unlist(strsplit(main, "\n")))
  # determine arrangement of plots
  # arrange like mfrow(nrows, ncols) plus a banner at the top
  n <- summary(x)$ncases
  if(is.null(nrows) && is.null(ncols)) {
    nrows <- as.integer(floor(sqrt(n)))
    ncols <- as.integer(ceiling(n/nrows))
  } else if(!is.null(nrows) && is.null(ncols))
    ncols <- as.integer(ceiling(n/nrows))
  else if(is.null(nrows) && !is.null(ncols))
    nrows <- as.integer(ceiling(n/ncols))
  else stopifnot(nrows * ncols >= length(x))
  nblank <- ncols * nrows - n
  # declare layout
  mat <- matrix(c(seq_len(n), numeric(nblank)),
                byrow=TRUE, ncol=ncols, nrow=nrows)
  heights <- rep.int(1, nrows)
  if(banner) {
    # Increment existing panel numbers
    # New panel 1 is the banner
    panels <- (mat > 0)
    mat[panels] <- mat[panels] + 1L
    mat <- rbind(rep.int(1,ncols), mat)
    heights <- c(0.1 * (1 + nlines), heights)
  }
  # initialise plot
  layout(mat, heights=heights)
  # plot banner
  if(banner) {
    opa <- par(mar=rep.int(0,4), xpd=TRUE)
    on.exit(par(opa))
    plot(numeric(0),numeric(0),type="n",ann=FALSE,axes=FALSE,
         xlim=c(-1,1),ylim=c(-1,1))
    cex <- resolve.defaults(list(...), list(cex.title=2))$cex.title
    text(0,0,main, cex=cex)
  }
  # plot panels
  npa <- do.call(par, parargs)
  if(!banner) on.exit(par(npa))
  with(x, ee=ee)
  # revert
  layout(1)
  return(invisible(NULL))
}


str.hyperframe <- function(object, ...) {
  d <- dim(object)
  x <- unclass(object)
  argh <- resolve.defaults(list(...), list(nest.lev=0, indent.str="  .."))
  cat(paste("'hyperframe':\t",
            d[1L], ngettext(d[1L], "row", "rows"),
            "and",
            d[2L], ngettext(d[2L], "column", "columns"),
            "\n"))
  nr <- d[1L]
  nc <- d[2L]
  if(nc > 0) {
    vname <- x$vname
    vclass <- x$vclass
    vtype  <- as.character(x$vtype)
    indentstring <- with(argh, paste(rep.int(indent.str, nest.lev), collapse=""))
    for(j in 1:nc) {
      tag <- paste("$", vname[j])
      switch(vtype[j],
             dfcolumn={
               desc <- vclass[j]
               if(nr > 0) {
                 vals <- object[1:min(nr,3),j,drop=TRUE]
                 vals <- paste(paste(format(vals), collapse=" "), "...")
               } else vals <- ""
             },
             hypercolumn=,
             hyperatom={
               desc <- "objects of class"
               vals <- vclass[j]
             })
      cat(paste(paste(indentstring, tag, sep=""),
                ":", desc, vals, "\n"))
    }
  }
  return(invisible(NULL))
}

subset.hyperframe <- function(x, subset, select, ...) {
  stopifnot(is.hyperframe(x))
  r <- if(missing(subset)) {
    rep_len(TRUE, nrow(x))
  } else {
      r <- eval(substitute(
        with(x, e, enclos=parent.frame()),
        list(e=substitute(subset))))
    if (!is.logical(r)) 
      stop("'subset' must be logical")
    r & !is.na(r)
  }
  vars <- if(missing(select)) { 
    TRUE
  } else {
    nl <- as.list(seq_len(ncol(x)))
    names(nl) <- names(x)
    eval(substitute(select), nl, parent.frame())
  }
  nama <- names(x)
  names(nama) <- nama
  vars <- nama[vars]
  z <- x[i=r, j=vars, ...]
  return(z)
}

head.hyperframe <- function (x, n = 6L, ...) {
  stopifnot(length(n) == 1L)
  n <- if(n < 0L) max(nrow(x) + n, 0L) else min(n, nrow(x))
  x[seq_len(n), , drop = FALSE]
}

tail.hyperframe <- function(x, n = 6L, ...) {
  stopifnot(length(n) == 1L)
  nrx <- nrow(x)
  n <- if(n < 0L) max(nrx + n, 0L) else min(n, nrx)
  sel <- seq.int(to = nrx, length.out = n)
  x[sel, , drop = FALSE]
}

edit.hyperframe <- function(name, ...) {
  x <- name
  isdf <- unclass(x)$vtype == "dfcolumn"
  if(!any(isdf)) {
    warning("No columns of editable data", call.=FALSE)
    return(x)
  }
  y <- x[,isdf]
  ynew <- edit(as.data.frame(y), ...)
  xnew <- x
  for(na in names(ynew)) xnew[,na] <- ynew[,na]
  losenames <- setdiff(names(y), names(ynew))
  for(na in losenames) xnew[,na] <- NULL
  return(xnew)
}
