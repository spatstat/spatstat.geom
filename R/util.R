##
##    util.R    miscellaneous utilities
##
##    $Revision: 1.269 $    $Date: 2026/02/24 00:50:07 $
##

## common invocation of matrixsample

rastersample <- function(X, Y=NULL, ..., phase=0L, scale=1L) {
  stopifnot(is.im(X) || is.mask(X))
  if(!is.null(Y)) {
    ## sample X onto raster of Y
    stopifnot(is.im(Y) || is.mask(Y))
    phase <- c((Y$yrow[1] - X$yrow[1])/X$ystep,
               (Y$xcol[1] - X$xcol[1])/X$xstep)
    scale <- c(Y$ystep/X$ystep,
               Y$xstep/X$xstep)
  } else if(missing(phase) && missing(scale)) {
    return(X)
  } else {
    ## create a coarsened version of raster of X
    ## phase and scale are expressed as multiples of the pixels of X
    if(!is.integer(phase)) {
      phase <- as.numeric(phase)
      if(all(phase %% 1 == 0))
        phase <- as.integer(phase)
    }
    if(!is.integer(scale)) {
      scale <- as.numeric(scale)
      if(all(scale %% 1 == 0))
        scale <- as.integer(scale)
    }
    phase <- ensure2vector(phase)
    scale <- ensure2vector(scale)
    if(is.integer(phase) && is.integer(scale)) {
      ## subsample the existing coordinates
      yrow <- X$yrow[seq(from=1L+phase[1L], by=scale[1L], to=length(X$yrow))]
      xcol <- X$xcol[seq(from=1L+phase[2L], by=scale[2L], to=length(X$xcol))]
    } else {
      ## fractionally displaced coordinates
      ystep <- scale[1L] * X$ystep
      xstep <- scale[2L] * X$xstep
      ystart <- X$yrow[1L] + phase[1L] * X$ystep
      xstart <- X$xcol[2L] + phase[2L] * X$xstep
      yrow <- seq(from=ystart, by=ystep, to=X$yrow[length(X$yrow)])
      xcol <- seq(from=xstart, by=xstep, to=X$xcol[length(X$xcol)])
    }
    Y <- owinInternalMask(xy=list(x=xcol, y=yrow),
                          mask=matrix(FALSE, length(yrow), length(xcol)),
                          unitname=unitname(X))
  }
  if(is.im(X)) {
    # resample an image
    if(!is.im(Y))
      Y <- as.im(Y)
    Xtype <- X$type
    Xv    <- X$v
    # handle factor-valued image as integer
    if(Xtype == "factor") 
      Xv <- array(as.integer(Xv), dim=X$dim)
    # resample
    naval <- switch(Xtype,
                    factor=,
                    integer= NA_integer_, 
                    logical = as.logical(NA_integer_), 
                    real = NA_real_, 
                    complex = NA_complex_, 
                    character = NA_character_,
                    NA)
    Y$v <- matrixsample(Xv, Y$dim, phase=phase, scale=scale, na.value=naval)
    # inherit pixel data type from X
    Y$type <- Xtype
    if(Xtype == "factor") {
      lev <- levels(X)
      Y$v <- factor(Y$v, labels=lev, levels=seq_along(lev))
      dim(Y$v) <- Y$dim
    }
  } else {
    # resample a mask
    if(!is.mask(Y)) Y <- as.mask(Y)
    Y$m <- matrixsample(X$m, Y$dim, phase=phase, scale=scale, na.value=FALSE)
  }
  return(Y)
}

pointgrid <- function(W, ngrid) {
  W <- as.owin(W)
  masque <- as.mask(W, dimyx=ngrid)
  rxy <- rasterxy.mask(masque, drop=TRUE)
  xx <- rxy$x
  yy <- rxy$y
  return(ppp(xx, yy, W))
}

onecolumn <- function(m) {
  switch(markformat(m),
         none=stop("No marks provided"),
         vector=m,
         dataframe=m[,1, drop=TRUE],
         NA)
}


checkbigmatrix <- function(n, m, fatal=FALSE, silent=FALSE) {
  nm <- as.numeric(n) * as.numeric(m)
  if(nm <= spatstat.options("maxmatrix"))
    return(TRUE)
  whinge <- paste("Attempted to create binary mask with",
                  n, "*", m, "=", nm, "entries")
  if(fatal) stop(whinge, call.=FALSE)
  if(!silent) warning(whinge, call.=FALSE)
  return(FALSE)
}


## ........... progress reports .....................

progressreport <- local({

  Put <- function(name, value, state) {
    if(is.null(state)) {
      putSpatstatVariable(paste0("Spatstat.", name), value)
    } else {
      state[[name]] <- value
    }
    return(state)
  }
  Get <- function(name, state) {
    if(is.null(state)) {
      value <- getSpatstatVariable(paste0("Spatstat.", name))
    } else {
      value <- state[[name]] 
    }
    return(value)
  }
  Exists <- function(name, state) {
    if(is.null(state)) {
      answer <- existsSpatstatVariable(paste0("Spatstat.", name))
    } else {
      answer <- name %in% names(state)
    }
    return(answer)
  }
 

  IterationsPerLine <- function(charsperline, n, every, tick,
                                showtimeinline, showevery) {
    ## Calculate number of iterations that triggers a newline.
    ## A dot is printed every 'tick' iterations
    ## Iteration number is printed every 'every' iterations.
    ##
    ## Number of characters in each report of the iteration number
    chars.report <- max(1, ceiling(log10(n)))
    chars.punctu <- if(every == 1) nchar(', ') else 0
    chars.report <- chars.report + chars.punctu
    if(showtimeinline) {
      ## If showtimeinline=TRUE, the time remaining is shown in brackets
      ## every 'showevery' iterations, where showevery \in {1, every, n}.
      ## If showtimeinline=FALSE, either the time remaining is never shown,
      ## or time remaining + estimated finish are displayed on a separate line.
      chars.time <- nchar(' [12:00:00 remaining] ')
      timesperreport <- if(showevery == 1) every else
                        if(showevery == every) 1 else 0
      chars.report <- chars.report + timesperreport * chars.time
    }
    ## Total number of characters in a complete block between iteration numbers
    chars.ticks <- floor((every-1)/tick)
    chars.block <- chars.report + chars.ticks
    ## Number of whole blocks per line
    nblocks <- max(1, floor(charsperline/chars.block))
    ## Number of iterations per line
    nperline <- nblocks * every
    ## Adjust
    leftover <- charsperline - nblocks * chars.block
    if(leftover > 0)
      nperline <- nperline + min(leftover * tick, every - 1, showevery - 1)
    ## iteration number that triggers newline
    return(nperline)
  }
  
  progressreport <- function(i, n,
                             every=min(100,max(1, ceiling(n/100))),
                             tick=1,
                             nperline=NULL,
                             charsperline=getOption("width"),
                             style=spatstat.options("progress"),
                             showtime=NULL,
                             state=NULL,
                             formula=(time ~ i),
                             savehistory=FALSE) {
    missevery <- missing(every)
    nperline.fixed <- !is.null(nperline)
    showtime.optional <- is.null(showtime)
    if(showtime.optional) showtime <- FALSE # initialise only
    if(i > n) {
      warning(paste("progressreport called with i =", i, "> n =", n))
      return(invisible(NULL))
    }
    if(style == "tk" && !requireNamespace("tcltk")) {
      warning("tcltk is unavailable; switching to style='txtbar'", call.=FALSE)
      style <- "txtbar"
    }
    if(is.null(state) && style != "tty")
      stop(paste("Argument 'state' is required when style =",sQuote(style)),
           call.=FALSE)
    ## determine model for extrapolation of time
    if(missing(formula)) formula <- NULL
    linear <- is.null(formula)
    if(!linear) {
      if(!inherits(formula, "formula"))
        stop(paste("Argument", sQuote("formula"), "should be a model formula"),
             call.=FALSE)
      savehistory <- TRUE
    }
    ## get current time
    if(savehistory || style == "tty")
      now <- proc.time()
    if(savehistory) {
      ahora <- as.numeric(now[3])
      if(i == 1) {
        state <- Put("History", data.frame(i=i, time=ahora), state)
      } else {
        history <- Get("History", state)
        history <- rbind(history, data.frame(i=i, time=ahora))
        state <- Put("History", history, state)
      }
    }
    ## display progress
    fallback <- FALSE
    switch(style,
           txtbar={
             if(i == 1) {
               ## initialise text bar
               state <- Put("ProgressBar",
                            txtProgressBar(1, n, 1, style=3),
                            state)
             } else {
               ## get text bar
               pbar <- Get("ProgressBar", state)
               ## update 
               setTxtProgressBar(pbar, i)
               if(i == n) {
                 close(pbar)
                 state <- Put("ProgressBar", NULL, state)
               } 
             }
           },
           tk={
             requireNamespace("tcltk")
             if(i == 1) {
               ## initialise text bar
               state <- Put("ProgressBar",
                            tcltk::tkProgressBar(title="progress",
                                                 min=0, max=n, width=300),
                            state)
             } else {
               ## get text bar
               pbar <- Get("ProgressBar", state)
               ## update 
               tcltk::setTkProgressBar(pbar, i,
                                       label=paste0(round(100 * i/n), "%"))
               if(i == n) {
                 close(pbar)
                 state <- Put("ProgressBar", NULL, state)
               } 
             }
           },
           tty={
             if(i == 1 || !Exists("ProgressData", state)) {
               ## Initialise stuff
               starttime   <- now
               lastnewline <- 0
               if(missevery && every > 1 && n > 10) 
                 every <- niceround(every)
               showevery <- if(showtime) every else n
               if(!nperline.fixed) 
                 nperline <- IterationsPerLine(charsperline, n, every, tick,
                                               showtime, showevery)
             } else {
               ## Extract information from previous state
               pd <- Get("ProgressData", state)
               if(is.null(pd))
                 stop(paste("progressreport called with i =", i,
                            "before i = 1"))
               every        <- pd$every
               tick         <- pd$tick
               nperline     <- pd$nperline
               lastnewline  <- pd$lastnewline
               starttime    <- pd$starttime
               showtime     <- pd$showtime
               showevery    <- pd$showevery
               showtime.optional <- pd$showtime.optional
               nperline.fixed    <- pd$nperline.fixed
               if(i < n) {
                 if(showtime || showtime.optional) {
                   ## estimate time remaining
                   elapsed <- now - starttime
                   elapsed <- unname(elapsed[3])
                   if(linear) {
                     rate <- elapsed/(i-1)
                     remaining <- rate * (n-i)
                   } else {
                     fit <- try(lm(formula, data=history))
                     ok <- !inherits(fit, "try-error") && !anyNA(coef(fit))
                     if(ok) {
                       pred <- suppressWarnings(
                         predict(fit, newdata=data.frame(i=c(i, i+1, n)))
                       )
                       ok <- all(diff(pred) >= 0)
                     }
                     if(ok) {
                       ## predictions of model 
                       remaining <- pred[3] - pred[1]
                       rate <- pred[2] - pred[1]
                     } else {
                       ## linear extrapolation
                       fallback <- TRUE
                       rate <- elapsed/(i-1)
                       remaining <- rate * (n-i)
                     }
                   }
                   if(!showtime) {
                     ## Currently not showing the time remaining.
                     ## Change this if:
                     if(rate > 20) {
                       ## .. more than 20 seconds until next iteration
                       showtime <- TRUE
                       showevery <- 1
                     } else if(remaining > 180) {
                       ## ... more than 3 minutes remaining
                       showtime <- TRUE
                       showevery <- every
                       aminute <- ceiling(60/rate)
                       if(aminute < showevery) 
                         showevery <- min(niceround(aminute), showevery)
                     }
                     # update number of iterations per line
                     if(showtime && !nperline.fixed) {
                       showtimeinline <- (remaining < 600)
                       nperline <- IterationsPerLine(charsperline,
                                                     n, every, tick,
                                                     showtimeinline,
                                                     showevery)
                     }
                   }
                 }
               }
             }
             ## determine whether newline is required
             offset <- if(lastnewline == 0 && every != 1) 6 else 0
             do.newline <- ((i - lastnewline + offset) %% nperline == 0)
             ## Finally, print the report
             if(i == n) {
               cat(paste0("\n", n, ".\n"))
             } else if(every == 1 || i <= 3) {
               cat(paste0(i, ",", if(do.newline) "\n" else " "))
             } else {
               if(i %% every == 0) 
                 cat(i)
               else if(i %% tick == 0)
                 cat(".")
               if(do.newline)
                 cat("\n")
             }
             if(showtime && i > 1 && i < n && (i %% showevery == 0)) {
               st <- paste(codetime(round(remaining)),
                           paste0("remaining",
                                  if(fallback) "(linear)" else ""))
               if(longwait <- (remaining > 600)) {
                 finishtime <- Sys.time() + remaining
                 st <- paste0(st, ", estimate finish ", round(finishtime))
                 do.newline <- TRUE
               }
               st <- paren(st, "[")
               brk <- if(longwait) "\n" else " "
               cat(paste0(brk, st, brk))
             }
             ## remember when the last newline occurred
             if(do.newline)
               lastnewline <- i
             ## save the current state
             state <- Put("ProgressData", 
                          list(every=every,
                               tick=tick,
                               nperline=nperline,
                               lastnewline=lastnewline,
                               starttime=starttime,
                               showtime=showtime,
                               showevery=showevery,
                               nperline.fixed=nperline.fixed,
                               showtime.optional=showtime.optional),
                          state)
             flush.console()
           },
           stop(paste("Unrecognised option for style:", dQuote(style)))
           )
    return(invisible(state))
  }

  progressreport
})

## .... special tweaks .........

multiply.only.finite.entries <- function(x, a) {
  # In ppm a potential value that is -Inf must remain -Inf
  # and a potential value that is 0 multiplied by NA remains 0
  y <- x
  ok <- is.finite(x) & (x != 0)
  y[ok] <- a * x[ok]
  return(y)
}
 
## print names and version numbers of libraries loaded

sessionLibs <- local({

  sessionLibs <- function() {
    a <- sessionInfo()
    d1 <- mangle(a$otherPkgs, "loaded")
    d2 <- mangle(a$loadedOnly, "imported")
    return(invisible(list(loaded=d1,imported=d2)))
  }

  mangle <- function(pkglist, type="loaded") {
    if(length(pkglist)) {
      b <- unlist(lapply(pkglist, getElement, name="Version"))
      b <- b[order(names(b))]
      g <- rbind(names(b), unname(b))
      d <- apply(g, 2, paste, collapse=" ")
    } else d <- NULL
    if(length(d) > 0) {
      cat(paste0("Libraries ", type, ":\n"))
      for(di in d) cat(paste("\t", di, "\n"))
    } else cat(paste0("Libraries ", type, ": none\n"))
    return(invisible(d))
  }
    
  sessionLibs
})



# ..................

prepareTitle <- function(main) {
  ## Count the number of lines in a main title
  ## Convert title to a form usable by plot.owin
  if(is.expression(main)) {
    nlines <- 1
  } else {
    main <- paste(main)
    ## break at newline 
    main <- unlist(strsplit(main, "\n"))
    nlines <- if(sum(nchar(main)) == 0) 0 else length(main)
  }
  return(list(main=main,
              nlines=nlines,
              blank=rep('  ', nlines)))
}

requireversion <- function(pkg, ver, fatal=TRUE) {
  pkgname <- deparse(substitute(pkg))
  pkgname <- gsub("\"", "", pkgname)
  pkgname <- gsub("'", "", pkgname)
  dfile <- system.file("DESCRIPTION", package=pkgname)
  if(nchar(dfile) == 0) {
    ## package is not installed
    if(!fatal) return(FALSE) else 
    stop(paste("Package", sQuote(pkgname), "is needed but is not installed"),
         call.=FALSE)
  }
  v <- read.dcf(file=dfile, fields="Version")
  ok <- (package_version(v) >= ver)
  if(!ok && fatal) 
    stop(paste("Package",
               sQuote(pkgname),
               "is out of date: version >=",
               ver,
               "is needed"),
         call.=FALSE)
  return(if(ok) invisible(TRUE) else FALSE)
}

spatstatDiagnostic <- function(msg) {
  cat("-----------------------------\n")
  cat(paste(" >>> Spatstat Diagnostic: ", msg, "<<<\n"))
  cat("-----------------------------\n")
  invisible(NULL)
}

allElementsIdentical <- function(x, entry=NULL) {
  if(length(x) <= 1) return(TRUE)
  if(is.null(entry)) {
    x1 <- x[[1]]
    for(i in 2:length(x))
      if(!identical(x[[i]], x1)) return(FALSE)
  } else {
    e1 <- x[[1]][[entry]]
    for(i in 2:length(x))
      if(!identical(x[[i]][[entry]], e1)) return(FALSE)
  }
  return(TRUE)
}

resolve.stringsAsFactors <- function(stringsAsFactors=NULL) {
  if(is.null(stringsAsFactors) || is.na(stringsAsFactors)) {
    if(getRversion() < "4.1.0") default.stringsAsFactors() else FALSE
  } else isTRUE(stringsAsFactors) 
}  
