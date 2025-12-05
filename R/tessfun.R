#'
#'  tessfun.R
#'
#'  Functions which are constant on each tile of a tessellation
#'
#' Copyright (c) 2024 Adrian Baddeley, Tilman Davies and Martin Hazelton
#' 

as.function.tess <- function(x, ..., values=NULL) {
  V <- x
  if(is.null(values)) {
    f <- function(x,y) { tileindex(x,y,V) }
  } else {
    if(is.data.frame(values)) values <- unlist(values)
    if(length(values) != nobjects(x))
      stop("Length of 'values' should equal the number of tiles", call.=FALSE)
    values <- unname(values)
    f <- function(x,y) { values[as.integer(tileindex(x,y,V))] }
  }
  g <- funxy(f, Window(V))
  class(g) <- unique(c("tessfun", class(g)))
  return(g)
}

as.tess.tessfun <- function(X) {
  get("V", envir=environment(X))
}

tessfunvalues <- function(f) {
  get("values", envir=environment(f)) %orifnull% seq_len(nobjects(as.tess(f)))
}

integral.tessfun <- function(f, domain=NULL, ...) {
  tes <- as.tess(f)
  val <- tessfunvalues(f)
  if(is.factor(val) || is.character(val))
    stop(paste("Cannot integrate a function which returns",
               if(is.factor(val)) "factor" else "character", "values"),
         call.=FALSE)
  if(!is.complex(val)) val <- as.numeric(val) # need real or complex values
  if(!is.null(domain)) {
    marks(tes) <- val
    tes <- intersect.tess(tes, domain)
    val <- unlist(marks(tes))
  }
  sum(tile.areas(tes) * val)
}

print.tessfun <- function(x, ...) {
  splat("Function which is constant on each tile of a tessellation")
  cat("\n")
  print.funxy(x, ...)
  cat("\n")
  print.tess(as.tess(x))
  cat("\n")
  values <- tessfunvalues(x)
  if(is.factor(values)) {
    splat("Function values are categorical, with levels")
    print(levels(values))
  } else {
    splat("Function values are of type", sQuote(typeof(values)))
    if(is.numeric(values)) 
      splat("Range of function values:", prange(signif(range(values), 4)))
  }
  invisible(NULL)
}
  
plot.tessfun <- function(x, ...) {
  xname <- short.deparse(substitute(x))
  tes <- as.tess(x)
  val <- tessfunvalues(x)
  do.call(plot.tess,
          resolve.defaults(list(quote(tes)),
                           list(...),
                           list(do.col=TRUE, values=val, main=xname)))
}

levelset.tessfun <- function(X, thresh, compare="<=", ...) {
  #' extract tessellation and values attached to each tile
  Tess     <- get("V", envir=environment(X))
  values   <- tessfunvalues(X)
  selected <- switch(compare,
                     "<"  = (values < thresh),
                     ">"  = (values > thresh),
                     "<=" = (values <= thresh),
                     ">=" = (values >= thresh),
                     "==" = (values == thresh),
                     "!=" = (values  != thresh),
                     stop(paste("unrecognised comparison operator",
                                sQuote(compare))))
  W <- Window(Tess)
  if(all(selected)) return(W)
  if(!any(selected)) return(emptywindow(Frame(W)))
  result <- do.call(union.owin, tiles(Tess)[selected])
  Frame(result) <- Frame(W)
  return(result)
}


