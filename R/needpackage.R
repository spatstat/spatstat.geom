#'
#'   needpackage.R
#'
#'   Check that a package is present
#'
#'   $Revision: 1.1 $ $Date: 2025/05/31 02:55:03 $

needpackage <- function(pkg, loaded=FALSE, fatal=TRUE, purpose=NULL) {
  if(!requireNamespace(pkg, quietly=TRUE)) {
    if(fatal)
      stop(paste("The package", sQuote(pkg), "is required", purpose),
           call.=FALSE)
    return(FALSE)
  }
  if(loaded && !isNamespaceLoaded(pkg)) {
    if(fatal) 
      stop(paste("The package", sQuote(pkg), "must be loaded",
                 paste0(purpose, ";"),
                 "please type", sQuote(paste0("library", paren(pkg)))),
           call.=FALSE)
    return(FALSE)
  }
  return(TRUE)
}

ensureModelSupport <- function(model, loaded=TRUE) {
  if(inherits(model, c("ppm", "kppm", "dppm", "slrm", "rppm")))
    needpackage("spatstat.model", loaded=loaded)
  if(inherits(model, "lppm"))
    needpackage("spatstat.linnet", loaded=loaded)
  return(TRUE)
}
