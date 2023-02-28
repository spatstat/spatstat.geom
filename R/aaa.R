#'
#'    aaa.R
#'
#'   Code that must be read before the rest of the R code in spatstat.geom
#' 
#'    $Revision: 1.2 $  $Date: 2023/02/28 03:36:53 $

.spEnv <- new.env()

putSpatstatVariable <- function(name, value) {
  assign(name, value, envir=.spEnv)
}
getSpatstatVariable <- function(name, default=NULL) {
  if(exists(name, envir=.spEnv)) get(name, envir=.spEnv) else default 
}
existsSpatstatVariable <- function(name) {
  exists(name, envir=.spEnv)
}

