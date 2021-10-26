##  spatstat.geom/R/First.R

.onLoad <- function(...) {
  reset.spatstat.options()
}

.onAttach <- function(libname, pkgname) {
  vs <- read.dcf(file=system.file("DESCRIPTION", package="spatstat.geom"),
                 fields="Version")
  vs <- as.character(vs)
  putSpatstatVariable("SpatstatGeomVersion", vs)
  packageStartupMessage(paste("spatstat.geom", vs))
  return(invisible(NULL))
}

  
