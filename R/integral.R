#'
#'   integral.R
#'
#'   generic
#'
#'   $Revision: 1.1 $ $Date: 2023/10/22 02:03:35 $

integral <- function(f, domain=NULL, ...) {
  UseMethod("integral")
}

