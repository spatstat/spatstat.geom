#'
#'   flaky.R
#'
#'  Code to defend against crashes in calls to external packages
#'
#'  $Revision: 1.1 $ $Date: 2023/07/09 04:13:48 $
#'

safeDevCapabilities <- function() {
  a <- try(dev.capabilities())
  if(!inherits(a, "try-error"))
    return(a)
  warning("dev.capabilities() caused an error!", call.=FALSE)
  return(NULL)
}
