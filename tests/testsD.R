#'
#'   Header for all (concatenated) test files
#'
#'   Require spatstat.geom
#'   Obtain environment variable controlling tests.
#'
#'   $Revision: 1.5 $ $Date: 2020/04/30 05:31:37 $

require(spatstat.geom)
FULLTEST <- (nchar(Sys.getenv("SPATSTAT_TEST", unset="")) > 0)
ALWAYS   <- TRUE
cat(paste("--------- Executing",
          if(FULLTEST) "** ALL **" else "**RESTRICTED** subset of",
          "test code -----------\n"))
#'
#'   tests/dominic.R
#'
#'   Additional tests for Dominic Schuhmacher's code
#'
#'   $Revision: 1.5 $  $Date: 2020/12/03 03:23:25 $

if(ALWAYS) {   # tests C code
local({
  X <- runifrect(10)
  Y <- runifrect(10)

  d  <- pppdist(X, Y, type="ace", show.rprimal=TRUE)
  a <- matchingdist(d, type="ace")
  b <- matchingdist(d, type="mat")

  d2 <- pppdist(X, Y, type="spa", ccode=FALSE)
  d2 <- pppdist(X, Y, type="spa", ccode=TRUE, auction=FALSE)
  d3 <- pppdist(X, Y, type="mat", ccode=TRUE, auction=FALSE)
  d4 <- pppdist(X[FALSE], Y[FALSE], matching=TRUE, type="spa")
  d4 <- pppdist(X[FALSE], Y[FALSE], matching=FALSE, type="spa")
  d4 <- pppdist(X[FALSE], Y[FALSE], matching=TRUE, type="ace")
  d4 <- pppdist(X[FALSE], Y[FALSE], matching=FALSE, type="ace")

  m  <- pppdist.mat(X, Y, q=Inf, cutoff=0.001)
  m2 <- pppdist.mat(X[FALSE], Y[FALSE], q=Inf, cutoff=0.001)
  m3 <- pppdist.mat(X[FALSE], Y[FALSE], q=2, cutoff=0.001)

})
}



#'
#'  tests/discarea.R
#'
#'   $Revision: 1.3 $ $Date: 2020/04/28 12:58:26 $
#'

if(ALWAYS) {
local({
  u <- c(0.5,0.5)
  B <- owin(poly=list(x=c(0.3, 0.5, 0.7, 0.4), y=c(0.3, 0.3, 0.6, 0.8)))
  areaGain(u, cells, 0.1, exact=TRUE)
  areaGain(u, cells, 0.1, W=NULL)
  areaGain(u, cells, 0.1, W=B)

  X <- cells[square(0.4)]
  areaLoss(X, 0.1, exact=TRUE)  # -> areaLoss.diri
  areaLoss(X, 0.1, exact=FALSE) # -> areaLoss.grid
  areaLoss.poly(X, 0.1)

  areaLoss(X, 0.1, exact=FALSE, method="distmap")          # -> areaLoss.grid
  areaLoss(X, c(0.1, 0.15), exact=FALSE, method="distmap") # -> areaLoss.grid
})
}
#'
#'   tests/duplicity.R
#'
#'  Tests of duplicated/multiplicity code
#'
#' $Revision: 1.8 $ $Date: 2020/04/28 12:58:26 $

if(ALWAYS) {
local({
   X <- ppp(c(1,1,0.5,1), c(2,2,1,2), window=square(3), check=FALSE)
   Y <- X %mark% factor(letters[c(3,2,4,3)])
   ZC <- X %mark% letters[c(3,2,4,3)]
   ZM <- Y %mark% matrix(c(3,2,4,3), 4, 2)
   ZD <- Y %mark% as.data.frame(marks(ZM))

   #' multiplicity
   m <- multiplicity(X)
   mf <- multiplicity(Y)
   mm <- multiplicity(ZM)
   mz <- multiplicity(ZD)
   mc <- multiplicity(ZC)
   ## default method
   kk <- c(1,2,3,1,1,2)
   mk <- multiplicity(kk)
   ml <- multiplicity(list(sin, cos, tan)[kk])
   mc <- multiplicity(c("sin", "cos", "tan")[kk])
   if(!identical(ml, mk))
     stop("multiplicity.default(<list>) disagrees with multiplicityNumeric")
   if(!identical(mc, mk))
     stop("multiplicity(<character>) disagrees with multiplicity(<numeric>)")
   ## data frame method
   df <- data.frame(x=c(1:4, 1,3,2,4, 0,0, 3,4),
                    y=factor(rep(letters[1:4], 3)))
   md <- multiplicity(df)

   ## uniquemap.ppp
   checkum <- function(X, blurb) {
     a <- uniquemap(X)
     if(any(a > seq_along(a)))
       stop(paste("uniquemap", blurb,
                  "does not respect sequential ordering"))
     return(invisible(NULL))
   }
   checkum(X, "<unmarked point pattern>")
   checkum(Y, "<multitype point pattern>")
   checkum(ZC, "<point pattern with character marks>")
   checkum(ZM, "<point pattern with matrix of marks>")
   checkum(ZD, "<point pattern with several columns of marks>")

   ## uniquemap.data.frame
   dfbase <- as.data.frame(replicate(3, sample(1:20, 10), simplify=FALSE))
   df <- dfbase[sample(1:10, 30, replace=TRUE), , drop=FALSE]
   #' faster algorithm for numeric values
   checkum(df, "<numeric data frame>")
   a <- uniquemap(df)
   #' general algorithm using 'duplicated' and 'match'
   dfletters <- as.data.frame(matrix(letters[as.matrix(df)], nrow=nrow(df)))
   checkum(dfletters, "<character data frame>")
   b <- uniquemap(dfletters)
   if(!isTRUE(all.equal(a,b)))
     stop("inconsistency between algorithms in uniquemap.data.frame")

   ## uniquemap.matrix
   M0 <- matrix(1:12, 3, 4)
   ii <- sample(1:3, 5, replace=TRUE)
   M4 <- M0[ii, , drop=FALSE]
   checkum(M4, "<integer matrix>")
   u4 <- uniquemap(M4)
   C4 <- matrix(letters[M4], 5, 4)
   uc4 <- uniquemap(C4)
   checkum(C4, "<character matrix>")
   if(!isTRUE(all.equal(u4, uc4)))
     stop("Inconsistency between algorithms in uniquemap.matrix")
   
   ## uniquemap.default
   a <- letters[c(1, 1:4, 3:2)]
   checkum(a, "<character>")
   checkum(as.list(a), "<list>")
   u1 <- uniquemap(a)
   u2 <- uniquemap(as.list(a))
   if(!isTRUE(all.equal(u1, u2)))
     stop("Inconsistency between algorithms in uniquemap.default")
})
}
