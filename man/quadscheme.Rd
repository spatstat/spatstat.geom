\name{quadscheme}
\alias{quadscheme}
\title{Generate a Quadrature Scheme from a Point Pattern}
\description{
  Generates a quadrature scheme (an object of class \code{"quad"})
  from point patterns of data and dummy points.
}
\usage{
 quadscheme(data, dummy, method="grid", \dots)
}
\arguments{
  \item{data}{
    The observed data point pattern.
    An object of class \code{"ppp"}
    or in a format recognised by \code{\link{as.ppp}()}
  }
  \item{dummy}{
    The pattern of dummy points for the quadrature. 
    An object of class \code{"ppp"}
    or in a format recognised by \code{\link{as.ppp}()}
    Defaults to \code{default.dummy(data, ...)}
  }
  \item{method}{
    The name of the method for calculating quadrature weights: either
    \code{"grid"} or \code{"dirichlet"}.
  }
  \item{\dots}{
    Parameters of the weighting method (see below)
    and parameters for constructing the dummy points if necessary.
  }
} 
\value{
  An object of class \code{"quad"} describing the quadrature scheme
  (data points, dummy points, and quadrature weights)
  suitable as the argument \code{Q} of the function \code{\link[spatstat.model]{ppm}()} for
  fitting a point process model.

  The quadrature scheme can be inspected using the
  \code{print} and \code{plot} methods for objects
  of class \code{"quad"}. 
}
\details{
  This is the primary method for producing a quadrature schemes
  for use by \code{\link[spatstat.model]{ppm}}.
  
  The function \code{\link[spatstat.model]{ppm}} fits a point process model to an
  observed point pattern using 
  the Berman-Turner quadrature approximation (Berman and Turner, 1992;
  Baddeley and Turner, 2000) to the pseudolikelihood of the model. 
  It requires a quadrature scheme consisting of 
  the original data point pattern, an additional pattern of dummy points,
  and a vector of quadrature weights for all these points.
  Such quadrature schemes are represented by objects of class
  \code{"quad"}. See \code{\link{quad.object}} for a description of this class.

  Quadrature schemes are created by the function
  \code{quadscheme}.
  The arguments \code{data} and \code{dummy} specify the data and dummy
  points, respectively. There is a sensible default for the dummy 
  points (provided by \code{\link{default.dummy}}).
  Alternatively the dummy points
  may be specified arbitrarily and given in any format recognised by
  \code{\link{as.ppp}}.
  There are also functions for creating dummy patterns
  including \code{\link{corners}},
  \code{\link{gridcentres}},
  \code{\link{stratrand}} and
  \code{\link{spokes}}.
 
  The quadrature region is the region over which we are
  integrating, and approximating integrals by finite sums.
  If \code{dummy} is a point pattern object (class \code{"ppp"})
  then the quadrature region is taken to be \code{Window(dummy)}.
  If \code{dummy} is just a list of \eqn{x, y} coordinates
  then the quadrature region defaults to the observation window
  of the data pattern, \code{Window(data)}.

  If \code{dummy} is missing, then a pattern of dummy points
  will be generated using \code{\link{default.dummy}}, taking account
  of the optional arguments \code{...}.
  By default, the dummy points are arranged in a
  rectangular grid; recognised arguments
  include \code{nd} (the number of grid points
  in the horizontal and vertical directions)
  and \code{eps} (the spacing between dummy points).
  If \code{random=TRUE}, a systematic random pattern
  of dummy points is generated instead.
  See \code{\link{default.dummy}} for details. 
  
  If \code{method = "grid"} then the optional arguments (for \code{\dots}) are 
  \code{(nd, ntile, eps)}.
  The quadrature region (defined above) is divided into
  an \code{ntile[1]} by \code{ntile[2]} grid of rectangular tiles.
  The weight for each
  quadrature point is the area of a tile divided by the number of
  quadrature points in that tile. 
 
  If \code{method="dirichlet"} then the optional arguments are
  \code{(exact=TRUE, nd, eps)}.
  The quadrature points (both data and dummy) are used to construct the
  Dirichlet tessellation. The quadrature weight of each point is the
  area of its Dirichlet tile inside the quadrature region.
  If \code{exact == TRUE} then this area is computed exactly
  using the package \code{deldir}; otherwise it is computed
  approximately by discretisation.
}
\section{Error Messages}{
  The following error messages need some explanation.
  (See also the list of error messages in \code{\link[spatstat.model]{ppm.ppp}}).
  \itemize{
    \item{\dQuote{Some tiles with positive area do not contain any
	quadrature points: relative error = X\%}}{
      This is not important unless the relative error is large.
      In the default rule for computing the
      quadrature weights, space is divided into rectangular tiles,
      and the number of quadrature points (data and dummy points) in
      each tile is counted. It is possible for a tile with non-zero area
      to contain no quadrature points; in this case, the quadrature
      scheme will contribute a bias to the model-fitting procedure.
      \bold{A small relative error (less than 2 percent) is not	important.}
      Relative errors of a few percent can occur because of the shape of
      the window.
      If the relative error is greater than about 5 percent, we
      recommend trying different parameters for the quadrature scheme,
      perhaps setting a larger value of \code{nd} to increase the number
      of dummy points. A relative error greater than 10 percent
      indicates a major problem with the input data. The quadrature
      scheme should be inspected by plotting and printing it.
      (The most likely cause of this problem is that the spatial coordinates
      of the original data were not handled correctly, for example,
      coordinates of the locations and the window boundary were incompatible.)
    }
    \item{\dQuote{Some tiles with zero area contain quadrature points}}{
      This error message is rare, and has no consequences.
      It is mainly of interest to
      programmers. It occurs when the area of a
      tile is calculated to be equal to zero, but a quadrature point
      has been placed in the tile. 
    }
  }
}
\references{
  Baddeley, A. and Turner, R.
  Practical maximum pseudolikelihood for spatial point patterns.
  \emph{Australian and New Zealand Journal of Statistics}
  \bold{42} (2000) 283--322.
 
  Berman, M. and Turner, T.R. 
  Approximating point process likelihoods with GLIM.
  \emph{Applied Statistics} \bold{41} (1992) 31--38.
}
\seealso{
  \code{\link[spatstat.model]{ppm}},
  \code{\link{as.ppp}},
  \code{\link{quad.object}},
  \code{\link{gridweights}},
  \code{\link{dirichletWeights}},
  \code{\link{corners}},
  \code{\link{gridcentres}},
  \code{\link{stratrand}},
  \code{\link{spokes}}
}
\examples{
  data(simdat)

  # grid weights
  Q <- quadscheme(simdat)
  Q <- quadscheme(simdat, method="grid")
  Q <- quadscheme(simdat, eps=0.5)         # dummy point spacing 0.5 units

  Q <- quadscheme(simdat, nd=50)           # 1 dummy point per tile
  Q <- quadscheme(simdat, ntile=25, nd=50) # 4 dummy points per tile

  # Dirichlet weights
  Q <- quadscheme(simdat, method="dirichlet", exact=FALSE)

  # random dummy pattern
  # D <- runifrect(250, Window(simdat))
  # Q <- quadscheme(simdat, D, method="dirichlet", exact=FALSE)

  # polygonal window
  data(demopat)
  X <- unmark(demopat)
  Q <- quadscheme(X)

  # mask window
  Window(X) <- as.mask(Window(X))
  Q <- quadscheme(X)
  
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{datagen}
