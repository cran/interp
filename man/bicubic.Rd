\name{bicubic}
\alias{bicubic}
\title{
  Bivariate Interpolation for Data on a Rectangular grid
}
\description{
This is a placeholder function for backward compatibility with packaga akima.

In its current state it simply calls the reimplemented Akima algorithm for
irregular grids applied to the regular gridded data given.

Later a reimplementation of the original algorithm for regular grids may follow.
}
\usage{
bicubic(x, y, z, x0, y0)
}

\arguments{
  \item{x}{
    a vector containing the \code{x} coordinates of the rectangular data grid.
  }
  \item{y}{
    a vector containing the \code{y} coordinates of the rectangular data grid.
  }
  \item{z}{
    a matrix containing the \code{z[i,j]} data values for the grid points (\code{x[i]},\code{y[j]}).
  }
  \item{x0}{
    vector of \code{x} coordinates used to interpolate at.
  }
  \item{y0}{
    vector of \code{y} coordinates used to interpolate at.
  }
}
\details{
  This function is a call wrapper for backward compatibility with package akima.

  Currently it applies Akimas irregular grid splines to regular grids, later 
  a FOSS reimplementation of his regular grid splines may replace this wrapper.
}
\value{
  This function produces a list of interpolated points:
  \item{x}{vector of \code{x} coordinates.}
  \item{y}{vector of \code{y} coordinates.}
  \item{z}{vector of interpolated data \code{z}.}
  If you need an output grid, see \code{\link{bicubic.grid}}.
}
\references{
  Akima, H. (1996) Rectangular-Grid-Data
  Surface Fitting that Has the Accuracy of a
  Bicubic Polynomial,
  J. ACM \bold{22}(3), 357-361
}
\note{
  Use \code{\link{interp}} for the general case of irregular gridded data!
}

\seealso{
  \code{\link{interp}}, \code{\link{bicubic.grid}}
  % maybe later:
  % \code{\link[rgeostat]{bilinear}}
}
\examples{
data(akima474)
# interpolate at the diagonal of the grid [0,8]x[0,10]
akima.bic <- bicubic(akima474$x,akima474$y,akima474$z,
                     seq(0,8,length=50), seq(0,10,length=50))
plot(sqrt(akima.bic$x^2+akima.bic$y^2), akima.bic$z, type="l")

}

\keyword{ dplot }

