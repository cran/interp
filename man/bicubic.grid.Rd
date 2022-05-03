\name{bicubic.grid}
\alias{bicubic.grid}
\title{
  Bicubic Interpolation for Data on a Rectangular grid
}
\description{
This is a placeholder function for backward compatibility with packaga akima.

In its current state it simply calls the reimplemented Akima algorithm for
irregular grids applied to the regular gridded data given.

Later a reimplementation of the original algorithm for regular grids may follow.
}
\usage{
bicubic.grid(x,y,z,xlim=c(min(x),max(x)),ylim=c(min(y),max(y)),
             nx=40,ny=40,dx=NULL,dy=NULL)
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
  \item{xlim}{
    vector of length 2 giving lower and upper limit for range \code{x}
    coordinates used for  output grid.
  }
  \item{ylim}{
    vector of length 2 giving lower and upper limit for range of \code{y}
    coordinates used for  output grid.
  }
  \item{nx}{
    output grid dimension in \code{x} direction.
    }
  \item{ny}{
    output grid dimension in \code{y} direction.
    }
  \item{dx}{
    output grid spacing in \code{x} direction, not used by default,
  overrides \code{nx} if specified.
    }
  \item{dy}{
    output grid spacing in \code{y} direction, not used by default,
    overrides \code{ny} if specified..
    }
}
\details{
  This function is a call wrapper for backward compatibility with package akima.

  Currently it applies Akimas irregular grid splines to regular grids, later 
  a FOSS reimplementation of his regular grid splines may replace this wrapper.
}
\value{
  This function produces a grid of interpolated points, feasible to be
  used directly with \code{\link{image}} and \code{\link{contour}}:
  \item{x}{vector of \code{x} coordinates of the output grid.}
  \item{y}{vector of \code{y} coordinates of the output grid.}
  \item{z}{matrix of interpolated data for the output grid.}
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
  \code{\link{interp}}, \code{\link{bicubic}}
  % maybe later:
  % \code{\link[rgeostat]{bilinear}}
}
\examples{
data(akima474)
# interpolate at a grid [0,8]x[0,10]
akima.bic <- bicubic.grid(akima474$x,akima474$y,akima474$z)
zmin <- min(akima.bic$z, na.rm=TRUE)
zmax <- max(akima.bic$z, na.rm=TRUE)
breaks <- pretty(c(zmin,zmax),10)
colors <- heat.colors(length(breaks)-1)
image(akima.bic, breaks=breaks, col=colors)
contour(akima.bic, levels=breaks,  add=TRUE)
}
\keyword{ dplot }

