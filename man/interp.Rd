\name{interp}
\alias{interp}
\title{
  Interpolation function
}
\description{
  This function currently implements piecewise linear interpolation
  (=barycentric interpolation).
}
\usage{
interp(x, y = NULL, z, xo = seq(min(x), max(x), length = nx),
       yo = seq(min(y), max(y), length = ny),
       linear = (method == "linear"), extrap = FALSE,
       duplicate = "error", dupfun = NULL,
       nx = 40, ny = 40, input="points", output = "grid",
       method = "linear", deltri = "shull", h=0,
       kernel="gaussian", solver="QR", degree=3,
       baryweight=TRUE, autodegree=FALSE, adtol=0.1,
       smoothpde=FALSE, akimaweight=TRUE, nweight=25)
}
\arguments{
  \item{x}{
    vector of \eqn{x}-coordinates of data points or a
    \code{SpatialPointsDataFrame} object.
    Missing values are not accepted.
}
  \item{y}{
    vector of \eqn{y}-coordinates of data points.
    Missing values are not accepted.

    If left as NULL indicates that \code{x} should be a
    \code{SpatialPointsDataFrame} and \code{z} names the variable of
    interest in this dataframe.
}
  \item{z}{
    vector of \eqn{z}-values at data points or a character variable
    naming the variable of interest in the
    \code{SpatialPointsDataFrame} \code{x}.

    Missing values are not accepted.

    \code{x}, \code{y}, and \code{z} must be the same length
    (execpt if \code{x} is a \code{SpatialPointsDataFrame}) and may
    contain no fewer than four points. The points of \code{x} and
    \code{y} should not be collinear, i.e, they should not fall on the same line
    (two vectors \code{x} and \code{y} such that \code{y = ax + b} for
    some \code{a}, \code{b} will not produce menaningful results).


    \code{interp} is
    meant for cases in which you have \eqn{x}, \eqn{y} values
    scattered over a plane and a \eqn{z} value for each.  If, instead,
    you are trying to evaluate a mathematical function, or get a
    graphical interpretation of relationships that can be described by a
    polynomial, try \code{\link{outer}}.
  }
  \item{xo}{
    If \code{output="grid"} (default): sequence of \eqn{x} locations for
    rectangular output grid, defaults to \code{nx} points between
    \code{min(x)} and \code{max(x)}.

    If \code{output="points"}: vector of \eqn{x} locations for output points.
  }
  \item{yo}{
    If \code{output="grid"} (default): sequence of \eqn{y} locations for
    rectangular output grid, defaults to \code{ny} points between
    \code{min(y)} and \code{max(y)}.

    If \code{output="points"}: vector of \eqn{y} locations for output
    points. In this case it has to be same length as \code{xo}.
  }
  \item{input}{
    text, possible values are \code{"grid"} (not yet implemented) and
    \code{"points"} (default).

    This is used to distinguish between regular and irregular gridded data.
  }
  \item{output}{
    text, possible values are \code{"grid"} (=default) and
    \code{"points"}.

    If \code{"grid"} is choosen then \code{xo} and \code{yo} are
    interpreted as vectors spanning a rectangular grid of points
    \eqn{(xo[i],yo[j])}, \eqn{i=1,...,nx}, \eqn{j=1,...,ny}. This
    default behaviour matches how \code{akima::interp} works.

    In the case of \code{"points"}  \code{xo} and \code{yo} have to be
    of same length and are taken as possibly irregular spaced output
    points \eqn{(xo[i],yo[i])}, \eqn{i=1,...,no} with
    \code{no=length(xo)}. \code{nx} and \code{ny} are ignored in this
    case. This case is meant as replacement for the pointwise
    interpolation done by \code{akima::interpp}. If the input \code{x}
    is a \code{SpatialPointsDataFrame} and \code{output="points"} then
    \code{xo} has to be a \code{SpatialPointsDataFrame}, \code{yo} will
    be ignored.
  }
  \item{linear}{
    logical, only for backward compatibility with \code{akima::interp},
    indicates if piecewise linear interpolation or Akima splines should be
    used. 

    Please use the new \code{method} argument instead!
  }
  \item{method}{
    text, possible methods are (currently only, more is under
    developement) \code{"linear"} (piecewise linear interpolation within
    the triangles of the Delauney triangulation, also
    referred to as barycentric interpolation based on barycentric 
    coordinates) and \code{"akima"} (a reimplementation for Akimas spline 
    algorithms for irregular gridded data with the accuracy of a 
    bicubic polynomial).

    This replaces the old \code{linear} argument of \code{akima::interp}.
  }
  \item{extrap}{
    logical, indicates if extrapolation outside the convex hull is
    intended, will not work for piecewise linear interpolation!
  }
  \item{duplicate}{
    character string indicating how to handle duplicate
    data points. Possible values are
    \describe{
      \item{\code{"error"}}{produces an error message,}
      \item{\code{"strip"}}{remove duplicate z values,}
      \item{\code{"mean"},\code{"median"},\code{"user"}}{calculate
	mean , median or user defined function (\code{dupfun}) of duplicate
	\eqn{z} values.}
    }
  }
  \item{dupfun}{
    a function, applied to duplicate points if
    \code{duplicate= "user"}.}
  \item{nx}{
    dimension of output grid in x direction
  }
  \item{ny}{
    dimension of output grid in y direction
  }
  \item{deltri}{
    triangulation method used, this argument will later be moved
    into a control set together with others related to the spline
    interpolation! Possible values are \code{"shull"} (default, sweep 
    hull algorithm) and \code{"deldir"} (uses package\code{deldir}).
  }
  \item{h}{bandwidth for partial derivatives estimation, compare \code{\link{locpoly}} for details}

  \item{kernel}{kernel for partial derivatives estimation, compare \code{\link{locpoly}} for details}

  \item{solver}{solver used in partial derivatives estimation, compare \code{\link{locpoly}} for details}
  \item{degree}{degree of local polynomial used for partial derivatives
    estimation, compare \code{\link{locpoly}} for details}
  \item{baryweight}{calculate three partial derivatives estimatiors and
    weight them}
  \item{autodegree}{try to reduce degree automatically}
  \item{adtol}{used for autodegree}
  \item{smoothpde}{Use a an averaged version of partial derivatives
    estimates, by default simple average of \code{nweight} estimates.

    Currently disabled by default (FALSE), underlying code still a bit 
    experimental.
  }

  \item{akimaweight}{apply Akima weighting scheme on partial derivatives
    estimations instead of simply averaging}
  \item{nweight}{size of search neighbourhood for weighting scheme,
    default: 25}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  a list with 3 components:
  \item{x,y}{
    If \code{output="grid"}:
    vectors of \eqn{x}- and \eqn{y}-coordinates of output grid, the same
    as the input
    argument \code{xo}, or \code{yo}, if present.  Otherwise, their
    default, a vector 40 points evenly spaced over the range of the
    input \code{x} and \code{y}.

    If \code{output="points"}: vectors of \eqn{x}- and \eqn{y}-coordinates of
    output points as given by \code{xo} and \code{yo}.
    }
  \item{z}{
    If \code{output="grid"}:
    matrix of fitted \eqn{z}-values.  The value \code{z[i,j]} is computed
    at the point \eqn{(xo[i], yo[j])}. \code{z} has
    dimensions \code{length(xo)} times \code{length(yo)}.

    If \code{output="points"}: a vector with the calculated z values for
    the output points as given by \code{xo} and \code{yo}.

    If the input was a \code{SpatialPointsDataFrame} a
    \code{SpatialPixelsDataFrame} is returned for \code{output="grid"}
    and a  \code{SpatialPointsDataFrame}  for \code{output="points"}.
  }
}
\references{
  Moebius, A. F. (1827) Der barymetrische Calcul. Verlag v. Johann
  Ambrosius Barth, Leipzig,
  https://books.google.at/books?id=eFPluv_UqFEC&hl=de&pg=PR1#v=onepage&q&f=false

  Franke, R., (1979). A critical comparison of some methods for interpolation
  of scattered data. Tech. Rep. NPS-53-79-003, Dept. of Mathematics, Naval
  Postgraduate School, Monterey, Calif.


  Akima, H. (1978). A Method of Bivariate Interpolation and
  Smooth Surface Fitting for Irregularly Distributed Data Points.
  ACM Transactions on Mathematical Software \bold{4}, 148-164.

  Akima, H. (1996). Algorithm 761: scattered-data surface fitting that has
  the accuracy of a cubic polynomial.
  ACM Transactions on Mathematical Software \bold{22}, 362--371.

}
\author{
  Albrecht Gebhardt <albrecht.gebhardt@aau.at>,
  Roger Bivand <roger.bivand@nhh.no>
}
\note{
Please note that this function tries to be a replacement for the interp() function 
from the akima package. So it should be call compatible for most applications. 
It also offers additional tuning parameters, usually the default settings will fit.
Please be aware that these additional parameters may change in the future as they are
still under development.
}

\seealso{
  \code{\link{interpp}}
}
\examples{
### Use all datasets from Franke, 1979:
data(franke)
## x-y irregular grid points:
oldseed <- set.seed(42)
ni <- 64
xi <- runif(ni,0,1)
yi <- runif(ni,0,1)
xyi <- cbind(xi,yi)
## linear interpolation
fi <- franke.fn(xi,yi,1)
IL <- interp(xi,yi,fi,nx=80,ny=80,method="linear")
## prepare breaks and colors that match for image and contour:
breaks <- pretty(seq(min(IL$z,na.rm=TRUE),max(IL$z,na.rm=TRUE),length=11))
db <- breaks[2]-breaks[1]
nb <- length(breaks)
breaks <- c(breaks[1]-db,breaks,breaks[nb]+db)
colors <- terrain.colors(length(breaks)-1)
image(IL,breaks=breaks,col=colors,main="Franke function 1",
      sub=paste("linear interpolation, ", ni,"points"))
contour(IL,add=TRUE,levels=breaks)
points(xi,yi)
## spline interpolation
fi <- franke.fn(xi,yi,1)
IS <- interp(xi,yi,fi,method="akima",
             kernel="gaussian",solver="QR")
## prepare breaks and colors that match for image and contour:
breaks <- pretty(seq(min(IS$z,na.rm=TRUE),max(IS$z,na.rm=TRUE),length=11))
db <- breaks[2]-breaks[1]
nb <- length(breaks)
breaks <- c(breaks[1]-db,breaks,breaks[nb]+db)
colors <- terrain.colors(length(breaks)-1)
image(IS,breaks=breaks,col=colors,main="Franke function 1",
      sub=paste("spline interpolation, ", ni,"points"))
contour(IS,add=TRUE,levels=breaks)
        points(xi,yi)
## regular grid:
nx <- 8; ny <- 8
xg<-seq(0,1,length=nx)
yg<-seq(0,1,length=ny)
xx <- t(matrix(rep(xg,ny),nx,ny))
yy <- matrix(rep(yg,nx),ny,nx)
xyg<-expand.grid(xg,yg)
## linear interpolation
fg <- outer(xg,yg,function(x,y)franke.fn(x,y,1))
IL <- interp(xg,yg,fg,input="grid",method="linear")
## prepare breaks and colors that match for image and contour:
breaks <- pretty(seq(min(IL$z,na.rm=TRUE),max(IL$z,na.rm=TRUE),length=11))
db <- breaks[2]-breaks[1]
nb <- length(breaks)
breaks <- c(breaks[1]-db,breaks,breaks[nb]+db)
colors <- terrain.colors(length(breaks)-1)
image(IL,breaks=breaks,col=colors,main="Franke function 1",
      sub=paste("linear interpolation, ", nx,"x",ny,"points"))
contour(IL,add=TRUE,levels=breaks)
points(xx,yy)
## spline interpolation
fg <- outer(xg,yg,function(x,y)franke.fn(x,y,1))
IS <- interp(xg,yg,fg,input="grid",method="akima",
             kernel="gaussian",solver="QR")
## prepare breaks and colors that match for image and contour:
breaks <- pretty(seq(min(IS$z,na.rm=TRUE),max(IS$z,na.rm=TRUE),length=11))
db <- breaks[2]-breaks[1]
nb <- length(breaks)
breaks <- c(breaks[1]-db,breaks,breaks[nb]+db)
colors <- terrain.colors(length(breaks)-1)
image(IS,breaks=breaks,col=colors,main="Franke function 1",
      sub=paste("spline interpolation, ", nx,"x",ny,"points"))
contour(IS,add=TRUE,levels=breaks)
        points(xx,yy)
set.seed(oldseed)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dplot }
\keyword{ math }
