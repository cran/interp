\name{convex.hull}
\title{Return the convex hull of a triangulation object}
\usage{convex.hull(tri.obj, plot.it=FALSE, add=FALSE,...)
       ConvexHull(x,y)
}
\alias{convex.hull}
\alias{ConvexHull}
\arguments{
  \item{tri.obj}{object of class \code{\link{triSht}}}
  \item{plot.it}{logical, if \code{TRUE} the convex hull of
    \code{tri.obj} will be plotted.}
  \item{add}{logical. if \code{TRUE} (and \code{plot.it=TRUE}), add to
    a current plot.}
  \item{...}{additional plot arguments}
  \item{x}{only for \code{ConvexHull()}: \code{x} coordinates for C++ call to \code{ConvexHull}}
  \item{y}{only for \code{ConvexHull()}: see \code{x}}
}
\description{
  Given a triangulation \code{tri.obj} of \eqn{n} points in the plane, this
  subroutine returns two vectors containing the coordinates
  of the nodes on the boundary of the convex hull.
  
  \code{ConvexHull} is an experimental C++ implementation of Grahams Scan
  without previous triangulation, should be much faster.
}
\note{
  In case that there are several collinear nodes on the convex hull
  \code{convex.hull} will return them all while \code{ConvexHull} will only
  give edge points.
}
\value{
  \item{x}{x coordinates of boundary nodes.}
  \item{y}{y coordinates of boundary nodes.}
}
\author{
  Albrecht Gebhardt <albrecht.gebhardt@aau.at>,
  Roger Bivand <roger.bivand@nhh.no>
}
\seealso{
  \code{\link{triSht}}, \code{\link{print.triSht}},
  \code{\link{plot.triSht}}, \code{\link{summary.triSht}}, \code{\link{triangles}}.
}
\examples{
## random points:
rand.tr<-tri.mesh(runif(10),runif(10))
plot(rand.tr)
rand.ch<-convex.hull(rand.tr, plot.it=TRUE, add=TRUE, col="red")
## use a part of the quakes data set:
data(quakes)
quakes.part<-quakes[(quakes[,1]<=-17 & quakes[,1]>=-19.0 &
                     quakes[,2]<=182.0 & quakes[,2]>=180.0),]
quakes.tri<-tri.mesh(quakes.part$lon, quakes.part$lat, duplicate="remove")
plot(quakes.tri)
convex.hull(quakes.tri, plot.it=TRUE, add=TRUE, col="red")
}
\keyword{spatial}
