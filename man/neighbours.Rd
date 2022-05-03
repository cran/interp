\name{neighbours}
\title{List of neighbours from a triangulation or voronoi object}
\author{A. Gebhardt}
\usage{neighbours(obj)
}
\alias{neighbours}
\arguments{
  \item{obj}{object of class \code{"triSht"} or \code{"voronoi.mosaic"}}
}
\description{Extract a list of neighbours from a triangulation or
  voronoi object
}
\value{nested list of neighbours per point
}
\seealso{\code{\link{triSht}}, \code{\link{print.triSht}}, \code{\link{plot.triSht}}, \code{\link{summary.triSht}}, \code{\link{triangles}}
}
\examples{
data(tritest)
tritest.tr<-tri.mesh(tritest$x,tritest$y)
tritest.nb<-neighbours(tritest.tr)
}
\keyword{spatial}
