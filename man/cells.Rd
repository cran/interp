\name{cells}
\alias{cells}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ extract info about voronoi cells }
\description{
  This function returns some info about the cells of a voronoi mosaic,
  including the coordinates of the vertices and the cell area.
}
\usage{
cells(voronoi.obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{voronoi.obj}{ object of class \code{voronoi}  }
}
\details{
  The function calculates the neighbourhood relations between the
  underlying triangulation and translates it into the neighbourhood
  relations between the voronoi cells.
}
\value{
  retruns a list of lists, one entry for each voronoi cell which contains
  \item{cell }{cell index}
  \item{center }{cell 'center'}
  \item{neighbours }{neighbour cell indices}
  \item{nodes}{2 times \code{nnb} matrix with vertice coordinates}
  \item{area}{cell area}
}

\author{ A. Gebhardt }
\note{ outer cells have \code{area=NA}, currently also \code{nodes=NA}
  which is not really useful -- to be done later }


\seealso{  \code{\link{voronoi.mosaic}},  \code{\link{voronoi.area}} }
\examples{
data(tritest)
tritest.vm <- voronoi.mosaic(tritest$x,tritest$y)
tritest.cells <- cells(tritest.vm)
# higlight cell 12:
plot(tritest.vm)
polygon(t(tritest.cells[[12]]$nodes),col="green")
# put cell area into cell center:
text(tritest.cells[[12]]$center[1],
     tritest.cells[[12]]$center[2],
     tritest.cells[[12]]$area)
}
\keyword{ spatial }% at least one, from doc/KEYWORDS
