\name{triSht2tri}
\alias{triSht2tri}
\title{
Converter to tripack objects
}
\description{
This function converts \code{triSht} objects (from this package) to \code{tri} objects (from tripack package).
}
\usage{
triSht2tri(t.triSht)
}
\arguments{
  \item{t.triSht}{
a class \code{triSht} object as returned by \code{tri.mesh} 
}
}
\note{
The converted objects are not fully compatible with \code{tripack} functions. Basic stuff (printing, plotting) works, \code{tripack::triangles} e.g. does not work.
Voronoi functions from package \code{tripack} are working correctly with translated objects.
}
\value{
	A class \code{tri} object, see tripack package.
}
\author{
A. Gebhardt
}
\seealso{
\code{\link{tri.mesh}},
\code{\link{triSht}}
}
