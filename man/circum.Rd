\name{circum}
\Rdversion{1.1}
\alias{circum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Determine the circumcircle (and some other characteristics) of a triangle
}
\description{
This function returns the circumcircle of a triangle and some additonal
values used to determine them.
}
\usage{
circum(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of three elements, giving the x coordinatres of the
    triangle nodes.
}
  \item{y}{Vector of three elements, giving the y coordinatres of the
    triangle nodes.
}
}
\details{
  This is an interface to the Fortran function CIRCUM found in TRIPACK.
}
\value{
  \item{ x }{ 'x' coordinate of center}
  \item{ y }{ 'y' coordinate of center}
  \item{ radius }{ circumcircle radius }
  \item{ signed.area }{ signed area of riangle (positive iff nodes are
    numbered counter clock wise) }
  \item{ aspect.ratio }{ ratio "radius of inscribed circle"/"radius of
    circumcircle", varies between 0 and 0.5

    0 means collinear points, 0.5 equilateral trangle.
    }
}
\references{
https://math.fandom.com/wiki/Circumscribed_circle#Coordinates_of_circumcenter, visited march 2022.
}
\author{
  A. Gebhardt
}
\note{
  This function is mainly intended to be used by \code{\link{circumcircle}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{circumcircle}}
}
\examples{
circum(c(0,1,0),c(0,0,1))

tr <- list()
tr$t1 <-list(x=c(0,1,0),y=c(0,0,1))
tr$t2 <-list(x=c(0.5,0.9,0.7),y=c(0.2,0.9,1))
tr$t3 <-list(x=c(0.05,0,0.3),y=c(0.2,0.7,0.1))
plot(0,0,type="n",xlim=c(-0.5,1.5),ylim=c(-0.5,1.5))
for(i in 1:3){
    x <- tr[[i]]$x
    y <- tr[[i]]$y
    points(x,y,pch=c("1","2","3"),xlim=c(-0.5,1.5),ylim=c(-0.5,1.5))
    cc =circum(x,y)
    lines(c(x,x[1]),c(y,y[1]))
    points(cc$x,cc$y)
    if(cc$signed.area<0)
      circles(cc$x,cc$y,cc$radius,col="blue",lty="dotted")
    else
      circles(cc$x,cc$y,cc$radius,col="red",lty="dotted")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }

