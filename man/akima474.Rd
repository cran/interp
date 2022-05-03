\name{akima474}
\alias{akima474}
\title{
Sample data from Akima's Bicubic Spline Interpolation code (TOMS 474)
}
\description{
\code{akima474} is a list with vector components \code{x}, \code{y} and a matrix \code{z} which
represents a smooth surface of \code{z} values at the points
of a regular grid spanned by the vectors \code{x} and \code{y}.

}
\references{
Hiroshi Akima, Bivariate Interpolation and
Smooth Surface Fitting Based
on Local Procedures [E2], 
Communications of ACM,
Vol. 17, No. 1, January 1974, pp. 26-30
}
\examples{
\dontrun{
library(rgl)
data(akima474)
# data
rgl.spheres(akima474$x,akima474$z , akima474$y,0.5,color="red")
rgl.bbox()
# bivariate linear interpolation
# interp:
akima474.li <- interp(akima474$x, akima474$y, akima474$z, 
                   xo=seq(min(akima474$x), max(akima474$x), length = 100),
                   yo=seq(min(akima474$y), max(akima474$y), length = 100))
# interp surface:
rgl.surface(akima474.li$x,akima474.li$y,akima474.li$z,color="green",alpha=c(0.5))
# interpp:
akima474.p <- interpp(akima474$x, akima474$y, akima474$z,
                    runif(200,min(akima474$x),max(akima474$x)),
                    runif(200,min(akima474$y),max(akima474$y)))
# interpp points:
rgl.points(akima474.p$x,akima474.p$z , akima474.p$y,size=4,color="yellow")

# bivariate spline interpolation
# data
rgl.spheres(akima474$x,akima474$z , akima474$y,0.5,color="red")
rgl.bbox()
# bivariate cubic spline interpolation
# interp:
akima474.si <- interp(akima474$x, akima474$y, akima474$z, 
                   xo=seq(min(akima474$x), max(akima474$x), length = 100),
                   yo=seq(min(akima474$y), max(akima474$y), length = 100),
                   linear = FALSE, extrap = TRUE)
# interp surface:
rgl.surface(akima474.si$x,akima474.si$y,akima474.si$z,color="green",alpha=c(0.5))
# interpp:
akima474.sp <- interpp(akima474$x, akima474$y, akima474$z,
                    runif(200,min(akima474$x),max(akima474$x)),
                    runif(200,min(akima474$y),max(akima474$y)),
                   linear = FALSE, extrap = TRUE)
# interpp points:
rgl.points(akima474.sp$x,akima474.sp$z , akima474.sp$y,size=4,color="yellow")

}
}
\keyword{datasets}
% Converted by Sd2Rd version 0.2-a3.
