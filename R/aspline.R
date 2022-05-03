aspline <- function(x, y=NULL, xout, n = 50, ties=mean, method="improved", degree=3) {

    if (! method %in% c("original", "improved")) stop(paste("unknown method:", method))

    x <- xy.coords(x, y) # -> (x,y) numeric of same length
    y <- x$y
    x <- x$x
    nx <- length(x)
    if(any(na <- is.na(x) | is.na(y))) {
	ok <- !na
	x <- x[ok]
	y <- y[ok]
	nx <- length(x)
    }
    if (!identical(ties, "ordered")) {
	if (length(ux <- unique(x)) < nx) {
	    if (missing(ties))
		warning("Collapsing to unique x values")
	    y <- as.vector(tapply(y,x,ties))# as.v: drop dim & dimn.
	    x <- sort(ux)
	    nx <- length(x)
	} else {
	    o <- order(x)
	    x <- x[o]
	    y <- y[o]
	}
    }
    if (nx < 1)
	stop("need at least one non-NA value to interpolate (single value yields constant result)!")
    if (missing(xout)) {
	if (n <= 0)
	    stop("aspline requires n >= 1")
	xout <- seq(x[1], x[nx], length = n)
    }
    aSpline(x,y,xout,method,degree)
}
