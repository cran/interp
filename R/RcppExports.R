# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

interpDeltri <- function(x, y, zD, t, input = "points", output = "grid") {
    .Call('interp_interpDeltri', PACKAGE = 'interp', x, y, zD, t, input, output)
}

interpShull <- function(x, y, xD, yD, zD, linear = TRUE, input = "points", output = "grid") {
    .Call('interp_interpShull', PACKAGE = 'interp', x, y, xD, yD, zD, linear, input, output)
}

locpoly.partderiv.grid <- function(x, y, xD, yD, zD, kernel = "gaussian", h = as.numeric( c(0.25,0.25)), solver = "QR", degree = 3L, smoothpde = FALSE, akimaweight = FALSE, nweight = 25L) {
    .Call('interp_partDerivGrid', PACKAGE = 'interp', x, y, xD, yD, zD, kernel, h, solver, degree, smoothpde, akimaweight, nweight)
}

locpoly.partderiv.points <- function(x, y, xD, yD, zD, kernel = "gaussian", h = as.numeric( c(0.25,0.25)), solver = "QR", degree = 3L, smoothpde = FALSE, akimaweight = FALSE, nweight = 25L) {
    .Call('interp_partDerivPoints', PACKAGE = 'interp', x, y, xD, yD, zD, kernel, h, solver, degree, smoothpde, akimaweight, nweight)
}

nearest.neighbours <- function(x, y) {
    .Call('interp_nearestNeighbours', PACKAGE = 'interp', x, y)
}

shull.deltri <- function(x, y) {
    .Call('interp_shullDeltri', PACKAGE = 'interp', x, y)
}

triFind <- function(nT, xT, yT, i1, i2, i3, x, y) {
    .Call('interp_triFind', PACKAGE = 'interp', nT, xT, yT, i1, i2, i3, x, y)
}

left <- function(x1, y1, x2, y2, x0, y0, eps = 1E-16) {
    .Call('interp_left', PACKAGE = 'interp', x1, y1, x2, y2, x0, y0, eps)
}

on <- function(x1, y1, x2, y2, x0, y0, eps = 1E-16) {
    .Call('interp_on', PACKAGE = 'interp', x1, y1, x2, y2, x0, y0, eps)
}

inHull <- function(triObj, x, y, eps = 1E-16) {
    .Call('interp_inHull', PACKAGE = 'interp', triObj, x, y, eps)
}

onHull <- function(triObj, x, y, eps = 1E-16) {
    .Call('interp_onHull', PACKAGE = 'interp', triObj, x, y, eps)
}

