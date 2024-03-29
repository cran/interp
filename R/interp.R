interpp <- function(x, y=NULL, z, xo, yo=NULL, linear = TRUE, extrap = FALSE,
                    duplicate = "error", dupfun = NULL,
                    deltri = "shull"){
    interp(x,y,z,xo,yo,linear,extrap,duplicate,dupfun,deltri,
           input="points", output="points")
}


interp <- function(x, y=NULL, z,
                   xo = seq(min(x), max(x), length = nx),
                   yo = seq(min(y), max(y), length = ny),
                   linear = (method=="linear"), extrap = FALSE,
                   duplicate = "error", dupfun = NULL,
                   nx = 40, ny = 40, input = "points", output = "grid",
                   method="linear", deltri="shull",
                   h=0, kernel="gaussian", solver="QR", degree=3,
                   baryweight=TRUE,
                   autodegree=FALSE,
                   adtol=0.1,
                   smoothpde=FALSE,
                   akimaweight=TRUE,
                   nweight=25,
                   na.rm=FALSE)
    {
        if(method=="linear")
            linear <- TRUE
        ## handle sp data, save coordinate and value names
        is.sp <- FALSE
        sp.coord <- NULL
        sp.z <- NULL
        sp.proj4string <- NULL
        if(is.null(y)&&is.character(z)){
            if((inherits(x,"SpatialPointsDataFrame") ||
                inherits(x,"SpatialPixelsDataFrame")) &&
                requireNamespace("sp", quietly=TRUE)) {
                sp.coord <- dimnames(sp::coordinates(x))[[2]]
                sp.z <- z
                sp.proj4string <- x@proj4string
                z <- x@data[,z]
                y <- sp::coordinates(x)[,2]
                x <- sp::coordinates(x)[,1]
                if(na.rm){
                    z <- z[!is.na(z)]
                    x <- z[!is.na(z)]
                    y <- z[!is.na(z)]
                }
                is.sp <- TRUE
                xo = seq(min(x), max(x), length = nx)
                yo = seq(min(y), max(y), length = ny)
            } else
                stop("either x,y,z are numerical or x is SpatialPointsDataFrame / SpatialPixelsDataFrame and z a name of a data column in x")
        }
        if(!(all(is.finite(x)) && all(is.finite(y)) && all(is.finite(z))))
            stop("missing values and Infs not allowed")
        drx <- diff(range(x))
        dry <- diff(range(y))
        if(drx == 0 || dry == 0)
            stop("all data collinear")    # other cases caught in Fortran code
        if(drx/dry > 10000 || drx/dry < 0.0001)
            stop("scales of x and y are too dissimilar")
        n <- length(x)
        nx <- length(xo)
        ny <- length(yo)
        if(input=="points" && (length(y) != n || length(z) != n))
            stop("Lengths of x, y, and z do not match")
        if(input=="points"){
           dups_found <- isTRUE(anyDuplicated(cbind(x, y), MARGIN=1) != 0L)
           if (dups_found) {
               if(duplicate == "error") {
                   stop("duplicate data points: need to set 'duplicate = ..' ")
               }   else { ## duplicate != "error"
        
                   xy <- paste(x, y, sep = ",") # trick for 'duplicated' (x,y)-pairs
                   i <- match(xy, xy)
                   if(duplicate == "user")
                       dupfun <- match.fun(dupfun)#> error if it fails
        
                   ord <- !duplicated(xy)
                   if(duplicate != "strip") {
                       centre <- function(x)
                           switch(duplicate,
                                  mean = mean(x),
                                  median = median(x),
                                  user = dupfun(x))
                       z <- unlist(lapply(split(z,i), centre))
                   } else {
                       z <- z[ord]
                   }
                   x <- x[ord]
                   y <- y[ord]
                   n <- length(x)
               }
           }
        }
        if(input=="grid"){
          if(!is.sp){
            nxi <- length(x)
            nyi <- length(y)
            
            x<- c(matrix(rep(x,nyi),nrow=nxi,ncol=nyi,byrow=FALSE))
            y<- c(matrix(rep(y,nxi),nrow=nxi,ncol=nyi,byrow=TRUE))
            z<- c(z)
          } 
        }
        if(method=="linear"||method=="akima"||method=="bilinear"){

            if(deltri=="deldir"){
                if(!linear)
                    stop("method=\"akima\" (linear=FALSE) is not implemented for deltri=\"deldir\"!")
                triangles <- triang.list(deldir(x=x,y=y,z=z))
                ans <- interpDeltri(xo,yo,z,triangles,input,output)
            } else if(deltri=="shull"){
                if(method=="bilinear")
                    if(input=="grid")
                        ans <- bilinear(xo,yo,x,y,z)
                    else
                        stop("method=\"bilinear\" needs input=\"grid\"!")
                if(method=="linear"||method=="akima")    
                    ans <- interpShull(xo,yo,x,y,z,linear,input,output,
                                       kernel,h,
                                       solver,degree,
                                       baryweight,
                                       autodegree,adtol,
                                       smoothpde,akimaweight,nweight)
                if(ans$err==-13){
                  ## retry with jitter
                  ans <- interpShull(xo,yo,jitter(x,1e-3),jitter(y,1e-3),z,linear,input,output,
                                     kernel,h,
                                     solver,degree,
                                     baryweight,
                                     autodegree,adtol,
                                     smoothpde,akimaweight,nweight)
                } else if(ans$err==-14){
                  ## retry with jitter (or TODO: rescale)
                  ans <- interpShull(xo,yo,jitter(x,1e-3),jitter(y,1e-3),z,linear,input,output,
                                     kernel,h,
                                     solver,degree,
                                     baryweight,
                                     autodegree,adtol,
                                     smoothpde,akimaweight,nweight)
                }
                
                if(output=="points") # back to vector from matrix:
                    ans$z <- c(ans$z)
            } else
                stop(paste("unknown triangulation method", deltri))
        } else
            stop(paste("method=\"",method,"\" not implemented!",sep=""))
        ## prepare return value
        if (is.sp && requireNamespace("sp", quietly=TRUE)) {
            zm <- nx
            zn <- ny
            zvec <- c(ans$z)
            xvec <- c(matrix(rep(ans$x,zn),nrow=zm,ncol=zn,byrow=FALSE))
            yvec <- c(matrix(rep(ans$y,zm),nrow=zm,ncol=zn,byrow=TRUE))
            nona <- !is.na(zvec)
            ret <- data.frame(xvec[nona],yvec[nona],zvec[nona])
            names(ret) <- c(sp.coord[1],sp.coord[2],sp.z)
            sp::coordinates(ret) <- sp.coord
            ret@proj4string <- sp.proj4string
            sp::gridded(ret) <- TRUE
        } else {
            if(output=="grid")
                ret <- list(x=ans$x,y=ans$y,z=matrix(ans$z,nx,ny))
            else
                ret <- list(x=ans$x,y=ans$y,z=ans$z)
        }
        ret
    }
