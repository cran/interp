cells<-function(voronoi.obj)
{
  if(!inherits(voronoi.obj,"voronoi"))
    stop("voronoi.obj must be of class \"voronoi\"")

  tri <- voronoi.obj$tri
  
  nnabs <- integer(tri$n)
  nptr <- integer(tri$n)
  nptr1 <- integer(tri$n)
  nbnos <- integer(tri$n)

  ret <- NULL
  for(i in 1:tri$n){
    vs <- voronoi.findvertices(i, voronoi.obj)
    if(length(vs)>0){
      center <- c(tri$x[i],tri$y[i])
      neighbours <- sort(c(arcs(tri)[arcs(tri)[,1]==i,2],arcs(tri)[arcs(tri)[,2]==i,1]))
      nodes <- rbind(voronoi.obj$x[vs],voronoi.obj$y[vs])
      rownames(nodes) <- c("x","y")
      area <- voronoi.polyarea( voronoi.obj$x[vs], voronoi.obj$y[vs])
      ret[[i]] <- list(cell=i,center=center,
                       neighbours=neighbours,
                       nodes=nodes,area=area)
    } else {
      center <- c(tri$x[i],tri$y[i])
      neighbours <- sort(c(arcs(tri)[arcs(tri)[,1]==i,2],arcs(tri)[arcs(tri)[,2]==i,1]))
      nodes <- NA # should better return at least the non-dummy nodes
      area <- NA
      ret[[i]] <- list(cell=i,center=center,
                       neighbours=neighbours,
                       nodes=nodes,area=area)
    }
  }
  
  ret
}

