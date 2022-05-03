neighbours<-function(obj)
{
  if(!inherits(obj,"triSht") && !inherits(obj,"voronoi.mosaic"))
    stop("obj must be of class \"triSht\" or \"voronoi.mosaic\"")
  if(inherits(obj,"triSht")){
      n <- obj$n
      ret<-rep(NULL,obj$n)
  } else if(inherits(obj,"voronoi.mosiac")) {
      n <- obj$tri$n
      ret<-rep(NULL,obj$tri$n)
  }
  
  for (i in 1:n)
  {
      if(inherits(obj,"triSht"))
          ret[i]<-list( sort(c(arcs(obj)[arcs(obj)[,1]==i,2],arcs(obj)[arcs(obj)[,2]==i,1])))
      else if(inherits(obj,"voronoi.mosiac"))
          ret[i]<-list( sort(c(arcs(obj$tri)[arcs(obj$tri)[,1]==i,2],arcs(obj$tri)[arcs(obj$tri)[,2]==i,1])))
  }
  ret
}
