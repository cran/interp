"voronoi.mosaic" <- function(x,y=NULL,duplicate="error")
  {

    dummy.node<-function(x0,y0,x1,y1,x2,y2,d)
      {
        # determine a direction orthogonal to p1--p2
        #
        #              p_1
        #               |
        #               |d
        #    p_0 ------>+ - - - - -> dummy_node
        #          r    |
        #               V
        #              p_2------->
        #                    n
        # two versions, r and n
        #

        dx<-  x2-x1
        dy<-  y2-y1
        nx<- -dy
        ny<-  dx

        rx<-(x1+x2)/2-x0
        ry<-(y1+y2)/2-y0

        lr<-sqrt(rx^2+ry^2)
        ln<-sqrt(nx^2+ny^2)
        # choose the numerically better version
        if(lr > ln)
          {
            vx<-rx/lr
            vy<-ry/lr

            if(in.convex.hull(ret$tri,x0,y0))
              d <- d
            else
              d <- -d
          }
        else
          {
            vx<-nx/ln
            vy<-ny/ln
            eps<-1e-7
            if(in.convex.hull(ret$tri,(x1+x2)/2+eps*vx,(y1+y2)/2+eps*vy))
              d <- - d
            else
              d <- d
          }
        list(x=x0+d*vx,y=y0+d*vy)
      }

    if(inherits(x,"tri")){
        if(!is.null(x$tlist))
            stop("this \"tri\" object has been created with tripack::tri.mesh,\n  recreate it with interp::tri.mesh!\n  The $call element gives a hint how it was created.")
        tri.obj <- x
    } else
        tri.obj<-tri.mesh(x=x,y=y,duplicate=duplicate)

    nt<-tri.obj$nt

    ret<-list(x=tri.obj$cclist[,"x"],
              y=tri.obj$cclist[,"y"],
              node=(tri.obj$cclist[,"area"]>0),
              area=tri.obj$cclist[,"area"],
              ratio=tri.obj$cclist[,"ratio"],
              radius=tri.obj$cclist[,"r"],
              n1=tri.obj$trlist[,"j1"],
              n2=tri.obj$trlist[,"j2"],
              n3=tri.obj$trlist[,"j3"],
              p1=tri.obj$trlist[,"i1"],
              p2=tri.obj$trlist[,"i2"],
              p3=tri.obj$trlist[,"i3"],
              tri=tri.obj)


    ret$dummy.x<-integer(0)
    ret$dummy.y<-integer(0)
    dummy.cnt<-0
    dmax<-max(diff(range(ret$x)),diff(range(ret$y)))
    n<-length(ret$x)
    # add dummy nodes on the border of the triangulation
    for (i in 1:n)
      {
        if(ret$node[i])
          # Triangle i has positive area.
          {
            # Find neighbour triangles
            tns<-sort(c(ret$n1[i],ret$n2[i],ret$n3[i]))
            ins <- order(c(ret$n1[i],ret$n2[i],ret$n3[i]))
            tn1<-tns[1]
            tn2<-tns[2]
            tn3<-tns[3]
            # Handle special cases on the border:
            # (This should better be done in the FORTRAN code!)
            if(any(tns==0))
              {
                if(tns[2]!=0)
                  {
                    # Only one edge of i coincides with border.
                    # Determine nodes of triangle i
                    tr<-c(ret$p1[i],ret$p2[i],ret$p3[i])
                    # Which of these nodes are border nodes (2)?
                    ns<-tr[on.convex.hull(ret$tri,
                                          ret$tri$x[tr],
                                          ret$tri$y[tr])]
                    if(length(ns)==2)
                      {
                        # 2 points on hull
                        i1<-ns[1]
                        i2<-ns[2]
                        # Find a dummy node
                        pn<-dummy.node(ret$x[i],ret$y[i],
                                       ret$tri$x[i1],ret$tri$y[i1],
                                       ret$tri$x[i2],ret$tri$y[i2],
                                       dmax)
                        dummy.cnt<- dummy.cnt+1
                        ret$dummy.x[dummy.cnt]<-pn$x
                        ret$dummy.y[dummy.cnt]<-pn$y
                        # update neighbour relation
                        # (negative index indicates dummy node)
                        if(ret$n1[i]==0) ret$n1[i]<- -dummy.cnt
                        if(ret$n2[i]==0) ret$n2[i]<- -dummy.cnt
                        if(ret$n3[i]==0) ret$n3[i]<- -dummy.cnt
                      }
                    # Other cases:
                    #   1 point on hull -- should not happen at all
                    #   3 points on hull -- should not happen here
                    #     see "else" tree
                  }
                else
                  {
                    # Two edges of i coincide with border.
                    # (= 3 points on hull )
                    # that means this triangle forms one corner of
                    # the convex hull
                    # Find out which edge of triangle i is not
                    # on the border: (check if midpoints of edges lay
                    # on hull)
                    tr<-c(ret$p1[i],ret$p2[i],ret$p3[i])
                    edge<-list(from=tr[c(1,2,3)],to=tr[c(2,3,1)])
                    mx <- (ret$tri$x[edge$from]+ret$tri$x[edge$to])/2
                    my <- (ret$tri$y[edge$from]+ret$tri$y[edge$to])/2
                    eonb <- on.convex.hull(ret$tri,mx,my)
                    # Find two dummy nodes
                    for (id in 1:3){
                      if (eonb[id]){
                        pn<-dummy.node(ret$x[i],ret$y[i],
                                       ret$tri$x[edge$from[id]],
                                       ret$tri$y[edge$from[id]],
                                       ret$tri$x[edge$to[id]],
                                       ret$tri$y[edge$to[id]],
                                       dmax)
                        dummy.cnt<- dummy.cnt+1
                        ret$dummy.x[dummy.cnt]<-pn$x
                        ret$dummy.y[dummy.cnt]<-pn$y
                        # update neighbour relation
                        # (negative index indicates dummy node)
                        if(ret$n1[i]==0) ret$n1[i]<- -dummy.cnt
                        else
                          if(ret$n2[i]==0) ret$n2[i]<- -dummy.cnt
                          else
                            if(ret$n3[i]==0) ret$n3[i]<- -dummy.cnt
                      }
                    }
                  }
              }
          }
        else
          {
            # A triangle i with area 0:
            # This can't happen on the border (already removed in FORTRAN code!).
            # Do nothing.
            tmp<-0
          }
      }
    ret$call <- match.call()
    class(ret) <- "voronoi"
    ret
  }
