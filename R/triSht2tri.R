## cyclic rotate a vector moving its minimum element first
## assumption: elements are unique
minrot =function(x){
    n=length(x)
    mpos = which.min(x)
    c(1:n,1:n)[mpos:(mpos+n-1)]
}

## operates on trlist structure from tri.mesh
## rotates first 3 columns so that minimum element is in first row
## row 4 to 6 and 7 to 9 are rotatetd accoringly (but later ignored)
## then create all three cyclic permutations of the first three elements
## repeatedly and add those two rows
minsort = function(t){
    ret<-NULL
    for(i in 1:dim(t)[1]){
        mr =minrot(t[i,1:3])
        ret=rbind(ret,t[i,c(mr,mr+3,mr+6)])
        nr = c(mr[2:3],mr[1])
        ret=rbind(ret,t[i,c(nr,nr+3,nr+6)])
        nr = c(nr[2:3],nr[1])
        ret=rbind(ret,t[i,c(nr,nr+3,nr+6)])
    }
    ret
}


triSht2tri = function(t.triSht){
    ms = minsort(t.triSht$trlist)
    ## this will contain a list of adjacency vectors for each
    ## triangulation node:
    msl=NULL
    ## number of data points
    nd=t.triSht$n

    for(i in 1:nd){
        ## get columns 2 and 3 matching "i" in column 1
        tmp=ms[ms[,1]==i,2:3]
        ## this can be a matrix or a vector, handle it appropriately
        if(is.vector(tmp))
            msl[[i]]=c(tmp)
        else{
            ## collect the adjacencies in lst
            lst=NULL
            ## sometimes start end end nodes are unique and only once found,
            ## sometimes twice, if twice remove one row to get
            ## also unique start and end nodes
            ctmp=table(tmp)
            if(any(ctmp==1)){
                ## nothing
            } else {
                ## each time we remove one row from tmp
                ## we have to ensure it stays a matrix
                ## even in the case of only one row, this
                ## is the reason for all the odim... stuff here and below
                odim=dim(tmp)
                tmp=tmp[-1,]
                dim(tmp)=c(odim[1]-1,odim[2])
            }
            ## start/stop nodes
            st=as.integer(names(ctmp)[ctmp==1])
            ## start node
            p1=which(tmp[,1]==st[1] | tmp[,1]==st[2])
            ## stop node
            p2=which(tmp[,2]==st[1] | tmp[,2]==st[2])
            ## take the first pair if adjencent neighbour nodes
            lst=tmp[1,]
            ## remove them from tmp
            odim=dim(tmp)
            tmp=tmp[-1,]
            dim(tmp)=c(odim[1]-1,odim[2])
            llst=length(lst)
            repeat{
                ## until tmp is empty
                if(length(c(tmp))/2==0){
                    break
                }
                ## find a pair of nodes matching left to lst
                lft=which(tmp[,2]==lst[1])
                if(length(lft)==1){
                    ## append it left to lst
                    lst=c(tmp[lft,1],lst)
                    ## remove it from tmp
                    odim=dim(tmp)
                    tmp=tmp[-lft,]
                    dim(tmp)=c(odim[1]-1,odim[2])
                    llst=length(lst)
                }
                ## find a pair of nodes matching right to lst
                rgt=which(tmp[,1]==lst[llst])
                if(length(rgt)==1){
                    ## append it right to lst
                    lst=c(lst,tmp[rgt,2])
                    ## remove it from tmp
                    odim=dim(tmp)
                    tmp=tmp[-rgt,]
                    dim(tmp)=c(odim[1]-1,odim[2])
                    llst=length(lst)
                }
            }
            ## TODO: use start and stop nodes to check
            ## if results are correct.

            ## store the adjacency list for node i
            msl[[i]]=lst
            ## do some cleaning
            attr(msl[[i]],"names") <- NULL
        }
        ## check for counter clock wise ordering of adjacency list,
        ## otherwise correct it.
        if(!interp::left(t.triSht$x[i],t.triSht$y[i],t.triSht$x[msl[[i]][1]],t.triSht$y[msl[[i]][1]],t.triSht$x[msl[[i]][2]],t.triSht$y[msl[[i]][2]]))
            msl[[i]]=msl[[i]][length(msl[[i]]):1]
    }
    ## prepare TRIPACK data structures
    tlist=NULL; tlptr=NULL;tlend=NULL
    for(l in msl){
        last=length(tlist)
        ## concat lists
        tlist=append(tlist,l)
        ## save end points of lists
        tlend=append(tlend,length(tlist))
        ## find pointers from list end points to previous
        ## elements of the list, a little bit tricky
        tlptr=append(c(tlptr,last+length(l)),last+(1:(length(l)-1)))
    }
    ## check if list end points ar on convex hull
    ## if yes mark them with negative sign
    ch = interp::on.convex.hull(t.triSht,t.triSht$x[tlist[tlend]],t.triSht$y[tlist[tlend]])
    tlist[tlend[ch]]=-tlist[tlend[ch]]
    ## collect output
    obj=list(
        n=t.triSht$n,
        x=t.triSht$x,
        y=t.triSht$y,
        tlist=tlist,
        tlptr=tlptr,
        tlend=tlend,
        tlnew=tlend[length(tlend)]+1,
        nc=0,
        lc=0,
        call="converted from interp",
        adjlist=msl
    )
    class(obj)<-"tri"
    obj
}

