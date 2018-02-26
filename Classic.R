setwd("C:/Users/batagelj/work/Delphi/Cluse/Cluse/data/RelCon")

hClust <- function(D,method="max"){
  orDendro <- function(i){if(i<0) return(-i)
    return(c(orDendro(m[i,1]),orDendro(m[i,2])))}
  numL <- nrow(D); numLm <- numL-1
# each unit is a cluster; compute inter-cluster dissimilarity matrix
  diag(D) <- Inf
  print(D); flush.console()
  N <- rownames(D)
  active <- 1:numL; m <- matrix(nrow=numLm,ncol=2)
  node <- rep(0,numL); h <- numeric(numLm); w <- rep(1,numL)
  for(k in 1:numLm){
  # determine the closest pair of clusters (p,q)
    ind <- active[sapply(active,function(i) which.min(D[i,active]))]
    dd <- sapply(active,function(i) min(D[i,active]))
    print(active); print(ind); print(dd)
    pq <- which.min(dd)
  # join the closest pair of clusters
    p<-active[pq]; q <- ind[pq]; h[k] <- D[p,q]
    cat(pq,'join ',p,'/',N[p],' and ',q,'/',N[q],' at level ',D[p,q],'\n')
    if(node[p]==0) m[k,1] <- -p else m[k,1] <- node[p]
    if(node[q]==0) m[k,2] <- -q else m[k,2] <- node[q]
    active <- setdiff(active,p)
  # determine dissimilarities to the new cluster
    for(s in setdiff(active,q)){
      if(method=="max") D[q,s] <- max(D[q,s],D[p,s]) else
      if(method=="min") D[q,s] <- min(D[q,s],D[p,s]) else
      if(method=="ward") { ww <- w[p]+w[q]+w[s]
        D[q,s] <- ((w[q]+w[s])*D[q,s] + (w[p]+w[s])*D[p,s] - w[s]*h[k])/ww
      } else {cat('unknown method','\n'); return(NULL)}
      D[s,q] <- D[q,s]
    }
    w[q] <- w[q]+w[p]; node[[q]] <- k
  }
  hc <- list(merge=m,height=h,order=orDendro(numLm),labels=rownames(D),
    method=method,call=NULL,dist.method=NULL,leaders=NULL)
  class(hc) <- "hclust"
  return(hc)
}


setwd("C:/Users/batagelj/work/clamix/relC/someTy")
a <- scan("SomeTy.dis"); n <- round((sqrt(8*length(a)+1)-1)/2)
D <- matrix(0, nrow=n, ncol=n)
D[lower.tri(D,diag=TRUE)] <- a
D <- D+t(D); diag(D) <- Inf

# h <- hClust(D,method="ward")
h <- hClust(D,method="max")
plot(h,hang=-1)

