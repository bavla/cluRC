# Clustering with relational constraint based on dissimilarity matrix

listRel <- function(R){
  for(i in seq_along(R)) cat(names(R)[i],i,":",R[[i]],"\n")
}

relConM <- function(D,method="max",strategy="tolerant"){
  orDendro <- function(i){if(i<0) return(-i)
    return(c(orDendro(m[i,1]),orDendro(m[i,2])))}
  if(strategy %in% c("tolerant", "leader", "strict")){
    tol <- strategy=="tolerant"; nos <- strategy!="strict"} else 
    {cat("*** Error - Unknown strategy:", strategy,"\n"); return(NULL)}
  cat("Clustering with relational constraint based on dissimilarity matrix\n")
  cat("by Vladimir Batagelj, March 2018\n")
  cat("Method:",method,"  Strategy:",strategy,"\n")
  if(class(D)!="matrix"){cat("*** Error - D should be a matrix\n"); return(NULL)}
  print(paste("Started:",Sys.time()))
# each unit is a cluster; compute inter-cluster dissimilarity matrix
  diag(D) <- Inf; numL <- nrow(D); numLm <- numL-1
  active <- 1:numL; m <- matrix(nrow=numLm,ncol=2)
  node <- rep(0,numL); h <- numeric(numLm); w <- rep(1,numL)
  for(k in 1:numLm){
  # determine the closest pair of clusters (p,q)
    n <- length(active); ind <- rep(Inf,n); dd <- rep(Inf,n)
    for(a in seq_along(active)) {i <- active[a]
      for(j in Ro[[i]]) if(j>0) if(D[i,j] < dd[a]) {dd[a] <- D[i,j]; ind[a] <- j}}
    pq <- which.min(dd)
    if((length(pq)==0)|is.null(pq)) break
    dpq <- dd[pq]
  # join the closest pair of clusters
    p<-active[pq]; q <- ind[pq]; 
    if(is.infinite(q)){
      cat("several components\n")
      dpq <- h[k-1]*1.1; p <- active[1]
      for(q in active[2:length(active)]){
        if(node[p]==0) m[k,1] <- -p else m[k,1] <- node[p]
        if(node[q]==0) m[k,2] <- -q else m[k,2] <- node[q]
        w[p] <- w[q]+w[p]; node[[k]] <- k; h[k] <- dpq
        p <- k; k <- k+1 
      }
      break
    }
    h[k] <- dpq 
    if(node[p]==0) m[k,1] <- -p else m[k,1] <- node[p]
    if(node[q]==0) m[k,2] <- -q else m[k,2] <- node[q]
    active <- setdiff(active,q)
    Rop <- setdiff(Ro[[p]],q); Rip <- setdiff(Ri[[p]],q) 
    Roq <- setdiff(Ro[[q]],p); Riq <- setdiff(Ri[[q]],p)
    for(s in Riq) if(s>0) Ro[[s]] <- setdiff(Ro[[s]],q)
    r <- p; Ror <- 0; Rir <- Rip
    if(tol) Ror <- Rop else
      for(s in Rop) if(s>0) Ri[[s]] <- setdiff(Ri[[s]],p) 
    Ror <- union(Ror,Roq)
    for(s in Roq) if(s>0) Ri[[s]] <- setdiff(union(Ri[[s]],r),q)
    if(nos){
      Rir <- union(Rir,Riq)
      for(s in Riq) if(s>0) Ro[[s]] <- union(Ro[[s]],r)
    }
    Ro[[r]] <- Ror; Ri[[r]] <- Rir; Ro[[q]] <- 0; Ri[[q]] <- 0
  # determine dissimilarities to the new cluster
    for(s in setdiff(active,p)){
      if(method=="max") D[p,s] <- max(D[p,s],D[q,s]) else
      if(method=="min") D[p,s] <- min(D[p,s],D[q,s]) else
      if(method=="ward") { ww <- w[p]+w[q]+w[s]
        D[p,s] <- ((w[q]+w[s])*D[q,s] + (w[p]+w[s])*D[p,s] - w[s]*dpq)/ww
      } else {cat('unknown method','\n'); return(NULL)}
      D[s,p] <- D[p,s]
    }
    w[p] <- w[q]+w[p]; node[[p]] <- k
  }
  hc <- list(merge=m,height=h,order=orDendro(numLm),labels=rownames(D),
    method="relConM",call=NULL,dist.method=method,leaders=NULL)
  class(hc) <- "hclust"
  print(paste("Finished:",Sys.time()))
  return(hc)
}

# a <- scan("C:/Users/batagelj/work/Delphi/Cluse/Cluse/data/RelCon/SomeTy.dis")
# s <- length(a); n <- round((-1+sqrt(1+8*s))/2); nm <- n-1
# D <- matrix(0, nrow=n, ncol=n); D[lower.tri(D,diag=TRUE)] <- a; D <- D+t(D)
# netRel <- "C:/Users/batagelj/work/Delphi/Cluse/Cluse/data/RelCon/SomeTyXY.net"
# R <- read_Pajek_net(netRel,3)
# rownames(D) <- colnames(D) <- names(R)
# for(i in seq_along(R)) if(is.null(R[[i]])) R[[i]] <- 0
# Ro <- R; Ri <- vector("list",length(R))
# names(Ri)<-names(R)
# for(i in 1:length(R)) for(j in R[[i]]) if(j>0) Ri[[j]] <- union(Ri[[j]],i)
# for(i in seq_along(Ri)) if(is.null(Ri[[i]])) Ri[[i]] <- 0
# sRi <- Ri; sRo <- Ro; sd <- d; sD <- D

# Ri <- sRi; Ro <- sRo; d <- sd; D <- sD
# res <- relCon(D,strategy="leader")
# plot(res,hang=-1,main="leader / maximum")
