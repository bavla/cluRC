# Clustering with relational constraint based on a dictionary

key <- function(i,j) 
  if(i<j) return(as.character(i*np+j)) else return(as.character(j*np+i))

listRel <- function(R){
  for(i in seq_along(R)) cat(names(R)[i],i,":",R[[i]],"\n")
}

relConH <- function(method="max",strategy="tolerant"){
  orDendro <- function(i){if(i<0) return(-i)
    return(c(orDendro(m[i,1]),orDendro(m[i,2])))}
  idx <- function(i,j) {if (i<j) return(n*(i-1) - i*(i-1)/2 + j-i) else
    return(n*(j-1) - j*(j-1)/2 + i-j)}
  if(strategy %in% c("tolerant", "leader", "strict")){
    tol <- strategy=="tolerant"; nos <- strategy!="strict"} else 
    {cat("*** Error - Unknown strategy:", strategy,"\n"); return(NULL)}
  cat("Clustering with relational constraint based on a dictionary\n")
  cat("by Vladimir Batagelj, March 2018\n")
  cat("Method:",method,"  Strategy:",strategy,"\n")
  print(paste("Started:",Sys.time()))
# each unit is a cluster; compute inter-cluster dissimilarity matrix
  numL <- attr(hD,"Size"); numLm <- numL-1
  active <- 1:numL; m <- matrix(nrow=numLm,ncol=2)
  node <- rep(0,numL); h <- numeric(numLm); w <- rep(1,numL)
  for(k in 1:numLm){
  # determine the closest pair of clusters (p,q)
    nn <- length(active); ind <- rep(Inf,nn); dd <- rep(Inf,nn)
    for(a in seq_along(active)) {i <- active[a]
    #  if((length(Ro[[i]])==1)&&(Ro[[i]][1]==0)) dd[a] <- Inf else
      for(j in Ro[[i]]) if((j>0)&&(i!=j)) { dij <- get(key(i,j),envir=hD)
        if(dij < dd[[a]]) {dd[[a]] <- dij; ind[[a]] <- j}}}
    pq <- which.min(dd)
    if((length(pq)==0)|is.null(pq)) break
    dpq <- dd[[pq]]
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
    Sq <- setdiff(union(Roq,Riq),0)
    S <- setdiff(union(union(Rop,Rip),Sq),0)
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
    for(s in S){
      if(method=="max") assign(key(p,s),max(
         get0(key(p,s),envir=hD,inherits=FALSE,ifnotfound=0),
         get0(key(q,s),envir=hD,inherits=FALSE,ifnotfound=0)),envir=hD) else
      if(method=="min") assign(key(p,s),min(
         get0(key(p,s),envir=hD,inherits=FALSE,ifnotfound=Inf),
         get0(key(q,s),envir=hD,inherits=FALSE,ifnotfound=Inf)),envir=hD) else
      if(method=="ward") { ww <- w[p]+w[q]+w[s]
        D[idx(p,s)] <- ((w[q]+w[s])*D[idx(q,s)] + (w[p]+w[s])*D[idx(p,s)] - w[s]*dpq)/ww
      } else {cat('unknown method','\n'); return(NULL)}
    }
    for(s in union(Sq,p)) {ky <- key(q,s) 
      if(exists(ky,envir=hD,inherits=FALSE)) remove(list=ky,envir=hD,inherits=FALSE)}
    w[p] <- w[q]+w[p]; node[[p]] <- k
  }
  hc <- list(merge=m,height=h,order=orDendro(numLm),labels=attr(hD,"Labels"),
    method="cluRelH",call=NULL,dist.method=method,leaders=NULL)
  class(hc) <- "hclust"
  print(paste("Finished:",Sys.time()))
  return(hc)
}

# # Create initial weights - dissimilarities
# # n is number of nodes
# Ri <- sRi; Ro <- sRo; D <- sD
# n <- nrow(D); np <- n+1; hD <- new.env()
# for(i in 1:length(R)) for(j in R[[i]])  assign(key(i,j),D[i,j],envir=hD)
# attr(hD,"Size") <- n; attr(hD,"Labels") <- names(R)
# res <- relConH(strategy="tolerant")
# plot(res,hang=-1,main="tolerant / maximum")

