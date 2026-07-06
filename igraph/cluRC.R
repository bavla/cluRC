# cluRC - clustering with relational constraint
# March 2018 and July 2026
# by Vladimir Batagelj
#
# library(igraph); library(sf); library(tmap); library(spdep); library(pals)

z <- function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)


varCutree <- function(R,var,vmin,vmax){
  mark <- function(t,c){
    if(t<0) part[-t] <<- c else {mark(R$merge[t,1],c); mark(R$merge[t,2],c)}
  }
  n <- length(var); part <- rep(999999,n); nclust <- 0
  value <- cbind(var,rep(0,n))
  for(i in 1:(n-1)){
    j <- R$merge[i,1]; if(j==0) next 
    a <- ifelse(j<0,value[-j,1],value[j,2])
    j <- R$merge[i,2]; if(j==0) next 
    b <- ifelse(j<0,value[-j,1],value[j,2])
    value[i,2] <- a+b
  }
  value[n,2] <- 0
  for(i in 1:(n-1)){
    if(value[i,2]<=vmax) next
    l <- R$merge[i,1]; r <- R$merge[i,2] 
    if(l==0) a <- 0 else a <- ifelse(l<0,value[-l,1],value[l,2])
    if(r==0) b <- 0 else b <- ifelse(r<0,value[-r,1],value[r,2])
    if(min(a,b)>vmax) next
    if(a<=vmax) if(a>=vmin) {nclust <- nclust+1; mark(l,nclust)} else mark(l,0)
    if(b<=vmax) if(b>=vmin) {nclust <- nclust+1; mark(r,nclust)} else mark(r,0)    
  }  
  return(list(part=part,value=value))
}

derivedTree <- function(R,type='rank'){
  if (type == 'total') { c <- 0; ex <- expression(a+b+R$height[i]) }
  else if (type == 'count') { c <- 1; ex <- expression(a+b) }
  else  { c <- 0; ex <- expression(1+max(a,b)) }
  nm <- length(R$height)
  h <- rep(c,nm)
  for (i in 1:nm){
    j <- R$merge[i,1]; a <- ifelse(j<0,c,h[j])
    j <- R$merge[i,2]; b <- ifelse(j<0,c,h[j])
    h[i] <- eval(ex)
  }
  return(h)
}

# Clustering with relational constraint based on the class dist
cluRCdist <- function(N,D,method="max",strategy="tolerant"){
  orDendro <- function(i){if(i<0) return(-i)
    return(c(orDendro(m[i,1]),orDendro(m[i,2])))}
  idx <- function(i,j) {if (i<j) return(n*(i-1) - i*(i-1)/2 + j-i) else
    return(n*(j-1) - j*(j-1)/2 + i-j)}
  if(strategy %in% c("tolerant", "leader", "strict")){
    tol <- strategy=="tolerant"; nos <- strategy!="strict"} else 
    {cat("*** Error - Unknown strategy:", strategy,"\n"); return(NULL)}
  cat("Clustering with relational constraint based on the class dist\n")
  cat("by Vladimir Batagelj, March 2018 / July 2026\n")
  cat("Method:",method,"  Strategy:",strategy,"\n")
  if(class(D)!="dist"){cat("*** Error - D should be of class dist\n"); return(NULL)}
  print(paste("Started:",Sys.time()))
# each unit is a cluster; compute inter-cluster dissimilarity matrix
  n <- numL <- attr(D,"Size"); numLm <- numL-1
  Ro <- neighborhood(N,order=1,mode="out",V(N))
  Ri <- neighborhood(N,order=1,mode="in",V(N))
  for(i in 1:length(Ro)) { li <- length(Ro[[i]]) 
    if(li<=1) Ro[[i]] <- 0 else Ro[[i]] <- as.vector(Ro[[i]][2:li]) }
  for(i in 1:length(Ri)) { li <- length(Ri[[i]]) 
    if(li<=1) Ri[[i]] <- 0 else Ri[[i]] <- as.vector(Ri[[i]][2:li]) }
  active <- 1:numL; m <- matrix(nrow=numLm,ncol=2)
  node <- rep(0,numL); h <- numeric(numLm); w <- rep(1,numL)
  for(k in 1:numLm){
  # determine the closest pair of clusters (p,q)
    nn <- length(active); ind <- rep(Inf,nn); dd <- rep(Inf,nn)
    for(a in seq_along(active)) {i <- active[a]
    #  if((length(Ro[[i]])==1)&&(Ro[[i]][1]==0)) dd[a] <- Inf else      
      for(j in Ro[[i]]) if((j>0)&&(i!=j)) if(D[idx(i,j)] < dd[[a]]) {
        dd[[a]] <- D[idx(i,j)]; ind[[a]] <- j} }
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
    # cat("active =",active,"\n")
    Rop <- setdiff(Ro[[p]],q); Rip <- setdiff(Ri[[p]],q) 
    Roq <- setdiff(Ro[[q]],p); Riq <- setdiff(Ri[[q]],p)
    # cat("Rop:",Rop,"  Rip:",Rip,"\n")
    for(s in Riq) if(s>0) Ro[[s]] <- setdiff(Ro[[s]],q)
    r <- p; Ror <- 0; Rir <- Rip
    if(tol) Ror <- Rop else
      for(s in Rop) if(s>0) Ri[[s]] <- setdiff(Ri[[s]],p) 
    Ror <- union(Ror,Roq)
    # cat("Ror:",Ror,"\n")
    for(s in Roq) if(s>0) Ri[[s]] <- setdiff(union(Ri[[s]],r),q)
    if(nos){
      Rir <- union(Rir,Riq)
      for(s in Riq) if(s>0) Ro[[s]] <- union(Ro[[s]],r)
    }
    Ro[[r]] <- Ror; Ri[[r]] <- Rir; Ro[[q]] <- 0; Ri[[q]] <- 0
  # determine dissimilarities to the new cluster
    for(s in setdiff(active,p)){
      if(method=="max") D[idx(p,s)] <- max(D[idx(p,s)],D[idx(q,s)]) else
      if(method=="min") D[idx(p,s)] <- min(D[idx(p,s)],D[idx(q,s)]) else
      if(method=="ward") { ww <- w[p]+w[q]+w[s]
        D[idx(p,s)] <- ((w[q]+w[s])*D[idx(q,s)] + (w[p]+w[s])*D[idx(p,s)] - w[s]*dpq)/ww
      } else {cat('unknown method','\n'); return(NULL)}
    }
    w[p] <- w[q]+w[p]; node[[p]] <- k
  }
  hc <- list(merge=m,height=h,order=orDendro(numLm),labels=attr(D,"Labels"),
    method="cluRelD",call=NULL,dist.method=method,leaders=NULL)
  class(hc) <- "hclust"
  print(paste("Finished:",Sys.time()))
  return(hc)
}

# Clustering with relational constraint based on a network
cluRCnet <- function(N,method="max",strategy="tolerant"){
  orDendro <- function(i){if(i<0) return(-i)
    return(c(orDendro(m[i,1]),orDendro(m[i,2])))}
  key <- function(i,j)
    if(i<j) return(as.character(i*np+j)) else return(as.character(j*np+i))
  if(strategy %in% c("tolerant", "leader", "strict")){
    tol <- strategy=="tolerant"; nos <- strategy!="strict"} else
    {cat("*** Error - Unknown strategy:", strategy,"\n"); return(NULL)}
  cat("Clustering with relational constraint based on a dictionary\n")
  cat("by Vladimir Batagelj, March 2018 / July 2026\n")
  cat("Method:",method,"  Strategy:",strategy,"\n")
  print(paste("Started:",Sys.time()))
# network -> dict
  hD <- new.env(); n <- vcount(N); np <- n+1
  for(i in 1:ecount(N)) {e <- ends(N,i,names=FALSE)[1,]
    assign(key(e[1],e[2]),E(N)$weight[i],envir=hD)}
  attr(hD,"Size") <- n; attr(hD,"Labels") <- V(N)$name
  Ro <- neighborhood(N,order=1,mode="out",V(N))
  Ri <- neighborhood(N,order=1,mode="in",V(N))
  for(i in 1:length(Ro)) { li <- length(Ro[[i]])
    if(li<=1) Ro[[i]] <- 0 else Ro[[i]] <- as.vector(Ro[[i]][2:li]) }
  for(i in 1:length(Ri)) { li <- length(Ri[[i]])
    if(li<=1) Ri[[i]] <- 0 else Ri[[i]] <- as.vector(Ri[[i]][2:li]) }
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
    for(s in Roq) if(s>0) {Ris <- union(Ri[[s]],r); Ri[[s]] <- setdiff(Ris,q)}
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
      ### razdelaj
      #  D[idx(p,s)] <- ((w[q]+w[s])*D[idx(q,s)] + (w[p]+w[s])*D[idx(p,s)] - w[s]*dpq)/ww
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

