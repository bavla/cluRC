# clustAna()
#
# analysis of clusters
# Vladimir Batagelj, July 16, 2018

clustAna <- function(P,clu,all=FALSE){
  k <- max(clu); nVar <- ncol(P)
  mP <- matrix(0,nrow=k,ncol=nVar)
  for(i in 1:k) mP[i,] <- colMeans(P[clu==i,,drop=FALSE])
  rownames(mP) <- paste("C",1:k,sep=""); colnames(mP) <- colnames(P)
  if(all){mP <- rbind(mP,colMeans(P)); rownames(mP)[k+1] <- "all"}
  return(mP)
}

# clustAna(P,mtK6,all=TRUE)
# clustAna(ZP,mtK6)