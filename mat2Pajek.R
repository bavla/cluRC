mat2Pajek <- function(A,file,twomode=FALSE){
  library(Matrix)
  tit <- attr(A,"header")$description; net <- file(file,"w")
  cat("MATLAB to Pajek\n",tit,"\nTransformed: ",date(),"\n",sep="")
  cat("% MATLAB to Pajek\n% ",tit,"\n% Transformed: ",date(),"\n",sep="",file=net)
  M <- A$A; T <- as(M,'dgTMatrix'); D <- A$local.info
  nr <- T@Dim[1]; nc <- T@Dim[2]
  cat("*network fromMATLAB\n",file=net)
  if(twomode){
    cat("*vertices",nr+nc,nr,"\n",file=net)
    if(!is.null(T@Dimnames[[1]])){N <- T@Dimnames[[1]]
      for(i in seq_along(N)) cat(i,' "',N[i],'"\n',sep="",file=net)}
    if(!is.null(T@Dimnames[[2]])){N <- T@Dimnames[[2]]
      for(i in seq_along(N)) cat(i+nr,' "',N[i],'"\n',sep="",file=net)}
  } else {
    cat("*vertices",nr,"\n",file=net)
    if(!is.null(T@Dimnames[[1]])){N <- T@Dimnames[[1]]
      for(i in seq_along(N)) cat(i,' "',N[i],'"\n',sep="",file=net)}
  }
  cat("*arcs\n ",file=net)
  if(twomode) cat(paste(T@i+1,nr+T@j+1,T@x,"\n"),file=net) else
    cat(paste(T@i+1,T@j+1,T@x,"\n"),file=net)
  if(!is.null(D)){
    nr <- nrow(D); nc <- ncol(D)
    for(j in 1:nc){
      cat("\n*vector V",j,"\n*vertices ",nr,"\n ",sep="",file=net)
      cat(paste(D[,j],"\n"),file=net)
    }
  }
  close(net)
}

