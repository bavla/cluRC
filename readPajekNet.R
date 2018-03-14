read_Pajek_net <- function(f,skip=0){
# reads a network from Pajek's net file; skip initial comments lines
   L <- readLines(f)
   st <- grep("\\*",L)
   S <- unlist(strsplit(L[st[1]]," "))
   n <- as.integer(S[2]); n1 <- st[1]+1; n2 <- st[2]-1
   m1 <- st[2]+1; m2 <- length(L); m <- m2-m1+1
   Names <- unlist(strsplit(L[n1:n2],'"'))[3*(1:n)-1]
   R <- vector(mode="list",length=n)
   S <- unlist(strsplit(L[m1:m2],'[[:space:]]+'))
   b <- as.integer(S[4*(1:m)-2]); e <- as.integer(S[4*(1:m)-1])
   for(k in 1:m) R[[b[[k]]]] <- union(R[[b[[k]]]],e[[k]])
   names(R) <- Names
   return(R)
}
