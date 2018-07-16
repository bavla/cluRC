read_Pajek_clu <- function(f,skip=1){
# reads a partition from Pajek's clu file; skip *vertices and comments lines
  read.table(f,skip=skip,colClasses=c("integer"),header=FALSE)$V1
}

