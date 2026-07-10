# Authors in the field of clustering and BM

<img src="https://github.com/bavla/Nets/blob/master/netsWeight/workIP.jpg" width=100 />

## Converting ACiAn network from Pajek to igraph format
```
> setwd("C:/Users/vlado/DL/data/italy")
> source("https://raw.githubusercontent.com/bavla/cluRC/refs/heads/master/igraph/cluRC.R")
> library(igraph)

> netFile <- "C:/Users/vlado/docs/papers/2026/ifcs/slides/BM/ACiAn.net"
> N <- read_graph(netFile,format="pajek")
> N$name <- "cluRC example from Clustering and blockmodeling"
> N$ref <- paste0("P. Doreian, V. Batagelj, A. Ferligoj: Advances in ", 
+   "Network Clustering and Blockmodeling. Wiley 2020.")
> N$by <- "Vladimir Batagelj"
> N$date <- date()
> saveRDS(N,file="ACiAnNet.rds")
```
## Removing loops and converting weight to dissimilarity
```
> N <- simplify(N) # remove loops
> E(N)$cite <- E(N)$weight
> summary(E(N)$weight)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000116 0.0018215 0.0045662 0.0083398 0.0095238 2.0000000 
> wmax <- max(E(N)$weight)
> E(N)$weight <- 1 - E(N)$weight/wmax
> N
IGRAPH a2e04d6 DNW- 62143 642419 -- cluRC example from Clustering and blockmodeling
+ attr: name (g/c), ref (g/c), by (g/c), date (g/c), id (v/c), name (v/c), x
| (v/n), y (v/n), z (v/n), weight (e/n), cite (e/n)
+ edges from a2e04d6 (vertex names):
 [1] [ANONYMO_ ->SMITHSON_F HARRASSO_H->SMITHSON_F FORBES_M  ->GESELL_A  
 [4] THOMPSON_W->BOSE_R     THOMPSON_W->CLATWORT_W THOMPSON_W->CONNOR_W  
+ ... omitted several edges
```
## Analysis

```
showClu <- function(cl){
  C <- delete_vertices(N, which(V(N)$p!=cl))
  xy <- layout_with_fr(C)
  plot(C,layout=xy,vertex.size=10,vertex.label.cex=0.6,edge.width=10*E(C)$weight,
    main=paste0("Cluster ",cl))
} 
```

```
> N <- readRDS(file="ACiAnNet.rds")
> summary(E(N)$weight)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000116 0.0018382 0.0046296 0.0085531 0.0096154 2.8788444 
> N <- simplify(N) # remove loops
> N <- delete.edges(N, which(E(N)$weight<0.05))
> N <- delete_vertices(N, which(degree(N)==0))
> vcount(N)
> E(N)$cite <- E(N)$weight
> wmax <- 1.01*max(E(N)$weight)
> E(N)$weight <- 1 - E(N)$weight/wmax
> r <- cluRCnet(N,strategy="leader",step=100)
> E(N)$weight <- E(N)$cite

> n <- length(r$order)
> rC <- varCutree(r,rep(1,n),10,60)
> V(N)$p <- rC$part
> table(rC$part)
     0      1      2      3      4      6      7      8      9     10     11     12     13     14     31 999999 
  2652     60     60     60     22     15     16     10     10     10     10     10     25     34     60   2889
> # inspect selected clustering
> showClu(3)
```
 <img width="839" alt="Cluster4" src="https://github.com/user-attachments/assets/f3bcd3b9-4335-408c-bd24-a5f2728eb5cf" />

 
```
> cl <- 4
> C <- delete_vertices(N, which(V(N)$p!=cl))
> C <- delete_vertex_attr(delete_vertex_attr(C,"x"),"y")
> Pt <- tkplot(C,800,800,edge.curved=0,vertex.size=10,vertex.label.cex=0.6)
> coor <- tk_coords(Pt,norm=F) # save new coordinates
> # tkplot window is still active
> tk_close(Pt)
> V(C)$x <- coor[,1]; V(C)$y <- coor[,2]
> saveRDS(C,file="AciA-4.rds")
 
```
 
```
 
```
 


<hr />

[Datasets](../README.md), [cluRC](../../README.md)

