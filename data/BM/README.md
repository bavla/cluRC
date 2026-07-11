# Authors in the field of clustering and BM


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
> N <- readRDS(file="ACiAnNet.rds")
> N
IGRAPH a2e04d6 DNW- 62143 642419 -- cluRC example from Clustering and blockmodeling
+ attr: name (g/c), ref (g/c), by (g/c), date (g/c), id (v/c), name (v/c), x
| (v/n), y (v/n), z (v/n), weight (e/n), cite (e/n)
+ edges from a2e04d6 (vertex names):
 [1] [ANONYMO_ ->SMITHSON_F HARRASSO_H->SMITHSON_F FORBES_M  ->GESELL_A  
 [4] THOMPSON_W->BOSE_R     THOMPSON_W->CLATWORT_W THOMPSON_W->CONNOR_W  
+ ... omitted several edges
> summary(E(N)$weight)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000116 0.0018382 0.0046296 0.0085531 0.0096154 2.8788444 
> N <- simplify(N) # remove loops
> N <- delete_edges(N, which(E(N)$weight<0.05))
> N <- delete_vertices(N, which(degree(N)==0))
> vcount(N)
> E(N)$cite <- E(N)$weight
> wmax <- 1.01*max(E(N)$weight)
> E(N)$weight <- 1 - E(N)$weight/wmax
> r <- cluRCnet(N,strategy="leader",step=100)
Clustering with relational constraint based on a dictionary
by Vladimir Batagelj, March 2018 / July 2026
Method: max   Strategy: leader 
[1] "Started: 2026-07-10 20:56:16.17524"
Fri Jul 10 20:56:55 2026 dictionary
Fri Jul 10 20:56:56 2026 out neighbors
Fri Jul 10 20:56:58 2026 in neighbors
Fri Jul 10 20:57:05 2026  n = 100 
Fri Jul 10 20:57:12 2026  n = 200 
Fri Jul 10 20:57:19 2026  n = 300 
Fri Jul 10 20:57:25 2026  n = 400
...
Fri Jul 10 20:58:47 2026  n = 2600 
Fri Jul 10 20:58:47 2026  n = 2700 
several components 2732 3 Inf 
Create clustering
[1] "Finished: 2026-07-10 20:58:47.797469"
> E(N)$weight <- E(N)$cite
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
> n <- length(r$order)
> rC <- varCutree(r,rep(1,n),10,60)
> V(N)$p <- rC$part
> table(rC$part)
     0      1      2      3      4      6      7      8      9     10     11     12     13     14     31 999999 
  2652     60     60     60     22     15     16     10     10     10     10     10     25     34     60   2889
> # inspect selected clustering
> showClu(3)
```
 <img width="600" alt="Cluster4" src="https://github.com/user-attachments/assets/f3bcd3b9-4335-408c-bd24-a5f2728eb5cf" />

 
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
> C <- readRDS(file="AciA-4.rds")
> plot(C,vertex.size=10,vertex.label.cex=0.6,edge.width=10*E(C)$weight,
+     main=paste0("ACiAn cluster 4")) 
```
<img width="500" alt="ACiAnCluster4" src="https://github.com/user-attachments/assets/fc951eda-438f-4cce-accd-7f0185eea596" />
<img width="500" alt="ACiAnCluster7" src="https://github.com/user-attachments/assets/deec2f65-6c2c-4e90-9ab8-1746353c5fa1" />

<img width="500" alt="ACiAnCluster13" src="https://github.com/user-attachments/assets/2aed7f60-7e35-4482-8b9e-9ea8313a6e3d" />
<img width="500" alt="ACiAnCluster14" src="https://github.com/user-attachments/assets/1a294e3f-110e-4ea1-81fa-d50292a45dec" />


<hr />

[Datasets](../README.md), [cluRC](../../README.md)

