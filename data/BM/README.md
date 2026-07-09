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
 [7] THOMPSON_W->HUNT_G     THOMPSON_W->LEHMANN_E  THOMPSON_W->MADOW_W   
[10] THOMPSON_W->SKIBINSK_M THOMPSON_W->WALD_A     THOMPSON_W->YOUDEN_W  
[13] WALD_A    ->WHITE_H    WALD_A    ->BOORMAN_S  WALD_A    ->BREIGER_R 
[16] WALD_A    ->COLEMAN_J  WALD_A    ->WILLIAMS_O WALD_A    ->KNOKE_D   
[19] WALD_A    ->CARRINGT_P WALD_A    ->MARSH_D    WALD_A    ->WASSERMA_S
+ ... omitted several edges
```
## Analysis

```

> r <- cluRCnet(N,strategy="leader")
 
 
```
 
```
 
```
 
```
 
```
 
```
 
```
 


<hr />

[Datasets](../README.md), [cluRC](../../README.md)

