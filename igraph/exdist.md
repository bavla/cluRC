# cluRCdis

## Clustering with RC and complete dissimilarity matrix

For illustration, we will use the data about [Italian provinces](/data/IT/README.md).
The relational constraint is given as an igraph network N. Standardized numerical node properties are given in matrix U (same order of nodes in N and U).
```
> setwd("C:/data/italy")  # use your working directory
> source("https://raw.githubusercontent.com/bavla/cluRC/refs/heads/master/igraph/cluRC.R")
> library(igraph); library(sf); library(tmap); library(spdep); library(pals); library(paletteer)
> source("./cluRCdist.R")
> N <- readRDS(file=url("https://github.com/bavla/cluRC/raw/refs/heads/master/data/IT/ItalyBESsel22.rds"))
> U <- readRDS(file=url("https://github.com/bavla/cluRC/raw/refs/heads/master/data/IT/BES22selStd.rds"))
> r <- cluRCdist(N,dist(U))
Clustering with relational constraint based on the class dist
by Vladimir Batagelj, March 2018 / July 2026
Method: max   Strategy: tolerant 
[1] "Started: 2026-07-06 03:17:30.01135"
[1] "Finished: 2026-07-06 03:17:30.174706"
> plot(r,hang=-1,cex=0.6,main="BES 2022 / RC tolerant max",lwd=1.2) 
```

<img width="1000" alt="dendroRCdTolMax" src="https://github.com/user-attachments/assets/d761c1bc-33d4-44ae-8cf8-91f7d92093e6" />

## Clusters


```
> rc <- cutree(r,k=16)
> table(rc)
rc
 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 
10 29  1  2  3  1 10  9 17  1 13  1  7  1  1  1 
> t <- table(rc)
> rc[rc %in% which(t==1)] <- 999
> p <- as.integer(factor(rc))
> names(p) <- names(rc)
> P <- st_read("./shape/prov24/ProvCM01012024_WGS84.shp")
> P$Clustering <- p
> tm_shape(P) +
+   tm_polygons("Clustering",tm_scale_categorical(values=as.vector(cols25(20))),
+   fill.legend = tm_legend(position = tm_pos_in("right", "top")) ) +
+   tm_title_out("BES 2022 / RC dis: tolerant max",
+     position = tm_pos_out("center", "top"))  
```
<img width="839" alt="ITdistCut16" src="https://github.com/user-attachments/assets/42866ce7-d9b6-4036-8653-37ed5d8ff350" />

```
> t <- transTree(r,type="rank")
> plot(t,hang=-1,cex=0.6,main="BES 2022 / RC tolerant max / rank",lwd=1.2)
> tc <- cutree(t,k=16)
> (h <- table(tc))
tc
 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 
10  6  1 18  5  2  3  1 10  9 24 14  1  1  1  1 
```
<img width="1000" alt="dendroRCdTolMaxRank" src="https://github.com/user-attachments/assets/6ca4f49f-cad7-45f0-a98d-04d1a3aba430" />

```
> tc[tc %in% which(h==1)] <- 999
> p <- as.integer(factor(tc))
> names(p) <- names(tc)
> P$Clustering <- p
> cols <- c(as.vector(paletteer_d("RColorBrewer::Set2")),"yellow","lightcyan","darkred")
> tm_shape(P) +
+   tm_polygons("Clustering",tm_scale_categorical(values=cols),
+   fill.legend = tm_legend(position = tm_pos_in("right", "top")) ) +
+   tm_title_out("BES 2022 / RC dis: tolerant max / rank",
+     position = tm_pos_out("center", "top")) 
```
<img width="824" alt="ITdistCut16rank" src="https://github.com/user-attachments/assets/dbf644f3-66ed-4b9c-a8a9-03c95cb4320c" />


## varCutree

### count
```
> r$height <- R
> rC <- varCutree(r,rep(1,107),2,30)
> table(rC$part)
 0  1  2  3  4  5  6  7 
 5 10 29 19 25 14  3  2 
> P$Clustering <- rC$part
> # cols <- c("darkred",as.vector(cols25(25))) 
> cols <- c("darkred",as.vector(paletteer_d("RColorBrewer::Set2")))
> tm_shape(P) +
+   tm_polygons("Clustering",tm_scale_categorical(values=cols),
+   fill.legend = tm_legend(position = tm_pos_in("right", "top")) ) +
+   tm_title_out("BES 2022 / RC dis: tolerant max / var count 2-30",
+     position = tm_pos_out("center", "top")) 
```
<img width="796" alt="ITdistVC30" src="https://github.com/user-attachments/assets/4ec62ee8-a521-4da0-8c6a-b8d7431f1f4d" />

### population

```
> rC <- varCutree(r,V(N)$pop21,1500000,15000000)
> table(rC$part)
 0  1  2  3  4  5  6  7 
 7 10 29 19 25 14  2  1 
> P$Clustering <- rC$part
> cols <- c("darkred",as.vector(paletteer_d("RColorBrewer::Set2")))
> tm_shape(P) +
+   tm_polygons("Clustering",tm_scale_categorical(values=cols),
+   fill.legend = tm_legend(position = tm_pos_in("right", "top")) ) +
+   tm_title_out("BES 2022 / RC dis: tolerant max / var population 1.5M-15M",
+     position = tm_pos_out("center", "top"))  
```
<img width="789" alt="ITdistVP15M" src="https://github.com/user-attachments/assets/f6408d71-f0db-4dde-9268-c02149cd3295" />

### area

```
> rC <- varCutree(r,V(N)$area,3000,60000)
> table(rC$part) 
 0  1  2  3  4  5  6  7  8  9 10 
 7  6 18  5 10 17  7 19  1 14  3 
> P$Clustering <- rC$part
> cols <- c("darkred",as.vector(paletteer_d("RColorBrewer::Set2")),"yellow","lightcyan")
> tm_shape(P) +
+   tm_polygons("Clustering",tm_scale_categorical(values=cols),
+   fill.legend = tm_legend(position = tm_pos_in("right", "top")) ) +
+   tm_title_out("BES 2022 / RC dis: tolerant max / var area 3K-60K",
+     position = tm_pos_out("center", "top")) 
>  
```
<img width="821" alt="ITdistVA60K" src="https://github.com/user-attachments/assets/474fa501-c2f2-4294-9a43-0e6521eac38a" />
 


<hr />

[igraph](./README.md), [cluRC](../README.md)
