# cluRCdis

http://vladowiki.fmf.uni-lj.si/doku.php?id=notes:clu:cluster

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
<img width="778" alt="ITdistC" src="https://github.com/user-attachments/assets/3e140f6d-dc88-422d-b77a-3ee11d2d9306" />
 
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
 
```
 
```
 
```
 
```
 


<hr />

[igraph](./README.md), [cluRC](../README.md)
