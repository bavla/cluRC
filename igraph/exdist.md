# cluRCdis

http://vladowiki.fmf.uni-lj.si/doku.php?id=notes:clu:cluster

```
> setwd("C:/data/italy")  # use your working directory
> source("https://raw.githubusercontent.com/bavla/cluRC/refs/heads/master/igraph/cluRC.R")
> library(igraph); library(sf); library(tmap); library(spdep); library(pals)
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
 
```
 
```
 
```
 
```
 
```
 
```
 
```
 


<hr />

[igraph](./README.md), [cluRC](../README.md)
