# Europe


<img src="https://github.com/bavla/Nets/blob/master/netsWeight/workIP.jpg" width=100 />

http://vladowiki.fmf.uni-lj.si/doku.php?id=pro:relc:eu

Add ISO alpha-2 country codes to the The Hammond Almanac 1980 European countries data.
```
> setwd("C:/Users/vlado/work/Delphi/Cluse/Cluse/data/Europe")
> Evar <- read.csv2("Europe.csv",skip=2)
> Cc <- read.table("Ccodes.csv",header=TRUE,sep="|")
> Cc$ISOa2 <- trimws(Cc$ISOa2)
> E <- cbind(Evar$country,Cc$ISOa2,Evar[,1:24])
> colnames(E)[1:2] <- c("country","ISOa2")
> for(i in 3:26) E[,i] <- as.numeric(E[,i])
> head(E)
          country ISOa2  area population indent percUrban density popFirstCity incomPC percIndust birthRate deathRate
1         Albania    AL 11100       2650      1      34.0   239.0        6.415     647       61.5      33.3       8.1
2         Austria    AT 32375       7620      2      54.0   235.4       22.310    5491       61.0      11.0      13.5
3         Belgium    BE 11781       9890      2      87.1   839.0       10.829    7085       30.0      12.3      11.4
4        Bulgaria    BG 42823       8900      2      58.7   208.0       10.843    2799       43.0      16.1      10.7
5  Czechoslovakia    CS 49373      15100      2      66.7   306.0        7.689    4673       45.0      18.7      11.5
6         Denmark    DK 16629       5120      3      80.0   308.0       26.602    7607       35.0      12.2       9.9
  lifeExp inPHBed inPPhysic infMort illiter higEduc    roads vehicles   cars railway  radio     tv telephon newspaper
1      69     164       159    86.8      25  10.818  378.378    0.381  0.098  16.937  6.792  0.170    0.383        46
2      73      88       479    16.9       0  12.648  628.788   30.013 23.990 126.054 28.675 23.255   29.934       320
3      71     112       530    14.0       3  15.180  650.539   30.657 27.685 218.572 40.890 26.754   29.818       239
4      71     116       465    23.7       5  14.449  439.203    0.438  0.128  88.270 30.899 17.371    9.584       232
5      71      99       418    19.6       0  10.269  925.182   12.934 11.106 166.650 26.013 25.119   18.166       300
6      73     103       624     8.9       0  21.537 2405.436   31.367 26.250  95.015 36.152 31.973   48.926       341
> write.csv2(E,file="EuropePsy.csv")
```

Combine the geographical neighbors graph and European countries' data into a igraph network.
```
> library(igraph)
> netFile <- "https://github.com/bavla/cluRC/raw/refs/heads/master/data/Eu/EuropePsy.net"
> datFile <- "https://github.com/bavla/cluRC/raw/refs/heads/master/data/Eu/EuropePsy.csv"
> N <- simplify(read_graph(netFile,format="pajek"))
> V(N)$y <- 1-V(N)$y
> E <- read.csv2(datFile,row.names=1,skip=2)
> for(a in colnames(E)) vertex_attr(N,a) <- E[,a] 
> V(N)$country <- V(N)$name
> N$name <- "cluRC example from Psychometrika / Europe 1980"
> N$ref <- paste0("A. Ferligoj, V. Batagelj: Clustering with ", 
+   "relational constraints. Psychometrika 47 (4): 413-426 1982")
> N$by <- "Vladimir Batagelj"
> N$date <- date()
> saveRDS(N,file="EuropePsy.rds") 
```
 
```
> setwd("C:/Users/vlado/work/Delphi/Cluse/Cluse/data/Europe")
> source("https://raw.githubusercontent.com/bavla/cluRC/refs/heads/master/igraph/cluRC.R")
> library(igraph); library(httr)
> library(sf); library(tmap); library(spdep); library(pals); library(paletteer)
> rdsFile <- "https://github.com/bavla/cluRC/raw/refs/heads/master/data/Eu/EuropePsy.rds"
> N <- readRDS(file=url(rdsFile))
> vars <- as.matrix(as_data_frame(N,what="vertices")[,11:31]) 
> D <- dist(scale(vars))
> h <- hclust(D)
> plot(h,hang=-1,main="European countries 1980")
```

<img width="795" alt="Eu80HC" src="https://github.com/user-attachments/assets/8df4b91a-8b23-4d80-b93e-06cd0115b22f" />

 
```
> r <- cluRCdist(N,D)
Clustering with relational constraint based on the class dist
by Vladimir Batagelj, March 2018 / July 2026
Method: max   Strategy: tolerant 
[1] "Started: 2026-07-10 00:23:32.045258"
[1] "Finished: 2026-07-10 00:23:32.107539"
> plot(r,hang=-1,main="European countries 1980 / RC tolerant max")

```
<img width="846" alt="Eu80RCtolMax" src="https://github.com/user-attachments/assets/9e68c07e-905b-4e85-ac7b-73ea2fdc75be" />
 
```
> library(paletteer)
> cols <- as.vector(paletteer_d("RColorBrewer::Set2"))
> V(N)$name <- V(N)$ISOa2
> plot(N,vertex.size=8,vertex.label.cex=2,vertex.color=cols[tc])
> plot(N,vertex.size=14,vertex.label.cex=1,vertex.color=cols[tc],edge.width=3)
> plot(t,hang=-1,main="Eu 1980 / RC tolerant max / count")
> printClu(t,tc)
1 : Albania 
2 : Austria Czechoslovakia Hungary Italy Poland USSR 
3 : Belgium France East_Germany West_Germany Luxembourg Netherlands Switzerland United_Kingdom 
4 : Bulgaria Greece Romania Yugoslavia 
5 : Denmark Finland Iceland Norway Sweden 
6 : Ireland Portugal Spain  
```

<img width="837" alt="Eu80RCtolMaxC" src="https://github.com/user-attachments/assets/a1dea464-704d-4899-8a28-b0110b7f5bba" />

<img width="800" alt="Eu80netC6" src="https://github.com/user-attachments/assets/fc7b4e01-3bcf-42d8-a908-38d74c9f96d2" />


<hr />

[Datasets](../README.md), [cluRC](../../README.md)
