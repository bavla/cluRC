# Psychometrika

Data of the example from the paper  

A. Ferligoj, V. Batagelj: Some types of clustering with relational constraints. Psychometrika 48(1983)4, 541-552.

are available on the files [SomeTyXY.net](https://github.com/bavla/cluRC/blob/master/SomeTyXY.net) and [SomeTy.dis](https://github.com/bavla/cluRC/blob/master/SomeTy.dis). To use them with cluRC we have to transform them into an igraph network and a distance lower triangular matrix.

## Converting data into cluRC format


```
> library(igraph)
> netFile <- "https://raw.githubusercontent.com/bavla/cluRC/refs/heads/master/SomeTyXY.net"
> disFile <- "https://raw.githubusercontent.com/bavla/cluRC/refs/heads/master/SomeTy.dis"
> N <- read_graph(netFile,format="pajek")
> N$name <- "cluRC example from Psychometrika"
> N$ref <- paste0("A. Ferligoj, V. Batagelj: Some types of clustering with ", 
+   "relational constraints. Psychometrika 48(1983)4, 541-552.")
> N$by <- "Vladimir Batagelj"
> N$date <- date()
> N
IGRAPH 55b1b08 DNW- 10 17 -- cluRC example from Psychometrika
+ attr: name (g/c), ref (g/c), by (g/c), date (g/c), id (v/c), name (v/c), x (v/n), y (v/n),
| z (v/n), weight (e/n)
+ edges from 55b1b08 (vertex names):
 [1] b->a b->d c->a c->b c->e c->f d->h e->g f->d f->e f->g g->j h->f h->i h->j i->d j->g
> plot(N)
> L <- readLines(url(disFile))
> S <- lapply(strsplit(trimws(L),"\\s+"),as.numeric)
> n <- length(S); T <- matrix(0,nrow=n,ncol=n)
> rownames(T) <- colnames(T) <- V(N)$name
> for (i in 1:n) T[i,i:n] <- S[[i]]
> D <- as.dist(t(T))
> D
  a b c d e f g h i
b 5                
c 7 8              
d 4 1 4            
e 6 3 5 3          
f 6 4 7 2 6        
g 2 4 9 4 4 6      
h 4 5 3 8 6 8 4    
i 2 3 2 6 5 5 8 3  
j 3 4 5 3 7 7 2 4 5
> saveRDS(N,file="SomeTyNet.rds"); saveRDS(D,file="SomeTyDis.rds") 
```

<img width="476"  alt="SomeTyNet" src="https://github.com/user-attachments/assets/17ed97cf-c418-4247-a06a-3266f27aeb53" />

## Apply the strategies

Strategies: tolerant, leader, strict

Methods: max, min, ward

```
> source("https://raw.githubusercontent.com/bavla/cluRC/refs/heads/master/igraph/cluRC.R")
> library(igraph)

> N <- readRDS(file="SomeTyNet.rds"); D <- readRDS(file="SomeTyDis.rds")
> # tolerant
> r <- cluRCdist(N,D)
> plot(r,hang=-1,cex=1.5,main="Some types: tolerant max")
> r <- cluRCdist(N,D,method="ward",strategy="tolerant")
> plot(r,hang=-1,cex=1.5,main="Some types: tolerant ward")
```

<img width="500" alt="SomeTyTolMax" src="https://github.com/user-attachments/assets/172a0cb7-d1f4-4d29-9e1c-c86a45ffc5d5" />
<img width="500" alt="SomeTyTolWard" src="https://github.com/user-attachments/assets/b5328907-4740-4475-aff2-995c127da5cc" />

 
```
> # leader
> r <- cluRCdist(N,D,strategy="leader")
> plot(r,hang=-1,cex=1.5,main="Some types: leader max")
> # strict
> r <- cluRCdist(N,D,method="max",strategy="strict")
> plot(r,hang=-1,cex=1.5,main="Some types: strict max") 
```
<img width="500" alt="SomeTyLdrMax" src="https://github.com/user-attachments/assets/8882e2ac-2936-4305-ba64-69ba5221759d" />
<img width="500" alt="SomeTyStrMax" src="https://github.com/user-attachments/assets/1076e64f-76f2-4cda-af31-8e01bd679e54" />
 
```
 
```
 
```
 
```
 

<hr />

[Datasets](../README.md), [cluRC](../../README.md)
