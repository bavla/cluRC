# Italian provinces map



```
> library(igraph); library(sf); library(tmap); library(pals)
> N <- readRDS(file=url("https://github.com/bavla/cluRC/raw/refs/heads/master/data/IT/ItalyBESsel22.rds"))
> P <- st_read("./shape/prov24/ProvCM01012024_WGS84.shp")
> tm_shape(P) +
+   tm_polygons("COD_REG",tm_scale_categorical(values=as.vector(cols25(25)))) +
+   tm_layout(legend.outside=TRUE)
> P$region <- V(N)$region
> tm_shape(P) +
+   tm_polygons("region",tm_scale_categorical(values=as.vector(cols25(25)))) +
+   tm_layout(legend.outside=TRUE) 
```

<img width="825" height="794" alt="ITregions" src="https://github.com/user-attachments/assets/8bec775f-6b94-4a4d-ad49-c99788dd4c5c" />

```
 
```
 
```
 
```
 
```
 
```
 
```
 
```
 


<hr />

[Italy](./README.md), [Datasets](../README.md), [cluRC](../../README.md)
