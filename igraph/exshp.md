# Network from a shape file

```
===== Application =====


<code>
> setwd("C:/Users/vlado/DL/data/italy")
> library(sf)
> library(spdep)
> polygons <- st_read("./shape/prov24/ProvCM01012024_WGS84.shp")
Reading layer `ProvCM01012024_WGS84' from data source 
  `C:\Users\vlado\DL\data\italy\shape\prov24\ProvCM01012024_WGS84.shp' using driver `ESRI Shapefile'
Simple feature collection with 107 features and 12 fields
Geometry type: MULTIPOLYGON
Dimension:     XY
Bounding box:  xmin: 313279.3 ymin: 3933683 xmax: 1312016 ymax: 5220292
Projected CRS: WGS 84 / UTM zone 32N
> neighbors <- poly2nb(polygons,queen=TRUE)
Warning message:
In poly2nb(polygons, queen = TRUE) : neighbour object has 3 sub-graphs;
if this sub-graph count seems unexpected, try increasing the snap argument.
> summary(neighbors)
Neighbour list object:
Number of regions: 107 
Number of nonzero links: 476 
Percentage nonzero weights: 4.157568 
Average number of links: 4.448598 
3 disjoint connected subgraphs
Link number distribution:

 1  2  3  4  5  6  7  8  9 
 2 11 19 22 27 14  9  2  1 
2 least connected regions:
32 92 with 1 link
1 most connected region:
48 with 9 links
> rook <- poly2nb(polygons,queen=FALSE)
Warning message:
In poly2nb(polygons, queen = FALSE) : neighbour object has 3 sub-graphs;
if this sub-graph count seems unexpected, try increasing the snap argument.
> summary(rook)
Neighbour list object:
Number of regions: 107 
Number of nonzero links: 476 
Percentage nonzero weights: 4.157568 
Average number of links: 4.448598 
3 disjoint connected subgraphs
Link number distribution:

 1  2  3  4  5  6  7  8  9 
 2 11 19 22 27 14  9  2  1 
2 least connected regions:
32 92 with 1 link
1 most connected region:
48 with 9 links

> centroids <- st_centroid(st_geometry(polygons))
> coordinates <- st_coordinates(centroids)
> plot(st_geometry(polygons), border = "grey", col = "lightgray")
> plot(neighbors, coordinates, add = TRUE, col = "blue", lwd = 2, pch = 19)
</code>


{{work:data:nets:pics:italynet.png?600}}

{{work:data:nets:pics:italynet.pdf|italynet.pdf}}


<code>
> library(igraph)
> source("https://raw.githubusercontent.com/bavla/Rnet/master/R/igraph+.R")
> str(polygons)
Classes ‘sf’ and 'data.frame':  107 obs. of  13 variables:
 $ COD_RIP   : num  1 1 1 1 1 1 1 1 1 1 ...
 $ COD_REG   : num  1 1 1 1 1 1 2 7 7 7 ...
 $ COD_PROV  : num  1 2 3 4 5 6 7 8 9 10 ...
 $ COD_CM    : num  201 0 0 0 0 0 0 0 0 210 ...
 $ COD_UTS   : num  201 2 3 4 5 6 7 8 9 210 ...
 $ DEN_PROV  : chr  "-" "Vercelli" "Novara" "Cuneo" ...
 $ DEN_CM    : chr  "Torino" "-" "-" "-" ...
 $ DEN_UTS   : chr  "Torino" "Vercelli" "Novara" "Cuneo" ...
 $ SIGLA     : chr  "TO" "VC" "NO" "CN" ...
 $ TIPO_UTS  : chr  "Città metropolitana" "Provincia" "Provincia" "Provincia" ...
 $ Shape_Leng: num  593230 458711 276707 542032 356341 ...
 $ Shape_Area: num  6.83e+09 2.08e+09 1.34e+09 6.89e+09 1.51e+09 ...
 $ geometry  :sfc_MULTIPOLYGON of length 107; first list element: List of 1
  ..$ :List of 2
  .. ..$ : num [1:18330, 1:2] 411015 411070 411132 411266 411330 ...
  .. ..$ : num [1:229, 1:2] 378673 378671 378670 378675 378698 ...
  ..- attr(*, "class")= chr [1:3] "XY" "MULTIPOLYGON" "sfg"
 - attr(*, "sf_column")= chr "geometry"
 - attr(*, "agr")= Factor w/ 3 levels "constant","aggregate",..: NA NA NA NA NA NA NA NA NA NA ...
  ..- attr(*, "names")= chr [1:12] "COD_RIP" "COD_REG" "COD_PROV" "COD_CM" ...
> neighbors[[1]]
[1]  2  4  5  6  7 96
> neighbors[[2]]
[1]   1   3   6   7  18  96 103
> m <- 476; k <- 0; n <- 107
> tail <- head <- rep(0,m)
> for(i in 1:n) { L <- neighbors[[i]]; nL <- length(L)
+   for(j in 1:nL) { h <- L[j]
+     if(i <= h) { k <- k+1; tail[k] <- i; head[k] <- h } }
+ }
> k
[1] 238
> links <- data.frame(from=tail[1:k],to=head[1:k],weight=rep(1,k))
> head(links)
  from to weight
1    1  2      1
2    1  4      1
3    1  5      1
4    1  6      1
5    1  7      1
6    1 96      1

> xy <- as.matrix(coordinates)
> head(xy)
            X       Y
[1,] 377217.0 5000174
[2,] 438311.3 5041480
[3,] 465200.2 5045628
[4,] 387693.8 4925820
[5,] 435743.5 4969641
[6,] 473418.1 4963911
> nodes <- data.frame(name=1:n,ind=polygons$COD_PROV,lab=polygons$DEN_UTS,
+   short=polygons$SIGLA,x=xy[,1],y=xy[,2])
> N <- graph_from_data_frame(links,directed=FALSE,vertices=nodes)
> N
IGRAPH 8de3d89 UNW- 107 238 -- 
+ attr: name (v/c), ind (v/n), lab (v/c), short (v/c), x (v/n), y (v/n), weight (e/n)
+ edges from 8de3d89 (vertex names):
  [1] 1 --2   1 --4   1 --5   1 --6   1 --7   1 --96  2 --3   2 --6   2 --7   2 --18  2 --96  2 --103
 [13] 3 --12  3 --15  3 --18  3 --103 4 --5   4 --8   4 --9   5 --6   5 --9   6 --9   6 --10  6 --18 
 [25] 6 --33  7 --96  8 --9   9 --10  10--11  10--33  10--34  11--34  11--45  12--13  12--15  12--103
 [37] 12--104 13--14  13--97  13--104 14--16  14--17  14--21  14--22  14--97  15--16  15--18  15--19 
 [49] 15--98  15--104 16--17  16--19  16--97  16--104 17--19  17--20  17--22  17--23  18--33  18--98 
 [61] 19--20  19--33  19--34  19--98  20--23  20--29  20--34  20--35  20--36  20--38  21--22  21--25 
 [73] 22--23  22--24  22--25  23--24  23--28  23--29  24--25  24--26  24--28  25--26  25--30  25--93 
 [85] 26--27  26--28  26--93  27--28  27--29  27--30  27--93  28--29  29--38  30--31  30--93  31--32 
+ ... omitted several edges
> V(N)$name <- V(N)$short
> N$name <- "Italian provinces 2017-2025"
> N$by <- "Vladimir Batagelj"
> N$cdate <- date()
> N
IGRAPH 8de3d89 UNW- 107 238 -- Italian provinces 2017-2025
+ attr: name (g/c), by (g/c), cdate (g/c), name (v/c), ind (v/n), lab (v/c), short (v/c), x
| (v/n), y (v/n), weight (e/n)
+ edges from 8de3d89 (vertex names):
  [1] TO--VC TO--CN TO--AT TO--AL TO--AO TO--BI VC--NO VC--AL VC--AO VC--PV VC--BI VC--VB NO--VA NO--MI
 [15] NO--PV NO--VB CN--AT CN--IM CN--SV AT--AL AT--SV AL--SV AL--GE AL--PV AL--PC AO--BI IM--SV SV--GE
 [29] GE--SP GE--PC GE--PR SP--PR SP--MS VA--CO VA--MI VA--VB VA--MB CO--SO CO--LC CO--MB SO--BG SO--BS
 [43] SO--BZ SO--TN SO--LC MI--BG MI--PV MI--CR MI--LO MI--MB BG--BS BG--CR BG--LC BG--MB BS--CR BS--MN
 [57] BS--TN BS--VR PV--PC PV--LO CR--MN CR--PC CR--PR CR--LO MN--VR MN--RO MN--PR MN--RE MN--MO MN--FE
 [71] BZ--TN BZ--BL TN--VR TN--VI TN--BL VR--VI VR--PD VR--RO VI--BL VI--TV VI--PD BL--TV BL--UD BL--PN
 [85] TV--VE TV--PD TV--PN VE--PD VE--RO VE--UD VE--PN PD--RO RO--FE UD--GO UD--PN GO--TS PC--PR PC--LO
+ ... omitted several edges
> saveRDS(N,file="ItalyIgraph.rds")
> # N <- readRDS(file="ItalyIgraph.rds")
> write.graph.paj(N,file="ItalyIgraph.paj",coor=cbind(V(N)$x,V(N)$y),weight="weight")
</code>

<code>
> plot(N,vertex.size=2.5,vertex.label.cex=0.15)
</code>
{{work:data:nets:pics:italyigraph2.png}}

{{work:data:nets:pics:italyigraph2.pdf|italyigraph2.pdf}}

<code>
> N <- set_edge_attr(N,"color",value="red")
> N <- add_edges(N,c(
+  8, 90 ,   9, 90 ,  49, 90 ,  53, 90 ,  53, 91 ,  56, 91 ,
+ 81, 92 ,  81, 107,  83, 102,  80, 83 ,  80, 87 ), color="blue", weight=1)
> plot(N,vertex.size=2.5,vertex.label.cex=0.15)
> saveRDS(N,file="ItalyIgraphExt.rds")
</code>

{{work:data:nets:pics:italyigraphext.png}}
 
```
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
