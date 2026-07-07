# Data


```
> library(igraph)
> N <- readRDS(file=url("https://github.com/bavla/cluRC/raw/refs/heads/master/data/IT/ItalyBESsel22.rds"))
> N
IGRAPH b356491 UNW- 107 249 -- Italian provinces 2017-2025
+ attr: name (g/c), by (g/c), cdate (g/c), name (v/c), ind (v/n), lab (v/c), short (v/c),
| x (v/n), y (v/n), capital (v/c), region (v/c), macroReg (v/c), area (v/n), pop21 (v/n),
| popDensity (v/n), communes (v/n), estab (v/n), avIncome22 (v/n), hospBeds22 (v/n),
| wasteGen22 (v/n), npOrgs22 (v/n), tertEdu22 (v/n), wCouncil22 (v/n), cancerMT22 (v/n),
| weight (e/n), color (e/c)
+ edges from b356491 (vertex names):
 [1] TO--VC TO--CN TO--AT TO--AL TO--AO TO--BI VC--NO VC--AL VC--AO VC--PV VC--BI VC--VB NO--VA NO--MI
[15] NO--PV NO--VB CN--AT CN--IM CN--SV AT--AL AT--SV AL--SV AL--GE AL--PV AL--PC AO--BI IM--SV SV--GE
[29] GE--SP GE--PC GE--PR SP--PR SP--MS VA--CO VA--MI VA--VB VA--MB CO--SO CO--LC CO--MB SO--BG SO--BS
[43] SO--BZ SO--TN SO--LC MI--BG MI--PV MI--CR MI--LO MI--MB BG--BS BG--CR BG--LC BG--MB BS--CR BS--MN
+ ... omitted several edges
> sel <- c("popDensity","avIncome22","hospBeds22","wasteGen22","npOrgs22",
+    "tertEdu22","wCouncil22","cancerMT22")
> vars <- as.matrix(as_data_frame(N,what="vertices")[,sel])
> head(vars)
   popDensity avIncome22 hospBeds22 wasteGen22 npOrgs22 tertEdu22 wCouncil22 cancerMT22
TO        325    24380.8       37.8        480     63.4      32.1       34.3        8.2
VC         80    21272.3       32.1        538     85.6      25.1       33.4        7.8
NO        271    20729.6       38.4        524     66.2      22.2       34.0        7.2
CN         84    23946.6       32.6        521     79.5      22.4       29.7        7.3
AT        139    19912.3       24.5        450     77.0      26.5       29.8        7.3
AL        115    20768.5       37.2        489     69.8      24.7       29.9        8.7
> U <- scale(vars)
> head(U)
   popDensity  avIncome22  hospBeds22  wasteGen22    npOrgs22   tertEdu22  wCouncil22 cancerMT22
TO  0.1577776  1.10576909  0.72536189 -0.15409553 -0.19547329  0.81883301  0.17980312  0.1595494
VC -0.4926612  0.29583675 -0.03100367  0.53989347  1.25827515 -0.32156477 -0.04246954 -0.1768405
NO  0.0144156  0.15443405  0.80497932  0.37237888 -0.01211763 -0.79401528  0.10571223 -0.6814254
CN -0.4820418  0.99263651  0.03534419  0.33648290  0.85882175 -0.76143248 -0.95625718 -0.5973279
AT -0.3360249 -0.05851679 -1.03949109 -0.51305535  0.69511134 -0.09348521 -0.93156022 -0.5973279
AL -0.3997414  0.16456961  0.64574446 -0.04640758  0.22362536 -0.38673036 -0.90686325  0.5800368
> saveRDS(U,file="BES22selStd.rds")
> r <- hclust(d<-dist(U),method="ward.D")
> plot(r,hang=-1,cex=0.6,main="BES 2022 Ward / Free",lwd=1.2)

 
```
<img width="1000" alt="dendWardFree" src="https://github.com/user-attachments/assets/76d78199-9a9a-40a2-8157-8d4d6e0650f3" />

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
