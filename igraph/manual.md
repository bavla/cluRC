# cluRC manual

<img src="https://github.com/bavla/Nets/blob/master/netsWeight/workIP.jpg" width=100 />








## cluRCdist(N,D,method="max",strategy="tolerant")
Hierarchical clustering with a relational constraint network N (igraph) on a smaller set (up to some hundreds) of units, using a complete dissimilarity matrix D (dist).
method $\in \{$ "max", "min", "average'', "ward" $\}$;
strategy $\in \{$ "tolerant", "leader", "strict" $\}$.

## cluRCnet(N,method="max",strategy="tolerant")
Hierarchical clustering with a relational constraint network N (igraph) on a larger set of units; the dissimilarity between units is provided as the weight of links.
method $\in \{$ "max", "min", "average'', "ward" $\}$;
strategy $\in \{$ "tolerant", "leader", "strict" $\}$.

## varCutree(R,var,vmin,vmax) 
An improved version of the cutree function, which, given a clustering R (hclust) and a selected variable var, determines the clusters for which the sum of its values ​​per cluster lies in the interval (vmin, vmax). The remaining units are placed in cluster 0.
[Example](./exdist.md#varcutree)

## transTree(t,type="rank")
Converts the clustering t according to the selected type $\in \{$"count", "rank", "total" $\}$ by appropriately changing the heights (using the derivedTree helper function) and the order of the units, but preserving the tree structure.

## listClu(t,tc)

## printClu(t,tc)

## push(S,input), pop(S)
Stack operation on a list S.

## orDendro(m)

## showClu(cl)

```
 
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
