# cluRC manual



## cluRCdist(N,D,method="max",strategy="tolerant")
Hierarchical clustering with a relational constraint network N (igraph) on a smaller set (up to some hundreds) of units, using a complete dissimilarity matrix D (dist).
method ∈ { "max", "min", "average", "ward" };
strategy ∈ { "tolerant", "leader", "strict" }.

## cluRCnet(N,method="max",strategy="tolerant")
Hierarchical clustering with a relational constraint network N (igraph) on a larger set of units; the dissimilarity between units is provided as the weight of links.
method ∈ { "max", "min", "average", "ward" };
strategy ∈ { "tolerant", "leader", "strict" }.

For larger networks (tens of thousands of nodes) the computation can take several hours. To track the execution, include the step parameter (e.g. step=100) - every k*step iteration prints a message.

## varCutree(R,var,vmin,vmax) 
An improved version of the cutree function, which, given a clustering R (hclust) and a selected variable var, determines the clusters for which the sum of its values ​​per cluster lies in the interval (vmin, vmax). The remaining units are placed in cluster 0.
[Example](./exdist.md#varcutree)

## transTree(t,type="rank")
Converts the clustering t according to the selected type ∈ { "count", "rank", "total" } by appropriately changing the heights (using the derivedTree helper function) and the order of the units, but preserving the tree structure.

## listClu(t,tc)
For a given clustering t and partition tc, it creates a list of vectors with the names of the units belonging to each cluster.

## printClu(t,tc)
Prints a list of clusters created with the listClu function.

## push(S,input), pop(S)
Stack operation on a list S.

## orDendro(m)
Determines the ordering of units that is compatible with the given clustering tree m.

## showClu(cl)
Displays an image of a subgraph of graph N induced by the cluster cl of the partition p.

 

<hr />

[igraph](./README.md), [cluRC](../README.md)
