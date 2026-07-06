# cluRC
Clustering with relational constraint

  * [cluRC 2018 version](./2018.md)
  * [cluRc 2026 based on igraph](./igraph/README.md)
  * [IFCS 2026](./ifcs26/README.md)
  * [cluRC datasets](./data/README.md)
  * [cluRC documents](./doc/README.md)

Let $U$ be a set of units, each measured on several properties (attributes) $A$, from which we can compute the dissimilarity $d$ between pairs of units. We also have a relation $R$ on $U$ that indicates the compatibility of units. This defines a network $N = (U, R, A, d)$, with $U$ as the set of nodes,  $R$ as the set of links, and $d$ as the weight on links.

In clustering with a relational constraint, we aim to find a clustering / partition $\mathbf{C} = \{C_1,\ldots,C_k\}$ of the set of units $U$ that minimizes the criterion function, also known as the clustering "error", $P(\mathbf{C}) = \sum_{C \in \mathbf{C}} p(C)$, where $p(C)$ is a cluster error. Different $p(C)$ s can be found in the clustering literature. Here, each cluster $C \in \mathbf{C}$ must form a connected component of the graph $(U,R)$ in a chosen way -- strict, leader, tolerant, or two-way. These modes define the set of feasible clusterings $\Phi$. This can be expressed as an optimization problem [2]: find a clustering $\mathbf{C}^* \in  \Phi$ such that 
$$P(\mathbf{C}^*) = \min_{\mathbf{C} \in \Phi} P(\mathbf{C}) .$$ 
A more complex criterion function $P(\mathbf{C}) = \sum_{C_1, C_2 \in \mathbf{C}} q(C_1, C_2)$ is used in blockmodeling [6,5]. We can also cluster the links of a network [3].

\centerline{\includegraphics[width=120mm,viewport=80 70 1050 680,clip]{cut50_20.pdf}}

A typical example of tolerant clustering with a relational constraint is the regionalization problem, in which we combine basic administrative geographical units into a smaller number of clusters or regions composed of similar and geographically connected units. 

An example of leader clustering is identifying thematic clusters of authors in a weighted citation network between authors. The figure shows a cluster of authors from the field of traditional clustering, with Ward as the leader, based on the bibliography of clustering and block modeling. This cluster includes "classics" such as Rand, Johnson, Rohlf, Hartigan, Gordon, Hubert, and DeSarbo, among others.

Clustering with a relational constraint was introduced by Ferligoj and Batagelj [7,8]. They also proposed solution procedures based on a dissimilarity matrix. Both agglomerative and local optimization approaches can be adapted. However, these procedures are not suitable for networks with a very large number of units (over 10,000). Large networks are typically sparse, and by limiting ourselves to dissimilarities only between linked units and to criterion functions such as minimum, maximum, or average, we can develop much more efficient procedures [5].

Procedures for clustering with relational constraints are partially supported in Pajek, a program for analysis and visualization of large networks [4]. The Pajek or netsJSON format [1] is used to describe network data. cluRC is an R package that offers both types of procedures to users. In this talk, we will present these procedures and demonstrate their application to real-world data. The cluRC package is available on GitHub [9].

## References
[1] Batagelj, V., Pisanski, T., Savnik, I., Slavec, A., Bašić, N.: Towards a Format for Describing Networks. In
Proceedings of 28th International Multiconference IS 2025, 6.-10. October 2025, Ljubljana, Slovenia; SiKDD, 102--105. IJS, Ljubljana (2025).

[2]
Batagelj, V., Ferligoj, A.: Constrained clustering problems. In Advances in Data Science and Classification: Proceedings of the 6th Conference of the International Federation of Classification Societies (IFCS-98) Università “La Sapienza”, Rome, 21--24 July, 1998, pp. 137--144. Berlin, Heidelberg: Springer Berlin Heidelberg (1998).

[3]
Bodlaj, J., Batagelj, V.: Hierarchical link clustering algorithm in networks. Physical Review E 91(6): 062814  (2015).

[4]
De Nooy, W., Mrvar, A., Batagelj, V.: Exploratory social network analysis with Pajek: Revised and expanded edition for updated software. Vol. 46. Cambridge University Press (2018).

[5]
Doreian, P., Batagelj, V., Ferligoj, A., eds.: Advances in network clustering and blockmodeling. John Wiley \& Sons (2020).

[6]
Doreian, P., Batagelj, V., Ferligoj, A.: Generalized blockmodeling. No. 25. Cambridge University Press (2005).

[7]
Ferligoj, A., Batagelj, V.: Clustering with relational constraint. Psychometrika 47(4), 413--426 (1982).

[8]
Ferligoj, A., Batagelj, V.: Some types of clustering with relational constraints. Psychometrika 48(4),  541--552 (1983).

[9]
cluRC  repository, https://github.com/bavla/cluRC , last accessed 2026/1/18.

