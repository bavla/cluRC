# colDiss()
# Color plots of a dissimilarity matrix, without and with ordering
# http://ichthyology.usm.edu/courses/multivariate/coldiss.R
#
# License: GPL-2 
# Author: Francois Gillet, August 2009
# addapted by: Vladimir Batagelj, July 16, 2018
#

"colDiss" <- function(D, nc = 4, byrank = TRUE, diag = FALSE, order = NULL, ...)
{
	require(gclus)
	D <- D/max(D)
      spe.color <- dmat.color(1-D, byrank=byrank, cm.colors(nc))
	if(is.null(order)) spe.o <- order.single(1-D) else spe.o <- order
	speo.color <- spe.color[spe.o,spe.o]	
	if (diag) {plotcolors(speo.color, rlabels=attributes(D)$Labels[spe.o],  
	    dlabels=attributes(D)$Labels[spe.o], ...)
	} else plotcolors(speo.color, rlabels=attributes(D)$Labels[spe.o], ...)
}

# colDiss(D,label.cex=0.5,order=1:48,main="Unordered Dissimilarity Matrix")
# colDiss(D,label.cex=0.5,order=ward$order,main="Ordered Dissimilarity Matrix/Ward")


