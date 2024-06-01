################################################################
## Compute in Embedding Result for Each Cell Closest Neighbor ##
################################################################

## Libraries that need to be installed
library(spatstat.geom); library(matrixStats)
## Toy sample for cell embedding in 2D space 
cell_coor <- matrix(rnorm(20), 10, 2, dimnames=list(paste("cell", 1:10, sep=""), c("X", "Y")))
cell_coor # also see nearest neighbor function in spatstat.geom

## Compute distance matrix for all cells (items) in 2D space
dma <- spatstat.geom::pairdist(cell_coor)
rownames(dma) <- rownames(cell_coor)
colnames(dma) <- rownames(cell_coor)
dma

## Compute closest/nearest neighbor pairs (min distance)
dma[dma==0] <- NA
dma_log <- dma == matrixStats::rowMins(dma, na.rm=TRUE)
dma_log[is.na(dma_log)] <- FALSE
dma_log

## Return nearest neighbor pairs in a list
closest_pairs <- lapply(seq_along(rownames(dma_log)), function(x) unlist(dimnames(dma_log[x, dma_log[x,], drop=FALSE])))
closest_pairs

## Convert list to nearest neighbor pair character vector
closest_pairs <- sapply(closest_pairs, paste, collapse="-")
closest_pairs

