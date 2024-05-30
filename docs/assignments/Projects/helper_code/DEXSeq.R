#####################################
## Run DESseq on Howard et al Data ##
#####################################
## Helper code (May 25, 2024)

## Create flattened exon ranges 
library(GenomicFeatures)  # for the exonicParts() function
library(rtracklayer)        # for the makeTxDbFromGFF() function
txdb <- loadDb("./data/tair10.sqlite")
flattenedAnnotation2 <- exonicParts(txdb, linked.to.single.gene.only=TRUE)
names(flattenedAnnotation2) <- sprintf("%s:E%0.3d", flattenedAnnotation2$gene_id, flattenedAnnotation2$exonic_part)

## Exon-level read counting in parallel mode 
library(Rsamtools)
library(GenomicAlignments)
outpaths <- list.files('results/hisat2_mapping', pattern='sorted.bam$', full.names=TRUE)
file.exists(outpaths) # should be all true
bfl <- BamFileList(outpaths)
se2 <- summarizeOverlaps(
    flattenedAnnotation2, BamFileList(bfl), singleEnd=FALSE,
    fragments=TRUE, ignore.strand=TRUE )
assays(se2)$counts[1:4,]

## Exon-level read counting in parallel mode 
library(Rsamtools); library(GenomicAlignments); library(GenomicFeatures); library(BiocParallel)
outpaths <- list.files('results/hisat2_mapping', pattern='sorted.bam$', full.names=TRUE)
bfl <- BamFileList(outpaths, yieldSize = 50000, index = character())
multicoreParam <- MulticoreParam(workers = 4)
register(multicoreParam)
registered()
countExons <- bplapply(bfl, function(x) summarizeOverlaps(flattenedAnnotation2,
        x, ignore.strand = TRUE, singleEnd = FALSE, BPPARAM = multicoreParam))
countDFexons <- sapply(seq(along = countExons), function(x) assays(countExons[[x]])$counts)
rownames(countDFexons) <- names(rowRanges(countExons[[1]]))
colnames(countDFexons) <- names(bfl)
countDFexons[1:4,]
write.table(countDFexons, "results/countExons.xls", col.names = NA, quote = FALSE, sep = "\t")
countDF <- read.delim("results/countExons.xls", row.names = 1, check.names = FALSE)

## Building DEXSeqDataSet object

library(DEXSeq)
targets <- read.delim("targetsPE.txt", comment.char="#")
colData(se2)$condition <- factor(targets$Factor)
colData(se2)$libType <- factor(c("paired-end"))
dxd2 <- DEXSeqDataSetFromSE(se2, design= ~ sample + exon + condition:exon)
assays(dxd2)$counts[1:4,]
colData(dxd2)

## Note that the number of columns is now doubled from 18 to 36. The first 18 
## correspond to the number of reads mapping to each exon and the last 18 correspond 
## to the sum of the counts mapping to the rest of the exons from the same gene in 
## each sample. To check, run:
split( seq_len(ncol(dxd2)), colData(dxd2)$exon )

## To show only the counts for the actual sample, run: 
featureCounts(dxd2)[1:4,]

## The exon range information can be returned as follows:
rowRanges(dxd2)[1:4,]

## The experimental design table can be returned like this:
sampleAnnotation(dxd2)

## Normalisation

dxd2 <- estimateSizeFactors( dxd2 )

## Dispersion estimation
dxd2 <- estimateDispersions( dxd2 )
plotDispEsts( dxd2 )
```

## Testing for differential exon usage
dxd2 <- testForDEU( dxd2 )

##################################
## Continue in DEXseq vignette ###
##################################
## ...

