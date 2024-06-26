#####################################
## Run DEXseq on Howard et al Data ##
#####################################
## Helper code (May 25, 2024)

## Create flattened exon ranges 
library(GenomicFeatures)  # for the exonicParts() function
library(rtracklayer)        # for the makeTxDbFromGFF() function
txdb <- loadDb("./data/tair10.sqlite")
flattenedAnnotation2 <- exonicParts(txdb, linked.to.single.gene.only=TRUE)
names(flattenedAnnotation2) <- sprintf("%s:E%0.3d", flattenedAnnotation2$gene_id, flattenedAnnotation2$exonic_part)

## Exon-level read counting in serial mode 
library(Rsamtools)
library(GenomicAlignments)
targets <- read.delim("targetsPE.txt", comment.char="#")
outpaths <- list.files('results/hisat2_mapping', pattern='sorted.bam$', full.names=TRUE)
outpaths <- setNames(outpaths, gsub("^.*mapping/(.*)\\.sorted.bam", "\\1", outpaths))
outpaths <- outpaths[targets$SampleName] # This ensures proper colmun order
file.exists(outpaths) # should be all true
bfl <- BamFileList(outpaths)
se2 <- summarizeOverlaps(
    flattenedAnnotation2, BamFileList(bfl), singleEnd=FALSE,
    fragments=TRUE, ignore.strand=TRUE)
assays(se2)$counts[1:4,]

## Building DEXSeqDataSet object
library(DEXSeq)
colData(se2)$condition <- factor(targets$Factor)
colData(se2)$libType <- factor(c("paired-end"))
dxd2 <- DEXSeqDataSetFromSE(se2, design= ~ sample + exon + condition:exon)
assays(dxd2)$counts[1:4,]
colData(dxd2)

## Exon-level read counting in parallel mode 
library(Rsamtools); library(GenomicAlignments); library(GenomicFeatures); library(BiocParallel)
targets <- read.delim("targetsPE.txt", comment.char="#")
outpaths <- list.files('results/hisat2_mapping', pattern='sorted.bam$', full.names=TRUE)
outpaths <- setNames(outpaths, gsub("^.*mapping/(.*)\\.sorted.bam", "\\1", outpaths))
outpaths <- outpaths[targets$SampleName] # This ensures proper colmun order
file.exists(outpaths) # should be all true
bfl <- BamFileList(outpaths, yieldSize = 50000, index = character())
multicoreParam <- MulticoreParam(workers = 8)
register(multicoreParam)
registered()
countExons <- bplapply(bfl, function(x) summarizeOverlaps(flattenedAnnotation2,
        x, fragments=TRUE, ignore.strand = TRUE, singleEnd = FALSE, BPPARAM = multicoreParam))
countDFexons <- sapply(seq(along = countExons), function(x) assays(countExons[[x]])$counts)
rownames(countDFexons) <- names(rowRanges(countExons[[1]]))
colnames(countDFexons) <- names(bfl)
countDFexons[1:4,]
write.table(countDFexons, "results/countExons.xls", col.names = NA, quote = FALSE, sep = "\t")

## Building DEXSeqDataSet object from imported count file (coundDF)
library(DEXSeq)
countDF <- read.delim("results/countExons.xls", row.names = 1, check.names = FALSE)
countDF <- SummarizedExperiment(assays=list(counts=countDF), rowRanges=flattenedAnnotation2)
assays(countDF)$counts[1:4,]
colData(countDF)
targets <- read.delim("targetsPE.txt", comment.char="#")
colData(countDF)$condition <- factor(targets$Factor)
colData(countDF)$libType <- factor(c("paired-end"))
dxd3 <- DEXSeqDataSetFromSE(countDF, design= ~ sample + exon + condition:exon)
assays(dxd3)$counts[1:4,]
colData(dxd3)

## Note that the number of columns is now doubled from 18 to 36. The first 18 
## correspond to the number of reads mapping to each exon and the last 18 correspond 
## to the sum of the counts mapping to the rest of the exons from the same gene in 
## each sample. To check, run:
split( seq_len(ncol(dxd3)), colData(dxd3)$exon )

## To show only the counts for the actual sample, run: 
featureCounts(dxd3)[1:4,]

## The exon range information can be returned as follows:
rowRanges(dxd3)[1:4,]

## The experimental design table can be returned like this:
sampleAnnotation(dxd3)

## Normalisation
dxd3 <- estimateSizeFactors( dxd3 )
## Dispersion estimation
dxd3 <- estimateDispersions( dxd3 )

## Plot fit diagnostics
png("results/dexseq_plotdispests.png") # Writes plot to png file
plotDispEsts( dxd3 )
dev.off()

## Testing for differential exon usage
dxd3 <- testForDEU( dxd3 )
dxd3 <- estimateExonFoldChanges( dxd3, fitExpToVar="condition")
dxr4 <- DEXSeqResults( dxd3 )
dxr4
mcols(dxr4)$description

table ( dxr4$padj < 0.1 )
table ( tapply( dxr4$padj < 0.1, dxr4$groupID, any ) )

## Plot mean expressions versus log2 fold changes
png("results/dexseq_plotma.png") # Writes plot to png file
plotMA( dxr4, cex=0.8 )
dev.off()

## Plot exon/transcript abundance for slected genes

dxr4_filtered <- dxr4[dxr4$padj < 0.1 & !is.na(dxr4$padj),]
dxr4_filtered <- dxr4_filtered[order(dxr4_filtered$padj),]
dxr4_filtered[1:4, 1:8] # sorted by padj. Useful for gene selection
png("results/dexseq_AT2G35840.png")
plotDEXSeq( dxr4, "AT2G35840", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()

##################################
## Continue in DEXseq vignette ###
##################################

