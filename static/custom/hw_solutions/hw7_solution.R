##########################################
## Solutions for HW7 - RNA-Seq Analysis ##
##########################################

#####################################################
## A. Unstranded and strand-specific read counting ## 
#####################################################

############
## Task 1 ##
############

## Generate count tables for exons by genes (eByg) ranges of the following three strand modes:

## 1. Unstranded
unstranded <- summarizeOverlaps(eByg, bfl, mode="Union", 
                                               ignore.strand=TRUE, 
                    					       inter.feature=FALSE, 
                                               singleEnd=FALSE)
unstranded <- assays(unstranded)$counts
unstranded[1:4,]

## 2. Positive strand 
stranded_pos <- summarizeOverlaps(eByg, bfl, mode="Union", 
                                               ignore.strand=FALSE, 
					                           inter.feature=FALSE, 
                                               singleEnd=FALSE)
stranded_pos <- assays(stranded_pos)$counts
stranded_pos[1:4,]

## 3. Negative strand 
stranded_neg <- summarizeOverlaps(eByg, bfl, mode="Union", 
                                               ignore.strand=FALSE, 
                                               preprocess.reads=invertStrand,
					                           inter.feature=FALSE, 
                                               singleEnd=FALSE)
stranded_neg <- assays(stranded_neg)$counts
stranded_neg[1:4,]

############
## Task 2 ## 
############
## Test whether the two strand-specific count tables sum up to very similar
## values as the unstranded count table. For paired end data one can only
## expect very similar but not identical results. This is the case since the
## strand spec results allow to count reads that can only be ambiguously
## assigned in the non-strand specific case. The following calculates for
## each strand specific result the averaged mapping frequency in percent 
## relative to the unstranded result. Since both percentage values (here perc_pos
## and perc_net) are close to 50%, the chosen RNA-Seq data set is very likely 
## unstranded. Other approximation approaches would provide correct answers under this 
## homework task as well.
perc_pos <- mean((stranded_pos / unstranded) * 100, na.rm=TRUE)
perc_neg <- mean((stranded_neg / unstranded) * 100, na.rm=TRUE)
perc_pos 
perc_neg
  

############
## Task 3 ## 
############
## Utility (biological relevance) of the different strand-specific counting modes used under Task 1:
   ## (i)   If an RNA-Seq experiment has been performed in a stranded manner then 
   ##       the read counting should be performed for the corresponding strand. 
   ## (ii)  The strand-specific information allows to unambiguously assign
   ##       reads to overlapping genes that are encoded on opposite strands.
   ## (iii) Read counting for the opposite strand can be useful to discover anti-sense
   ##       regulation events.
   ## (iv)  By performing the read counting in a stranded and unstranded manner, one can
   ##       determine whether an RNA-Seq experiment is strand specific or not.


##################################################
## B. Read counting for different feature types ## 
##################################################

############
## Task 4 ## 
############
## Compute strand-specific count tables for the positive (sense) strand of the following feature types:

## 1. Genes
gene_ranges <- genes(txdb)
gene_ranges_countDF <- summarizeOverlaps(gene_ranges, bfl, mode="Union", ignore.strand=FALSE, inter.feature=FALSE, singleEnd=FALSE)
assays(gene_ranges_countDF)$counts[1:4,]

## 2. Exons
exon_ranges <- exons(txdb)
exon_ranges_countDF <- summarizeOverlaps(exon_ranges, bfl, mode="Union", ignore.strand=FALSE, inter.feature=FALSE, singleEnd=FALSE)
assays(exon_ranges_countDF)$counts[1:4,]

## 3. Exons by genes
exonByg_ranges <- exonsBy(txdb, by=c("gene"))
exonByg_ranges_countDF <- summarizeOverlaps(exonByg_ranges, bfl, mode="Union", ignore.strand=FALSE, inter.feature=FALSE, singleEnd=FALSE)
assays(exonByg_ranges_countDF)$counts[1:4,]

## 4. Introns by transcripts
intron_ranges <- intronsByTranscript(txdb, use.names=TRUE) 
intron_ranges_countDF <- summarizeOverlaps(intron_ranges, bfl, mode="Union", ignore.strand=FALSE, inter.feature=FALSE, singleEnd=TRUE)
assays(intron_ranges_countDF)$counts[1:4,]

## 5. 5'-UTRs by transcripts
fiveUTR_ranges <- fiveUTRsByTranscript(txdb, use.names=TRUE)
fiveUTR_ranges_countDF <- summarizeOverlaps(fiveUTR_ranges, bfl, mode="Union", ignore.strand=FALSE, inter.feature=FALSE, singleEnd=TRUE)
assays(fiveUTR_ranges_countDF)$counts[1:4,]


#####################
## C. DEG analysis ##
#####################

############
## Task 5 ## 
############

## Perform the DEG analysis with edgeR with both unstranded and sense strand count tables

## DEG analysis with unstranded count table
library(edgeR)
cmp <- readComp(stepsWF(sal)[['hisat2_mapping']], format="matrix", delim="-")
edgeDF_unstranded <- run_edgeR(countDF=unstranded, targets=targetsWF(sal)[['hisat2_mapping']], cmp=cmp[[1]], independent=FALSE, mdsplot="")

## DEG analysis with count table for positive strand
edgeDF_pos <- run_edgeR(countDF=stranded_pos, targets=targetsWF(sal)[['hisat2_mapping']], cmp=cmp[[1]], independent=FALSE, mdsplot="")


## Compare the DEG results of the two methods in two separate 4-way Venn diagrams

## (1) 4-way Venn diagram for unstranded count table
DEG_list_unstranded <- filterDEGs(degDF=edgeDF_unstranded, filter=c(Fold=2, FDR=40), plot=FALSE)
vennsetup <- overLapper(DEG_list_unstranded$Up[6:9], type="vennsets")
vennsetdown <- overLapper(DEG_list_unstranded$Down[6:9], type="vennsets")
pdf("results/DEGcount_unstranded.pdf")
vennPlot(list(vennsetup, vennsetdown), mymain="", mysub="", colmode=2, ccol=c("blue", "red"))
dev.off()

## (2) 4-way Venn diagram for sense strand count table
DEG_list_pos <- filterDEGs(degDF=edgeDF_pos, filter=c(Fold=2, FDR=40), plot=FALSE)
vennsetup <- overLapper(DEG_list_pos$Up[6:9], type="vennsets")
vennsetdown <- overLapper(DEG_list_pos$Down[6:9], type="vennsets")
pdf("results/DEGcount_pos.pdf")
vennPlot(list(vennsetup, vennsetdown), mymain="", mysub="", colmode=2, ccol=c("blue", "red"))
dev.off()


