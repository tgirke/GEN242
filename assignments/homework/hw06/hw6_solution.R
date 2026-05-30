##########################################
## Solution for HW7 - RNA-Seq Analysis  ##
##########################################

library(Rsubread)
library(GenomicFeatures)
library(systemPipeR)
library(edgeR)

## ── Input setup ───────────────────────────────────────────────
## Run after completing the RNA-Seq workflow through HISAT2 mapping.
## The sal object and sorted BAM files must be available.
## The TxDb object (tair10.sqlite) is created from the GFF annotation
## during the workflow and loaded here with loadDb().

txdb <- loadDb("./data/tair10.sqlite")

## Exons grouped by gene — used as annotation for featureCounts
## in Parts A and B (equivalent to exonsBy(txdb, by="gene"))
eByg <- exonsBy(txdb, by = "gene")

## Helper: convert GRangesList/GRanges to featureCounts annotation data frame
grl_to_df <- function(grl) {
    gr <- unlist(grl)
    data.frame(GeneID = names(gr),
               Chr    = as.character(seqnames(gr)),
               Start  = start(gr),
               End    = end(gr),
               Strand = as.character(strand(gr)),
               stringsAsFactors = FALSE)
}

bam_files <- getColumn(sal, step = "hisat2_mapping",
                       "outfiles", column = "samtools_sort_bam")
## Alternative if sal is not available:
# bam_files <- list.files("results/hisat2_mapping/",
#                         pattern = "sorted.bam$", full.names = TRUE)

########################################################
## A. Unstranded and Strand-Specific Read Counting    ##
########################################################

############
## Task A1 ##
############

## 1. Unstranded — strandSpecific = 0
fc_unstranded <- featureCounts(
    files                  = bam_files,
    annot.ext              = grl_to_df(eByg),
    isGTFAnnotationFile    = FALSE,
    useMetaFeatures        = TRUE,
    strandSpecific         = 0,
    isPairedEnd            = TRUE,
    countMultiMappingReads = FALSE,
    nthreads               = 4
)
unstranded <- fc_unstranded$counts
unstranded[1:4, ]

## 2. Sense (positive) strand — strandSpecific = 1
fc_stranded_pos <- featureCounts(
    files                  = bam_files,
    annot.ext              = grl_to_df(eByg),
    isGTFAnnotationFile    = FALSE,
    useMetaFeatures        = TRUE,
    strandSpecific         = 1,
    isPairedEnd            = TRUE,
    countMultiMappingReads = FALSE,
    nthreads               = 4
)
stranded_pos <- fc_stranded_pos$counts
stranded_pos[1:4, ]

## 3. Antisense (negative) strand — strandSpecific = 2
fc_stranded_neg <- featureCounts(
    files                  = bam_files,
    annot.ext              = grl_to_df(eByg),
    isGTFAnnotationFile    = FALSE,
    useMetaFeatures        = TRUE,
    strandSpecific         = 2,
    isPairedEnd            = TRUE,
    countMultiMappingReads = FALSE,
    nthreads               = 4
)
stranded_neg <- fc_stranded_neg$counts
stranded_neg[1:4, ]

############
## Task A2 ##
############
## Test whether the two strand-specific count tables sum to very similar
## values as the unstranded count table. For paired-end data one can only
## expect very similar but not identical results. The following calculates
## for each strand-specific result the averaged mapping frequency in percent
## relative to the unstranded result.

perc_pos <- mean((stranded_pos / unstranded) * 100, na.rm = TRUE)
perc_neg <- mean((stranded_neg / unstranded) * 100, na.rm = TRUE)
perc_pos
perc_neg
## Since both perc_pos and perc_neg are close to 50%, the RNA-Seq data set
## used here is very likely unstranded.

############
## Task A3 ##
############
## Biological utility of the different strand-specific counting modes:
##
## 1. When to use strand-specific counting:
##    If an RNA-Seq library has been prepared in a strand-specific manner
##    (e.g. using dUTP or ligation-based strand-specific library prep kits),
##    read counting should be performed using the corresponding strandSpecific
##    mode (1 or 2) to correctly assign reads to genes. Using unstranded
##    counting on stranded data would undercount reads and reduce sensitivity.
##
## 2. Overlapping genes on opposite strands:
##    Strand-specific counting allows unambiguous assignment of reads to
##    overlapping genes encoded on opposite strands. Without strand information,
##    reads from two overlapping antisense genes cannot be reliably assigned
##    to one or the other, leading to ambiguous or discarded counts.
##
## 3. Antisense strand counting:
##    Counting reads on the antisense strand (strandSpecific = 2) can reveal
##    natural antisense transcripts (NATs) and other antisense regulation
##    events such as gene silencing via antisense non-coding RNAs. These
##    are biologically important but invisible in unstranded or sense-only
##    counting.
##
## 4. Determining library strandedness:
##    By comparing the sense and antisense strand-specific count tables with
##    the unstranded count table, one can determine whether an RNA-Seq
##    experiment used a strand-specific protocol. If both perc_pos and
##    perc_neg are close to ~50%, the library is unstranded. If one value
##    dominates (e.g. ~100% sense, ~0% antisense), the library is
##    strand-specific and the appropriate strandSpecific mode should be
##    used for downstream analysis.

##################################################
## B. Read Counting for Different Feature Types ##
##################################################

############
## Task B  ##
############

## 1. Genes — eByg = exonsBy(txdb, by="gene"), loaded from tair10.sqlite
fc_genes <- featureCounts(
    files                  = bam_files,
    annot.ext              = grl_to_df(eByg),
    isGTFAnnotationFile    = FALSE,
    useMetaFeatures        = TRUE,
    strandSpecific         = 1,
    isPairedEnd            = TRUE,
    countMultiMappingReads = FALSE,
    nthreads               = 4
)
gene_ranges_countDF <- fc_genes$counts
gene_ranges_countDF[1:4, ]

## 2. Exons — individual exon level using exons(txdb)
exon_ranges <- exons(txdb)
fc_exons <- featureCounts(
    files                  = bam_files,
    annot.ext              = grl_to_df(exon_ranges),
    isGTFAnnotationFile    = FALSE,
    useMetaFeatures        = FALSE,
    strandSpecific         = 1,
    isPairedEnd            = TRUE,
    countMultiMappingReads = FALSE,
    nthreads               = 4
)
exon_ranges_countDF <- fc_exons$counts
exon_ranges_countDF[1:4, ]

## 3. Exons by genes — eByg is directly equivalent to exonsBy(txdb, by="gene")
##    so we reuse the gene-level result
exonByg_ranges_countDF <- gene_ranges_countDF
exonByg_ranges_countDF[1:4, ]

## 4. Introns by transcripts — derive from txdb, pass GRangesList directly
intron_ranges <- intronsByTranscript(txdb, use.names = TRUE)
fc_introns <- featureCounts(
    files                  = bam_files,
    annot.ext              = grl_to_df(intron_ranges),
    isGTFAnnotationFile    = FALSE,
    strandSpecific         = 1,
    isPairedEnd            = FALSE,
    countMultiMappingReads = FALSE,
    nthreads               = 4
)
intron_ranges_countDF <- fc_introns$counts
intron_ranges_countDF[1:4, ]

## 5. 5'-UTRs by transcripts — same approach as introns
fiveUTR_ranges <- fiveUTRsByTranscript(txdb, use.names = TRUE)
fc_fiveUTR <- featureCounts(
    files                  = bam_files,
    annot.ext              = grl_to_df(fiveUTR_ranges),
    isGTFAnnotationFile    = FALSE,
    strandSpecific         = 1,
    isPairedEnd            = FALSE,
    countMultiMappingReads = FALSE,
    nthreads               = 4
)
fiveUTR_ranges_countDF <- fc_fiveUTR$counts
fiveUTR_ranges_countDF[1:4, ]

#####################
## C. DEG Analysis ##
#####################

############
## Task C  ##
############

## Retrieve sample comparisons from the workflow sal object
cmp <- readComp(stepsWF(sal)[["hisat2_mapping"]],
                format = "matrix", delim = "-")
## Alternative if sal is not available:
# cmp <- readComp(file = "targetsPE.txt", format = "matrix", delim = "-")

## DEG analysis with unstranded count table
edgeDF_unstranded <- run_edgeR(
    countDF     = unstranded,
    targets     = targetsWF(sal)[["hisat2_mapping"]],
    cmp         = cmp[[1]],
    independent = FALSE,
    mdsplot     = ""
)

## DEG analysis with sense strand count table
edgeDF_pos <- run_edgeR(
    countDF     = stranded_pos,
    targets     = targetsWF(sal)[["hisat2_mapping"]],
    cmp         = cmp[[1]],
    independent = FALSE,
    mdsplot     = ""
)

## (1) 4-way Venn diagram for unstranded count table
DEG_list_unstranded <- filterDEGs(
    degDF  = edgeDF_unstranded,
    filter = c(Fold = 2, FDR = 40),
    plot   = FALSE
)
vennsetup   <- overLapper(DEG_list_unstranded$Up[6:9],   type = "vennsets")
vennsetdown <- overLapper(DEG_list_unstranded$Down[6:9], type = "vennsets")
pdf("results/DEGcount_unstranded.pdf")
vennPlot(list(vennsetup, vennsetdown),
         mymain = "", mysub = "", colmode = 2,
         ccol = c("blue", "red"))
dev.off()

## (2) 4-way Venn diagram for sense strand count table
DEG_list_pos <- filterDEGs(
    degDF  = edgeDF_pos,
    filter = c(Fold = 2, FDR = 40),
    plot   = FALSE
)
vennsetup   <- overLapper(DEG_list_pos$Up[6:9],   type = "vennsets")
vennsetdown <- overLapper(DEG_list_pos$Down[6:9], type = "vennsets")
pdf("results/DEGcount_pos.pdf")
vennPlot(list(vennsetup, vennsetdown),
         mymain = "", mysub = "", colmode = 2,
         ccol = c("blue", "red"))
dev.off()

###################################################################
## D. Bonus Challenge — edgeR Native Mode, Time Course Design    ##
###################################################################

## The Howard et al. (2013) experiment is a 2-factor time course:
##   treatment: mock (M) vs Pseudomonas-infected (A)
##   time: 1h, 2h, 12h post-inoculation
## An interaction model ~ treatment + time + treatment:time
## captures whether the treatment effect changes over time —
## the biologically meaningful question in a time course experiment.

## Step 1 — Prepare sample metadata
targets           <- read.delim("targetsPE.txt", comment.char = "#")
targets$treatment <- ifelse(grepl("^M", targets$Factor), "mock", "infected")
targets$time      <- gsub("^[MA]", "", targets$Factor)
targets$time      <- factor(targets$time,      levels = c("1", "2", "12"))
targets$treatment <- factor(targets$treatment, levels = c("mock", "infected"))
targets[, c("SampleName", "Factor", "treatment", "time")]

## Step 2 — Create DGEList and filter low-count genes
countDF <- unstranded[, targets$SampleName]
countDF[is.na(countDF)] <- 0
dge     <- DGEList(counts = countDF, group = targets$Factor)

## Keep genes with CPM >= 1 in at least 2 samples
keep <- rowSums(cpm(dge) >= 1) >= 2
dge  <- dge[keep, , keep.lib.sizes = FALSE]
message(sum(keep), " genes retained after low-count filtering")

## TMM normalization
dge <- calcNormFactors(dge, method = "TMM")

## Step 3 — Build interaction design matrix
design_tc <- model.matrix(~ treatment + time + treatment:time,
                          data = targets)
colnames(design_tc)
## Intercept, treatmentinfected, time2, time12,
## treatmentinfected:time2, treatmentinfected:time12

## Step 4 — Estimate dispersion and fit QL GLM
dge    <- estimateDisp(dge, design_tc, robust = TRUE)
plotBCV(dge)
fit_tc <- glmQLFit(dge, design_tc, robust = TRUE)

## Step 5 — Test interaction terms (treatment effect changes over time)
interaction_cols <- grep("treatmentinfected:time", colnames(design_tc))
qlf_interaction  <- glmQLFTest(fit_tc, coef = interaction_cols)

topTags(qlf_interaction, n = 20)
summary(decideTests(qlf_interaction, p.value = 0.05))

tc_results <- as.data.frame(topTags(qlf_interaction, n = Inf))
write.csv(tc_results, "results/edgeR_timecourse_interaction.csv")

## Step 6 — Comparison with pairwise result from Part C
##
## The interaction model identifies genes whose response to Pseudomonas
## infection changes significantly over time — e.g. genes strongly induced
## at 12h but not at 1h or 2h. This is the core biological question of a
## time course experiment and tests all samples jointly in a single model.
##
## The pairwise approach in run_edgeR (Part C) compares mock vs infected
## at each time point independently. It asks "which genes differ between
## conditions at each time point?" but does not model the temporal trajectory.
##
## The interaction model is generally more powerful for detecting dynamic
## treatment responses because it borrows information across time points.
## It typically returns a different — often larger — biologically relevant
## DEG set. The two approaches answer complementary questions:
##   - Pairwise: "Is this gene DE at time X?"
##   - Interaction: "Does this gene's response to infection change over time?"
## Choice of design should always be driven by the biological question.
