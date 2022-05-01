## ----setup_dir, echo=FALSE, include=FALSE, message=FALSE, warning=FALSE----
unlink(".SPRproject/", recursive = TRUE)


## pre code {

## white-space: pre !important;

## overflow-x: scroll !important;

## word-break: keep-all !important;

## word-wrap: initial !important;

## }


## ----style, echo = FALSE, results = 'asis'----------------
BiocStyle::markdown()
options(width=60, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")), 
    tidy.opts=list(width.cutoff=60), tidy=TRUE)


## ----setup, echo=FALSE, message=FALSE, wwarning=FALSE, eval=FALSE----
## suppressPackageStartupMessages({
##     library(systemPipeR)
## })


## srun --x11 --partition=gen242 --mem=20gb --cpus-per-task 8 --ntasks 1 --time 20:00:00 --pty bash -l

## module unload R; module load load R/4.1.2


## $ Rscript -e "systemPipeRdata::genWorkenvir(workflow='chipseq')"

## $ cd chipseq


## ----gen_workflow_envir, eval=FALSE-----------------------
## library(systemPipeRdata)
## genWorkenvir(workflow="chipseq")
## setwd("chipseq")


## ----load_systempiper, eval=TRUE, message=FALSE, warning=FALSE----
library(systemPipeR)


## ----load_custom_fct, eval=FALSE, message=FALSE-----------
## source("challengeProject_Fct.R")


## ----load_targets_file, eval=TRUE-------------------------
targetspath <- "targets_chipseq.txt"
targets <- read.delim(targetspath, comment.char = "#")
knitr::kable(targets)


## ----create_workflow, message=FALSE, eval=TRUE------------
library(systemPipeR)
sal <- SPRproject()
sal


## ----importWF, message=FALSE, eval=FALSE------------------
## sal <- importWF(sal, file_path = "systemPipeChIPseq.Rmd") ## Import all the Workflow steps
## sal <- runWF(sal) # Runs workflow


## ----load_SPR, message=FALSE, eval=TRUE, spr=TRUE---------
appendStep(sal) <- LineWise(code = {
                library(systemPipeR)
                }, step_name = "load_SPR")


## ----preprocessing, message=FALSE, eval=TRUE, spr=TRUE----
appendStep(sal) <- SYSargsList(
    step_name = "preprocessing",
    targets = "targets_chipseq.txt", dir = TRUE,
    wf_file = "preprocessReads/preprocessReads-se.cwl",
    input_file = "preprocessReads/preprocessReads-se.yml",
    dir_path = "param/cwl",
    inputvars = c(
        FileName = "_FASTQ_PATH1_",
        SampleName = "_SampleName_"
    ),
    dependency = c("load_SPR"))


## ----editing_preprocessing, message=FALSE, eval=TRUE------
yamlinput(sal, "preprocessing")$Fct


## ----fastq_report, eval=TRUE, message=FALSE, spr=TRUE-----
appendStep(sal) <- LineWise(
    code = {
        fastq <- getColumn(sal, step = "preprocessing", "targetsWF", column = 1)
        fqlist <- seeFastq(fastq = fastq, batchsize = 10000, klength = 8)
        pdf("./results/fastqReport.pdf", height = 18, width = 4*length(fqlist))
        seeFastqPlot(fqlist)
        dev.off()
    }, step_name = "fastq_report", 
    dependency = "preprocessing")


## ----bowtie2_index, eval=TRUE, spr=TRUE-------------------
appendStep(sal) <- SYSargsList(
    step_name = "bowtie2_index",
    dir = FALSE, targets = NULL,
    wf_file = "bowtie2/bowtie2-index.cwl",
    input_file = "bowtie2/bowtie2-index.yml",
    dir_path = "param/cwl",
    inputvars = NULL,
    dependency = c("preprocessing")
)


## ----bowtie2_alignment, eval=TRUE, spr=TRUE---------------
appendStep(sal) <- SYSargsList(
    step_name = "bowtie2_alignment",
    dir = TRUE,
    targets = "targets_chipseq.txt",
    wf_file = "workflow-bowtie2/workflow_bowtie2-se.cwl",
    input_file = "workflow-bowtie2/workflow_bowtie2-se.yml",
    dir_path = "param/cwl",
    inputvars = c(
        FileName = "_FASTQ_PATH1_",
        SampleName = "_SampleName_"
    ),
    dependency = c("bowtie2_index")
)


## ----bowtie2_alignment_check, eval=TRUE-------------------
cmdlist(sal, step="bowtie2_alignment", targets=1)


## ----align_stats, eval=TRUE, spr=TRUE---------------------
appendStep(sal) <- LineWise(
    code = {
        fqpaths <- getColumn(sal, step = "bowtie2_alignment", "targetsWF", column = "FileName")
        bampaths <- getColumn(sal, step = "bowtie2_alignment", "outfiles", column = "samtools_sort_bam")
        read_statsDF <- alignStats(args = bampaths, fqpaths = fqpaths, pairEnd = TRUE)
        write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
        }, 
    step_name = "align_stats", 
    dependency = "bowtie2_alignment")


## ----bam_IGV, eval=TRUE, spr=TRUE-------------------------
appendStep(sal) <- LineWise(
    code = {
        bampaths <- getColumn(sal, step = "bowtie2_alignment", "outfiles", 
                              column = "samtools_sort_bam")
        symLink2bam(
            sysargs = bampaths, htmldir = c("~/.html/", "somedir/"),
            urlbase="http://cluster.hpcc.ucr.edu/~<username>/",
            urlfile = "./results/IGVurl.txt")
    },
    step_name = "bam_IGV",
    dependency = "bowtie2_alignment",
    run_step = "optional"
)


## ----merge_bams, eval=TRUE, spr=TRUE----------------------
appendStep(sal) <- LineWise(
    code = {
        bampaths <- getColumn(sal, step = "bowtie2_alignment", "outfiles", column = "samtools_sort_bam")
        merge_bams <- mergeBamByFactor(args=bampaths, targetsDF = targetsWF(sal)[["bowtie2_alignment"]], overwrite=TRUE)
        updateColumn(sal, step = "merge_bams", position = "targetsWF") <- merge_bams
        writeTargets(sal, step = "merge_bams", file = "targets_merge_bams.txt", overwrite = TRUE)
    },
    step_name = "merge_bams",
    dependency = "bowtie2_alignment"
)


## ----writeTargetsRef, eval=TRUE, spr=TRUE-----------------
appendStep(sal) <- LineWise(
    code = {
        writeTargetsRef(infile = "targets_merge_bams.txt", 
                outfile = "targets_bam_ref.txt", silent = FALSE, overwrite = TRUE)
    },
    step_name = "writeTargetsRef",
    dependency = "merge_bams"
)


## ----call_peaks_macs_withref, eval=TRUE, spr=TRUE---------
appendStep(sal) <- SYSargsList(
    step_name = "call_peaks_macs_withref",
    targets = "targets_bam_ref.txt",
    wf_file = "MACS2/macs2-input.cwl",
    input_file = "MACS2/macs2-input.yml",
    dir_path = "param/cwl",
    inputvars = c(
        FileName1 = "_FASTQ_PATH1_",
        FileName2 = "_FASTQ_PATH2_",
        SampleReference = "_SampleName_"
    ),
    id = "SampleReference", 
    dependency = c("writeTargetsRef")
)


## ----annotation_ChIPseeker, eval=TRUE, spr=TRUE-----------
appendStep(sal) <- LineWise(
    code = {
        library(ChIPseeker); library(GenomicFeatures)
        peaks_files <- getColumn(sal, step = "call_peaks_macs_withref", "outfiles", column = "peaks_xls")
        txdb <- suppressWarnings(makeTxDbFromGFF(file="data/tair10.gff", format="gff", dataSource="TAIR", 
                                organism="Arabidopsis thaliana"))
        for(i in seq(along=peaks_files)) {
            peakAnno <- annotatePeak(peaks_files[i], TxDb=txdb, verbose=FALSE)
            df <- as.data.frame(peakAnno)
            outpaths <- paste0("./results/", names(peaks_files), "_ChIPseeker_annotated.xls")
            names(outpaths) <- names(peaks_files)
            write.table(df, outpaths[i], quote=FALSE, row.names=FALSE, sep="\t")
        }
        updateColumn(sal, step = "annotation_ChIPseeker", position = "outfiles") <- data.frame(outpaths)
    },
    step_name = "annotation_ChIPseeker",
    dependency = "call_peaks_macs_withref"
)



## ----ChIPseeker_plots, eval=TRUE, spr=TRUE----------------
appendStep(sal) <- LineWise(
    code = {
        peaks_files <- getColumn(sal, step = "call_peaks_macs_withref", "outfiles", column = "peaks_xls")
        peak <- readPeakFile(peaks_files[1])
        pdf("results/peakscoverage.pdf")
        covplot(peak, weightCol="X.log10.pvalue.")
        dev.off()
        pdf("results/peaksHeatmap.pdf")
        peakHeatmap(peaks_files[1], TxDb=txdb, upstream=1000, downstream=1000, 
                    color="red")
        dev.off()
        pdf("results/peaksProfile.pdf")
        plotAvgProf2(peaks_files[1], TxDb=txdb, upstream=1000, downstream=1000, 
                     xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency", 
                     conf=0.05)
        dev.off()
    },
    step_name = "ChIPseeker_plots",
    dependency = "annotation_ChIPseeker"
)


## ----count_peak_ranges, eval=TRUE, spr=TRUE---------------
appendStep(sal) <- LineWise(
    code = {
        library(GenomicRanges)
        bam_files <- getColumn(sal, step = "bowtie2_alignment", "outfiles", column = "samtools_sort_bam")
        args <- getColumn(sal, step = "call_peaks_macs_withref", "outfiles", column = "peaks_xls")
        outfiles <- paste0("./results/", names(args), "_countDF.xls")
        bfl <- BamFileList(bam_files, yieldSize=50000, index=character())
        countDFnames <- countRangeset(bfl, args, outfiles, mode="Union", ignore.strand=TRUE)
        updateColumn(sal, step = "count_peak_ranges", position = "outfiles") <- data.frame(countDFnames)
    },
    step_name = "count_peak_ranges",
    dependency = "call_peaks_macs_withref",
)


## ----show_counts_table, eval=TRUE-------------------------
countDF <- read.delim("results/AP1_1_countDF.xls")[1:200,]
colnames(countDF)[1] <- "PeakIDs"
library(DT)
datatable(countDF)


## ----diff_bind_analysis, eval=TRUE, spr=TRUE--------------
appendStep(sal) <- LineWise(
    code = {
        countDF_files <- getColumn(sal, step = "count_peak_ranges", "outfiles")
        outfiles <- paste0("./results/", names(countDF_files), "_peaks_edgeR.xls")
        names(outfiles) <- names(countDF_files)
        cmp <- readComp(file =stepsWF(sal)[["bowtie2_alignment"]],
                        format="matrix") 
        dbrlist <- runDiff(args=countDF_files, outfiles = outfiles, diffFct=run_edgeR, 
                           targets=targetsWF(sal)[["bowtie2_alignment"]], cmp=cmp[[1]], 
                           independent=TRUE, dbrfilter=c(Fold=2, FDR=1))
    },
    step_name = "diff_bind_analysis",
    dependency = "count_peak_ranges",
)


## ----go_enrich, eval=TRUE, spr=TRUE-----------------------
appendStep(sal) <- LineWise(
    code = {
        annofiles <- getColumn(sal, step = "annotation_ChIPseeker", "outfiles")
        gene_ids <- sapply(annofiles, 
                           function(x) unique(as.character
                                              (read.delim(x)[,"geneId"])), simplify=FALSE)
        load("data/GO/catdb.RData")
        BatchResult <- GOCluster_Report(catdb=catdb, setlist=gene_ids, method="all", 
                                        id_type="gene", CLSZ=2, cutoff=0.9, 
                                        gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)
        write.table(BatchResult, "results/GOBatchAll.xls", quote=FALSE, row.names=FALSE, sep="\t")
    },
    step_name = "go_enrich",
    dependency = "annotation_ChIPseeker",
)


## ----show_GO_table, eval=TRUE-----------------------------
BatchResult <- read.delim("results/GOBatchAll.xls")[1:50,]
library(DT)
datatable(BatchResult[,-10], options = list(scrollX = TRUE, autoWidth = TRUE))


## ----parse_peak_sequences, eval=TRUE, spr=TRUE------------
appendStep(sal) <- LineWise(
    code = {
        library(Biostrings); library(seqLogo); library(BCRANK)
        rangefiles <- getColumn(sal, step = "call_peaks_macs_withref", "outfiles")
        for(i in seq(along=rangefiles)) {
            df <- read.delim(rangefiles[i], comment="#")
            peaks <- as(df, "GRanges")
            names(peaks) <- paste0(as.character(seqnames(peaks)), "_", start(peaks),
                                   "-", end(peaks))
            peaks <- peaks[order(values(peaks)$X.log10.pvalue., decreasing=TRUE)]
            pseq <- getSeq(FaFile("./data/tair10.fasta"), peaks)
            names(pseq) <- names(peaks)
            writeXStringSet(pseq, paste0(rangefiles[i], ".fasta"))
            }
        },
    step_name = "parse_peak_sequences",
    dependency = "call_peaks_macs_withref"
)


## ----bcrank_enrich, eval=TRUE, spr=TRUE-------------------
appendStep(sal) <- LineWise(
    code = {
        library(Biostrings); library(seqLogo); library(BCRANK)
        rangefiles <- getColumn(sal, step = "call_peaks_macs_withref", "outfiles")
        set.seed(0)
        BCRANKout <- bcrank(paste0(rangefiles[1], ".fasta"), restarts=25,
                            use.P1=TRUE, use.P2=TRUE)
        toptable(BCRANKout)
        topMotif <- toptable(BCRANKout, 1)
        weightMatrix <- pwm(topMotif, normalize = FALSE)
        weightMatrixNormalized <- pwm(topMotif, normalize = TRUE)
        pdf("results/seqlogo.pdf")
        seqLogo(weightMatrixNormalized)
        dev.off()
        },
    step_name = "bcrank_enrich",
    dependency = "call_peaks_macs_withref"
)


## ----sessionInfo, eval=TRUE, spr=TRUE---------------------
appendStep(sal) <- LineWise(
    code = {
        sessionInfo()
        },
    step_name = "sessionInfo",
    dependency = "bcrank_enrich"
)


## ----runWF, eval=FALSE------------------------------------
## sal <- runWF(sal)


## ----runWF_cluster, eval=FALSE----------------------------
## resources <- list(conffile=".batchtools.conf.R",
##                   template="batchtools.slurm.tmpl",
##                   Njobs=18,
##                   walltime=120, ## minutes
##                   ntasks=1,
##                   ncpus=4,
##                   memory=1024, ## Mb
##                   partition = "short"
##                   )
## sal <- addResources(sal, c("bowtie2_alignment"), resources = resources)
## sal <- runWF(sal)


## ----plotWF, eval=TRUE------------------------------------
plotWF(sal)


## ----statusWF, eval=FALSE---------------------------------
## sal
## statusWF(sal)


## ----logsWF, eval=FALSE-----------------------------------
## sal <- renderLogs(sal)


## ----setup_clean, echo=FALSE, include=FALSE---------------
unlink(".SPRproject/", recursive = TRUE)

