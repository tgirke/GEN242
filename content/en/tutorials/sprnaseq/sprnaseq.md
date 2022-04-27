---
title: "RNA-Seq Workflow Template" 
author: "Author: First Last"
date: "Last update: 27 April, 2022" 
output:
  BiocStyle::html_document:
    toc_float: true
    code_folding: show
package: systemPipeR
vignette: |
  %\VignetteIndexEntry{WF: RNA-Seq Workflow Template}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
fontsize: 14pt
bibliography: bibtex.bib
weight: 8
type: docs
---

<script src="/rmarkdown-libs/htmlwidgets/htmlwidgets.js"></script>

<link href="/rmarkdown-libs/datatables-css/datatables-crosstalk.css" rel="stylesheet" />

<script src="/rmarkdown-libs/datatables-binding/datatables.js"></script>

<script src="/rmarkdown-libs/jquery/jquery-3.6.0.min.js"></script>

<link href="/rmarkdown-libs/dt-core/css/jquery.dataTables.min.css" rel="stylesheet" />
<link href="/rmarkdown-libs/dt-core/css/jquery.dataTables.extra.css" rel="stylesheet" />

<script src="/rmarkdown-libs/dt-core/js/jquery.dataTables.min.js"></script>

<link href="/rmarkdown-libs/crosstalk/css/crosstalk.min.css" rel="stylesheet" />

<script src="/rmarkdown-libs/crosstalk/js/crosstalk.min.js"></script>

<script src="/rmarkdown-libs/htmlwidgets/htmlwidgets.js"></script>

<link href="/rmarkdown-libs/datatables-css/datatables-crosstalk.css" rel="stylesheet" />

<script src="/rmarkdown-libs/datatables-binding/datatables.js"></script>

<script src="/rmarkdown-libs/jquery/jquery-3.6.0.min.js"></script>

<link href="/rmarkdown-libs/dt-core/css/jquery.dataTables.min.css" rel="stylesheet" />
<link href="/rmarkdown-libs/dt-core/css/jquery.dataTables.extra.css" rel="stylesheet" />

<script src="/rmarkdown-libs/dt-core/js/jquery.dataTables.min.js"></script>

<link href="/rmarkdown-libs/crosstalk/css/crosstalk.min.css" rel="stylesheet" />

<script src="/rmarkdown-libs/crosstalk/js/crosstalk.min.js"></script>

<script src="/rmarkdown-libs/htmlwidgets/htmlwidgets.js"></script>

<link href="/rmarkdown-libs/vizjs/plotwf.css" rel="stylesheet" />

<script src="/rmarkdown-libs/vizjs/viz.js"></script>

<script src="/rmarkdown-libs/vizjs/full.render.js"></script>

<script src="/rmarkdown-libs/dom_to_image/dom_to_image.js"></script>

<link id="plotwf_legend-1-attachment" rel="attachment" href="sprnaseq_files/plotwf_legend/plotwf_legend.svg"/>

<script src="/rmarkdown-libs/plotwf-binding/plotwf.js"></script>

<style type="text/css">
pre code {
white-space: pre !important;
overflow-x: scroll !important;
word-break: keep-all !important;
word-wrap: initial !important;
}
</style>

<!--
Rscript -e "rmarkdown::render('sprnaseq.Rmd', c('BiocStyle::html_document'), clean=F); knitr::knit('sprnaseq.Rmd', tangle=TRUE)"
-->

<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  document.querySelector("h1").className = "title";
});
</script>

<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  var links = document.links;  
  for (var i = 0, linksLength = links.length; i < linksLength; i++)
    if (links[i].hostname != window.location.hostname)
      links[i].target = '_blank';
});
</script>

<div style="text-align: right">

Source code download:    
\[ [.Rmd](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/sprnaseq/sprnaseq.Rmd) \]    
\[ [.R](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/sprnaseq/sprnaseq.R) \]

</div>

## Introduction

This report describes the analysis of the RNA-Seq data set from
Howard et al (2013). The corresponding FASTQ files were downloaded from
GEO (Accession: [SRP010938](http://www.ncbi.nlm.nih.gov/sra/?term=SRP010938)).
This data set contains 18 paired-end (PE) read sets from *Arabidposis thaliana*.
The details about all download steps are provided [here](https://girke.bioinformatics.ucr.edu/GEN242/assignments/projects/project_data/).

Users want to provide here additional background information about the design of their
RNA-Seq project.

### Experimental design

Typically, users want to specify here all information relevant for the
analysis of their NGS study. This includes detailed descriptions of
FASTQ files, experimental design, reference genome, gene annotations,
etc.

### Workflow environment

<font color="red">NOTE: this section</font> describes how to set up the proper
environment (directory structure) for running `systemPipeR` workflows. After
mastering this task the workflow run instructions <font color="red">can be
deleted</font> since they are not expected to be included in a final HTML/PDF
report of a workflow.

1.  If a remote system or cluster is used, then users need to log in to the
    remote system first. The following applies to an HPC cluster (*e.g.* HPCC
    cluster).
    
    A terminal application needs to be used to log in to a user’s cluster account. Next, one
    can open an interactive session on a computer node with `srun`. More details about
    argument settings for `srun` are available in this [HPCC
    manual](http://hpcc.ucr.edu/manuals_linux-cluster_jobs.html#partitions) or
    the HPCC section of this website
    [here](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/linux/linux/#job-submission-with-sbatch).
    Next, load the R version required for running the workflow with `module load`. Sometimes it may be necessary to
    first unload an active software version before loading another version, *e.g.* `module unload R`.

<!-- end list -->

``` sh
srun --x11 --partition=gen242 --mem=20gb --cpus-per-task 8 --ntasks 1 --time 20:00:00 --pty bash -l
# module unload R; module load load R/4.1.2
```

2.  Load a workflow template with the `genWorkenvir` function. This can be done
    from the command-line or from within R. However, only one of the two options needs to be used.

From command-line

``` sh
$ Rscript -e "systemPipeRdata::genWorkenvir(workflow='rnaseq')"
$ cd rnaseq
```

From R

``` r
library(systemPipeRdata)
genWorkenvir(workflow = "rnaseq")
setwd("rnaseq")
```

3.  Optional: if the user wishes to use another `Rmd` file than the template instance
    provided by the `genWorkenvir` function, then it can be copied or downloaded
    into the root directory of the workflow environment (*e.g.* with `cp` or `wget`).

4.  Now one can open from the root directory of the workflow the corresponding R Markdown script (*e.g.* systemPipeChIPseq.Rmd) using an R IDE, such as *nvim-r*, *ESS* or RStudio.
    Subsequently, the workflow can be run as outlined below. For learning purposes it is recommended to run workflows for the first time interactively. Once all workflow steps are understood and possibly
    modified to custom needs, one can run the workflow from start to finish with a
    single command using `runWF()`.

### Load packages

The `systemPipeR` package needs to be loaded to perform the analysis
steps shown in this report (H Backman and Girke 2016). The package allows users
to run the entire analysis workflow interactively or with a single command
while also generating the corresponding analysis report. For details
see `systemPipeR's` main [vignette](http://www.bioconductor.org/packages/devel/bioc/vignettes/systemPipeR/inst/doc/systemPipeR.html).

``` r
library(systemPipeR)
```

To apply workflows to custom data, the user needs to modify the *`targets`* file and if
necessary update the corresponding parameter (*`.cwl`* and *`.yml`*) files.
A collection of pre-generated *`.cwl`* and *`.yml`* files are provided in the *`param/cwl`* subdirectory
of each workflow template. They are also viewable in the GitHub repository of *`systemPipeRdata`* ([see
here](https://github.com/tgirke/systemPipeRdata/tree/master/inst/extdata/param/cwl)).
For more information of the structure of the *targets* file, please consult the documentation
[here](http://www.bioconductor.org/packages/devel/bioc/vignettes/systemPipeR/inst/doc/systemPipeR.html#25_structure_of_targets_file). More details about the new parameter files from systemPipeR can be found [here](http://www.bioconductor.org/packages/devel/bioc/vignettes/systemPipeR/inst/doc/systemPipeR.html#26_structure_of_the_new_param_files_and_construct_sysargs2_container).

### Import custom functions

Custom functions for the challenge projects can be imported with the source
command from a local R script (here [challengeProject\_Fct.R](https://raw.githubusercontent.com/tgirke/GEN242/main/content/en/tutorials/spchipseq/challengeProject_Fct.R)). Skip this step if such a script is not available. Alternatively, these
functions can be loaded from a custom R package.

``` r
source("challengeProject_Fct.R")
```

### Experiment definition provided by `targets` file

The `targets` file defines all FASTQ files and sample comparisons of the analysis workflow.
If needed the tab separated (TSV) version of this file can be downloaded from [here](https://github.com/tgirke/GEN242/tree/main/content/en/assignments/Projects/targets_files)
and the corresponding Google Sheet is [here](https://docs.google.com/spreadsheets/d/1DTgTGlZZscSPjlHOGdJC8QK4vvimN1BORjXKXzd_cfA/edit#gid=472150521).

``` r
targetspath <- "targetsPE.txt"
targets <- read.delim(targetspath, comment.char = "#")
DT::datatable(targets, options = list(scrollX = TRUE, autoWidth = TRUE))
```

<div id="htmlwidget-1" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-1">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"],["./data/SRR446027_1.fastq.gz","./data/SRR446028_1.fastq.gz","./data/SRR446029_1.fastq.gz","./data/SRR446030_1.fastq.gz","./data/SRR446031_1.fastq.gz","./data/SRR446032_1.fastq.gz","./data/SRR446033_1.fastq.gz","./data/SRR446034_1.fastq.gz","./data/SRR446035_1.fastq.gz","./data/SRR446036_1.fastq.gz","./data/SRR446037_1.fastq.gz","./data/SRR446038_1.fastq.gz","./data/SRR446039_1.fastq.gz","./data/SRR446040_1.fastq.gz","./data/SRR446041_1.fastq.gz","./data/SRR446042_1.fastq.gz","./data/SRR446043_1.fastq.gz","./data/SRR446044_1.fastq.gz"],["./data/SRR446027_2.fastq.gz","./data/SRR446028_2.fastq.gz","./data/SRR446029_2.fastq.gz","./data/SRR446030_2.fastq.gz","./data/SRR446031_2.fastq.gz","./data/SRR446032_2.fastq.gz","./data/SRR446033_2.fastq.gz","./data/SRR446034_2.fastq.gz","./data/SRR446035_2.fastq.gz","./data/SRR446036_2.fastq.gz","./data/SRR446037_2.fastq.gz","./data/SRR446038_2.fastq.gz","./data/SRR446039_2.fastq.gz","./data/SRR446040_2.fastq.gz","./data/SRR446041_2.fastq.gz","./data/SRR446042_2.fastq.gz","./data/SRR446043_2.fastq.gz","./data/SRR446044_2.fastq.gz"],["M1A","M1B","A1A","A1B","V1A","V1B","M6A","M6B","A6A","A6B","V6A","V6B","M12A","M12B","A12A","A12B","V12A","V12B"],["M1","M1","A1","A1","V1","V1","M6","M6","A6","A6","V6","V6","M12","M12","A12","A12","V12","V12"],["Mock.1h.A","Mock.1h.B","Avr.1h.A","Avr.1h.B","Vir.1h.A","Vir.1h.B","Mock.6h.A","Mock.6h.B","Avr.6h.A","Avr.6h.B","Vir.6h.A","Vir.6h.B","Mock.12h.A","Mock.12h.B","Avr.12h.A","Avr.12h.B","Vir.12h.A","Vir.12h.B"],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],["23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>FileName1<\/th>\n      <th>FileName2<\/th>\n      <th>SampleName<\/th>\n      <th>Factor<\/th>\n      <th>SampleLong<\/th>\n      <th>Experiment<\/th>\n      <th>Date<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"scrollX":true,"autoWidth":true,"columnDefs":[{"className":"dt-right","targets":6},{"orderable":false,"targets":0}],"order":[],"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

## Workflow steps

This tutorial will demonstrate how to build the workflow in an interactive mode,
appending each step. The workflow is constructed by connecting each step via
`appendStep` method. Each `SYSargsList` instance contains instructions for
processing a set of input files with a specific command-line or R software
and the paths to the corresponding outfiles generated by a particular command-line
software/step.

To create a workflow within *`systemPipeR`*, we can start by defining an empty
container and checking the directory structure:

``` r
library(systemPipeR)
sal <- SPRproject()
```

    ## Creating directory:  /home/tgirke/tmp/GEN242/content/en/tutorials/sprnaseq/data 
    ## Creating directory '/home/tgirke/tmp/GEN242/content/en/tutorials/sprnaseq/.SPRproject'
    ## Creating file '/home/tgirke/tmp/GEN242/content/en/tutorials/sprnaseq/.SPRproject/SYSargsList.yml'

``` r
sal
```

    ## Instance of 'SYSargsList': 
    ##  No workflow steps added

In case one wishes to run the entire workflow immediately, one can run it as follows:

``` r
sal <- importWF(sal, file_path = "systemPipeRNAseq.Rmd")  # Populates SYSargsList object with run instructions for all steps
sal <- runWF(sal)  # Runs workflow. This may take some time.
sal <- renderReport(sal)  # Renders report
```

### Required packages and resources

The `systemPipeR` package needs to be loaded (H Backman and Girke 2016).

``` r
appendStep(sal) <- LineWise(code = {
    library(systemPipeR)
}, step_name = "load_SPR")
```

### Read preprocessing

#### Read quality filtering and trimming

The function `preprocessReads` allows to apply predefined or custom
read preprocessing functions to all FASTQ files referenced in a
`SYSargsList` container, such as quality filtering or adapter trimming
routines. The paths to the resulting output FASTQ files are stored in the
`outfiles` slot of the `SYSargsList` object. The following example performs adapter trimming with
the `trimLRPatterns` function from the `Biostrings` package.
After the trimming step a new targets file is generated (here
`targets_trim.txt`) containing the paths to the trimmed FASTQ files.
The new targets file can be used for the next workflow step with an updated
`SYSargs2` instance, *e.g.* running the NGS alignments using the
trimmed FASTQ files.

Here, we are appending this step to the `SYSargsList` object created previously.
All the parameters are defined on the `preprocessReads/preprocessReads-pe.yml` and
`preprocessReads/preprocessReads-pe.cwl` files.

``` r
appendStep(sal) <- SYSargsList(step_name = "preprocessing", targets = "targetsPE.txt",
    dir = TRUE, wf_file = "preprocessReads/preprocessReads-pe.cwl",
    input_file = "preprocessReads/preprocessReads-pe.yml", dir_path = "param/cwl",
    inputvars = c(FileName1 = "_FASTQ_PATH1_", FileName2 = "_FASTQ_PATH2_",
        SampleName = "_SampleName_"), dependency = c("load_SPR"))
```

After, we can check the `trimLRPatterns` function in input parameter:

``` r
yamlinput(sal, "preprocessing")$Fct
```

    ## [1] "'trimLRPatterns(Rpattern=\"GCCCGGGTAA\", subject=fq)'"

After the preprocessing step, the `outfiles` files can be used to generate the new
targets files containing the paths to the trimmed FASTQ files. The new targets
information can be used for the next workflow step instance, *e.g.* running the
NGS alignments with the trimmed FASTQ files. The `appendStep` function is
automatically handling this connectivity between steps. Please check the `Alignments`
step for more details.

#### FASTQ quality report

The following `seeFastq` and `seeFastqPlot` functions generate and plot a series of useful
quality statistics for a set of FASTQ files including per cycle quality box
plots, base proportions, base-level quality trends, relative k-mer
diversity, length and occurrence distribution of reads, number of reads
above quality cutoffs and mean quality distribution. The results are
written to a PDF file named `fastqReport.pdf`.

``` r
appendStep(sal) <- LineWise(code = {
    fastq <- getColumn(sal, step = "preprocessing", "targetsWF",
        column = 1)
    fqlist <- seeFastq(fastq = fastq, batchsize = 10000, klength = 8)
    pdf("./results/fastqReport.pdf", height = 18, width = 4 *
        length(fqlist))
    seeFastqPlot(fqlist)
    dev.off()
}, step_name = "fastq_report", dependency = "preprocessing")
```

![](../results/fastqReport.png)

<div data-align="center">

Figure 1: FASTQ quality report for 18 samples

</div>

</br>

### Alignments

#### Read mapping with `HISAT2`

The following steps will demonstrate how to use the short read aligner `Hisat2`
(Kim, Langmead, and Salzberg 2015) in both interactive job submissions and batch submissions to
queuing systems of clusters using the *`systemPipeR's`* new CWL command-line interface.

The following steps will demonstrate how to use the short read aligner `Hisat2`
(Kim, Langmead, and Salzberg 2015). First, the `Hisat2` index needs to be created.

``` r
appendStep(sal) <- SYSargsList(step_name = "hisat2_index", dir = FALSE,
    targets = NULL, wf_file = "hisat2/hisat2-index.cwl", input_file = "hisat2/hisat2-index.yml",
    dir_path = "param/cwl", dependency = "load_SPR")
```

The parameter settings of the aligner are defined in the `workflow_hisat2-pe.cwl`
and `workflow_hisat2-pe.yml` files. The following shows how to construct the
corresponding *SYSargsList* object. Please note that the targets used in this
step are the outfiles from `preprocessing` step.

``` r
appendStep(sal) <- SYSargsList(step_name = "hisat2_mapping",
    dir = TRUE, targets = "preprocessing", wf_file = "workflow-hisat2/workflow_hisat2-pe.cwl",
    input_file = "workflow-hisat2/workflow_hisat2-pe.yml", dir_path = "param/cwl",
    inputvars = c(preprocessReads_1 = "_FASTQ_PATH1_", preprocessReads_2 = "_FASTQ_PATH2_",
        SampleName = "_SampleName_"), rm_targets_col = c("FileName1",
        "FileName2"), dependency = c("preprocessing", "hisat2_index"))
```

To double-check the command line for each sample, please use the following:

``` r
cmdlist(sal, step = "hisat2_mapping", targets = 1)
```

    ## $hisat2_mapping
    ## $hisat2_mapping$M1A
    ## $hisat2_mapping$M1A$hisat2
    ## [1] "hisat2 -S ./results/M1A.sam  -x ./data/tair10.fasta  -k 1  --min-intronlen 30  --max-intronlen 3000  -1 ./results/M1A_1.fastq_trim.gz -2 ./results/M1A_2.fastq_trim.gz --threads 4"
    ## 
    ## $hisat2_mapping$M1A$`samtools-view`
    ## [1] "samtools view -bS -o ./results/M1A.bam  ./results/M1A.sam "
    ## 
    ## $hisat2_mapping$M1A$`samtools-sort`
    ## [1] "samtools sort -o ./results/M1A.sorted.bam  ./results/M1A.bam  -@ 4"
    ## 
    ## $hisat2_mapping$M1A$`samtools-index`
    ## [1] "samtools index -b results/M1A.sorted.bam  results/M1A.sorted.bam.bai  ./results/M1A.sorted.bam "

#### Read and alignment stats

The following provides an overview of the number of reads in each sample
and how many of them aligned to the reference.

``` r
appendStep(sal) <- LineWise(code = {
    fqpaths <- getColumn(sal, step = "preprocessing", "targetsWF",
        column = "FileName1")
    bampaths <- getColumn(sal, step = "hisat2_mapping", "outfiles",
        column = "samtools_sort_bam")
    read_statsDF <- alignStats(args = bampaths, fqpaths = fqpaths,
        pairEnd = TRUE)
    write.table(read_statsDF, "results/alignStats.xls", row.names = FALSE,
        quote = FALSE, sep = "\t")
}, step_name = "align_stats", dependency = "hisat2_mapping")
```

The following shows the alignment statistics for a sample file provided by the `systemPipeR` package.

``` r
read.table("results/alignStats.xls", header = TRUE)[1:4, ]
```

    ##   FileName Nreads2x Nalign Perc_Aligned Nalign_Primary
    ## 1      M1A   115994 109977     94.81266         109977
    ## 2      M1B   134480 112464     83.62879         112464
    ## 3      A1A   127976 122427     95.66403         122427
    ## 4      A1B   122486 101369     82.75966         101369
    ##   Perc_Aligned_Primary
    ## 1             94.81266
    ## 2             83.62879
    ## 3             95.66403
    ## 4             82.75966

#### Create symbolic links for viewing BAM files in IGV

The `symLink2bam` function creates symbolic links to view the BAM alignment files in a
genome browser such as IGV without moving these large files to a local
system. The corresponding URLs are written to a file with a path
specified under `urlfile`, here `IGVurl.txt`.
Please replace the directory and the user name.

``` r
appendStep(sal) <- LineWise(code = {
    bampaths <- getColumn(sal, step = "hisat2_mapping", "outfiles",
        column = "samtools_sort_bam")
    symLink2bam(sysargs = bampaths, htmldir = c("~/.html/", "somedir/"),
        urlbase = "http://cluster.hpcc.ucr.edu/~<username>/",
        urlfile = "./results/IGVurl.txt")
}, step_name = "bam_urls", dependency = "hisat2_mapping", run_step = "optional")
```

### Read quantification

Reads overlapping with annotation ranges of interest are counted for
each sample using the `summarizeOverlaps` function (Lawrence et al. 2013). The read counting is
preformed for exonic gene regions in a non-strand-specific manner while
ignoring overlaps among different genes. Subsequently, the expression
count values are normalized by *reads per kp per million mapped reads*
(RPKM). The raw read count table (`countDFeByg.xls`) and the corresponding
RPKM table (`rpkmDFeByg.xls`) are written to separate files in the directory of
this project. Parallelization is achieved with the `BiocParallel` package, here
using 4 CPU cores.

#### Create a database for gene annotation

``` r
appendStep(sal) <- LineWise(code = {
    library(GenomicFeatures)
    txdb <- suppressWarnings(makeTxDbFromGFF(file = "data/tair10.gff",
        format = "gff", dataSource = "TAIR", organism = "Arabidopsis thaliana"))
    saveDb(txdb, file = "./data/tair10.sqlite")
}, step_name = "create_db", dependency = "hisat2_mapping")
```

#### Read counting with `summarizeOverlaps` in parallel mode using multiple cores

``` r
appendStep(sal) <- LineWise(code = {
    library(GenomicFeatures)
    library(BiocParallel)
    txdb <- loadDb("./data/tair10.sqlite")
    outpaths <- getColumn(sal, step = "hisat2_mapping", "outfiles",
        column = "samtools_sort_bam")
    eByg <- exonsBy(txdb, by = c("gene"))
    bfl <- BamFileList(outpaths, yieldSize = 50000, index = character())
    multicoreParam <- MulticoreParam(workers = 4)
    register(multicoreParam)
    registered()
    counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg,
        x, mode = "Union", ignore.strand = TRUE, inter.feature = FALSE,
        singleEnd = FALSE, BPPARAM = multicoreParam))
    countDFeByg <- sapply(seq(along = counteByg), function(x) assays(counteByg[[x]])$counts)
    rownames(countDFeByg) <- names(rowRanges(counteByg[[1]]))
    colnames(countDFeByg) <- names(bfl)
    rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts = x,
        ranges = eByg))
    write.table(countDFeByg, "results/countDFeByg.xls", col.names = NA,
        quote = FALSE, sep = "\t")
    write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names = NA,
        quote = FALSE, sep = "\t")
    ## Creating a SummarizedExperiment object
    colData <- data.frame(row.names = SampleName(sal, "hisat2_mapping"),
        condition = getColumn(sal, "hisat2_mapping", position = "targetsWF",
            column = "Factor"))
    colData$condition <- factor(colData$condition)
    countDF_se <- SummarizedExperiment::SummarizedExperiment(assays = countDFeByg,
        colData = colData)
    ## Add results as SummarizedExperiment to the workflow
    ## object
    SE(sal, "read_counting") <- countDF_se
}, step_name = "read_counting", dependency = "create_db")
```

When providing a `BamFileList` as in the example above, `summarizeOverlaps` methods
use by default `bplapply` and use the register interface from BiocParallel package.
If the number of workers is not set, `MulticoreParam` will use the number of cores
returned by `parallel::detectCores()`. For more information,
please check `help("summarizeOverlaps")` documentation.

Shows count table generated in previous step (`countDFeByg.xls`).
To avoid slowdowns of the load time of this page, ony 200 rows of the source
table are imported into the below `datatable` view .

``` r
countDF <- read.delim("results/countDFeByg.xls", row.names = 1,
    check.names = FALSE)[1:200, ]
DT::datatable(countDF, options = list(scrollX = TRUE, autoWidth = TRUE))
```

<div id="htmlwidget-2" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-2">{"x":{"filter":"none","vertical":false,"data":[["AT1G01010","AT1G01020","AT1G01030","AT1G01040","AT1G01046","AT1G01050","AT1G01060","AT1G01070","AT1G01073","AT1G01080","AT1G01090","AT1G01100","AT1G01110","AT1G01115","AT2G01008","AT2G01010","AT2G01020","AT2G01021","AT2G01023","AT3G01010","AT3G01015","AT3G01020","AT3G01030","AT3G01040","AT3G01050","AT3G01060","AT3G01070","AT3G01080","AT3G01085","AT3G01090","AT3G01100","AT3G01120","AT3G01130","AT3G01140","AT3G01142","AT3G01150","AT3G01160","AT3G01170","AT4G00005","AT4G00020","AT4G00026","AT4G00030","AT4G00040","AT4G00050","AT4G00060","AT4G00070","AT4G00080","AT4G00085","AT4G00090","AT4G00100","AT4G00110","AT4G00114","AT4G00120","AT4G00124","AT4G00130","AT4G00140","AT5G01010","AT5G01015","AT5G01020","AT5G01030","AT5G01040","AT5G01050","AT5G01060","AT5G01070","AT5G01075","AT5G01080","AT5G01090","AT5G01100","AT5G01110","AT5G01120","AT5G01130","AT5G01140","AT5G01150","AT5G01160","ATCG00010","ATCG00020","ATCG00030","ATCG00040","ATCG00050","ATCG00060","ATCG00070","ATCG00080","ATCG00090","ATCG00100","ATCG00110","ATCG00120","ATCG00130","ATCG00140","ATCG00150","ATCG00160","ATCG00170","ATCG00180","ATCG00190","ATCG00200","ATCG00210","ATCG00220","ATCG00230","ATCG00240","ATCG00250","ATCG00260","ATCG00270","ATCG00280","ATCG00290","ATCG00300","ATCG00310","ATCG00320","ATCG00330","ATCG00340","ATCG00350","ATCG00360","ATCG00370","ATCG00380","ATCG00390","ATCG00400","ATCG00410","ATCG00420","ATCG00430","ATCG00440","ATCG00450","ATCG00460","ATCG00470","ATCG00480","ATCG00490","ATCG00500","ATCG00510","ATMG00010","ATMG00020","ATMG00030","ATMG00040","ATMG00050","ATMG00060","ATMG00070","ATMG00080","ATMG00090","ATMG00100","ATMG00110","ATMG00120","ATMG00130","ATMG00140","ATMG00150","ATMG00160","ATMG00170","ATMG00180","ATMG00190","ATMG00200","NA","NA.1","NA.2","NA.3","NA.4","NA.5","NA.6","NA.7","NA.8","NA.9","NA.10","NA.11","NA.12","NA.13","NA.14","NA.15","NA.16","NA.17","NA.18","NA.19","NA.20","NA.21","NA.22","NA.23","NA.24","NA.25","NA.26","NA.27","NA.28","NA.29","NA.30","NA.31","NA.32","NA.33","NA.34","NA.35","NA.36","NA.37","NA.38","NA.39","NA.40","NA.41","NA.42","NA.43","NA.44","NA.45","NA.46","NA.47","NA.48","NA.49","NA.50","NA.51","NA.52","NA.53","NA.54"],[286,104,120,911,23,189,98,0,0,377,1167,871,11,0,20,28,49,1,0,0,0,0,0,381,125,890,7,158,2,337,106,2490,247,10,10,170,108,141,0,13,76,237,249,70,280,0,0,0,93,656,64,0,0,0,0,1,608,57,308,111,181,10,1,1,21,0,73,523,13,2,0,0,0,101,24,7014,6,43,3,1,31,92,48,35,2,287,157,199,61,15,48,60,123,2,45,45,2,54,54,0,722,1231,39,38,1,0,173,746,800,20,4,20,5,27,5,38,77,38,30,1,263,718,3399,21,54,3,2112,276,17,5,10,25,28,56,0,21,4,7,5,1,31,0,2,0,3,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null],[260,136,109,727,12,178,262,1,0,390,1195,478,10,0,10,52,15,1,1,0,0,0,0,507,143,1025,12,40,0,340,157,2764,143,21,21,152,257,141,0,30,127,164,197,161,409,0,0,0,78,440,119,0,0,0,0,0,513,28,241,205,156,3,0,2,21,0,109,1298,30,2,0,0,0,155,11,1866,10,40,4,1,62,93,53,15,4,308,225,64,74,18,75,84,186,5,25,155,2,52,11,0,463,944,31,67,0,3,128,576,642,9,4,28,28,25,5,157,404,110,19,9,562,1750,1286,12,6,1,737,60,21,3,11,34,34,95,0,47,2,7,2,3,41,0,6,0,1,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null],[364,139,167,1030,17,247,86,5,0,363,1323,961,16,0,10,16,33,2,0,0,0,0,0,693,216,935,9,317,0,573,177,3926,265,5,5,324,220,285,0,14,149,421,416,68,369,0,0,0,205,851,159,0,0,0,0,3,1024,11,342,159,274,17,1,4,12,0,84,1466,21,2,0,0,0,235,16,1890,2,13,1,0,41,97,42,9,0,234,82,111,47,12,20,41,99,3,13,21,0,51,50,0,426,588,18,27,0,0,113,471,450,14,4,5,0,6,2,34,72,33,13,0,215,772,1710,11,14,1,1165,171,15,8,10,19,51,72,0,27,3,8,6,7,35,0,2,0,1,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null],[181,131,59,627,13,184,88,0,0,454,1157,385,9,0,4,33,2,5,0,0,0,0,0,528,170,580,6,16,0,349,134,2183,119,36,36,176,295,176,0,62,71,160,182,140,365,0,0,0,97,407,141,0,0,0,0,0,551,19,196,217,108,3,0,0,21,0,102,1168,29,3,0,0,0,152,11,611,4,25,4,0,86,105,62,8,0,190,236,69,95,18,97,81,240,0,28,109,0,50,5,0,266,687,14,68,1,2,144,618,606,9,3,19,23,26,0,169,363,90,5,6,628,1828,909,11,15,0,623,41,18,4,3,38,76,179,0,51,9,3,11,0,49,0,1,0,1,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null],[568,174,136,962,16,226,33,10,0,476,1151,649,8,0,18,51,57,1,0,0,0,0,0,577,238,393,5,74,2,603,187,2674,370,9,9,250,231,230,0,25,92,201,289,35,472,2,0,0,116,669,96,0,0,0,0,0,980,8,313,171,264,7,0,0,9,1,77,1342,36,3,0,0,0,182,17,3808,2,13,4,1,53,116,50,11,1,275,148,94,41,12,57,64,175,1,31,60,0,24,24,0,383,780,10,36,0,0,97,558,575,21,3,11,2,22,2,40,70,37,36,1,218,812,1729,11,25,1,918,94,13,3,16,36,28,75,0,15,3,5,6,5,42,0,7,0,1,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null],[300,156,192,918,26,380,32,5,0,630,2312,1121,35,0,10,55,58,0,0,0,0,0,0,921,192,477,6,114,0,483,240,2866,316,6,6,478,157,395,0,46,75,302,342,28,543,0,0,0,185,1177,168,0,0,0,0,0,1188,14,561,136,526,14,1,1,20,0,101,2719,43,7,0,0,0,286,5,1520,6,16,0,0,63,108,59,4,0,380,121,62,56,18,72,92,270,1,8,59,2,9,6,0,401,576,13,34,0,0,82,653,654,12,8,11,5,4,1,92,156,37,7,0,599,1615,1443,18,5,4,828,57,12,2,17,31,51,71,1,29,4,6,8,7,55,0,2,0,1,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null],[255,148,74,862,19,524,8,11,0,437,1821,2090,5,0,33,43,130,2,0,0,0,0,0,188,188,41,1,89,0,443,186,3407,397,68,68,272,199,179,0,25,262,443,555,346,433,1,0,0,256,1746,118,0,0,0,0,3,1148,99,820,116,13,0,0,1,29,0,130,147,23,2,0,0,0,235,7,8203,3,9,0,0,20,54,21,5,1,185,68,176,22,6,18,18,54,3,15,20,5,83,89,0,427,891,10,20,0,0,97,597,674,10,4,11,1,16,5,26,36,17,10,0,111,387,3002,16,10,5,3061,515,34,21,14,24,30,56,0,25,5,12,10,2,35,0,6,0,3,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null],[135,131,26,618,14,619,4,13,0,747,1904,2156,11,0,23,55,73,4,1,0,0,0,0,221,235,133,3,88,0,537,243,4187,425,54,54,293,227,168,0,42,221,351,426,150,469,1,0,0,273,1687,103,0,0,0,0,0,1191,225,734,102,16,1,0,0,19,0,158,344,12,2,0,0,0,198,6,4264,5,9,0,0,26,32,21,5,2,209,80,105,33,8,30,51,80,0,21,10,0,17,17,0,306,675,4,7,0,0,67,483,536,11,3,6,1,9,4,27,42,12,14,1,184,642,2624,5,11,5,2099,201,47,28,12,23,24,73,1,34,3,10,6,5,48,0,6,0,3,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null],[514,114,23,880,23,382,6,29,0,266,1019,2443,11,0,31,18,64,1,1,0,0,0,1,218,181,41,0,197,0,440,238,3760,528,43,43,304,351,233,0,31,312,501,555,520,443,5,0,1,313,2312,98,0,0,0,0,2,1558,20,696,194,63,13,0,2,13,0,125,643,9,2,0,0,0,303,3,3679,1,5,0,0,14,27,17,0,0,139,26,94,16,4,15,12,43,0,5,12,1,41,44,0,287,464,5,7,0,2,61,343,338,3,6,3,0,4,3,16,34,6,5,1,80,292,1497,5,4,3,1972,398,36,22,22,31,32,56,0,33,2,4,8,6,55,0,3,0,1,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null],[318,104,73,639,21,414,3,8,0,350,1286,1429,13,0,12,90,128,10,0,0,1,0,0,165,146,40,0,77,0,352,148,2609,415,29,29,232,178,143,0,43,183,314,345,115,356,3,0,0,208,1331,107,0,0,0,0,0,1101,79,562,101,55,11,0,0,13,0,83,345,15,2,0,0,0,231,3,11266,6,28,2,0,15,38,18,18,0,368,199,227,48,16,35,47,76,1,52,18,0,11,11,1,742,1615,13,26,0,1,181,1073,1169,16,14,23,2,44,5,37,48,20,42,1,143,569,6908,26,34,2,5857,351,5,3,18,17,27,48,0,20,0,4,2,1,38,0,2,0,2,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null],[757,206,118,1632,24,622,2,28,0,352,1464,2733,28,0,33,18,39,1,0,0,1,0,1,289,259,34,1,129,0,557,331,4121,690,49,49,292,377,261,0,30,340,575,748,405,547,1,0,0,456,2631,112,1,0,0,0,2,1688,9,1103,215,26,2,2,1,20,0,181,602,17,1,0,0,0,321,7,2800,1,10,0,2,12,51,27,6,0,133,37,76,22,5,24,22,44,0,12,9,1,30,31,0,281,493,6,8,0,1,43,318,316,6,0,1,1,5,1,16,22,10,19,1,93,320,1138,5,13,8,1692,289,22,14,11,13,40,59,1,28,1,7,12,4,61,0,2,0,2,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null],[551,212,214,1552,36,962,10,14,0,765,2319,2760,40,0,14,35,55,1,0,0,0,0,1,316,403,52,2,46,0,647,385,3672,545,85,85,348,348,158,0,45,302,486,632,243,604,3,0,0,409,2447,129,0,0,0,0,0,1651,74,1324,190,25,1,0,1,44,1,206,219,18,1,0,0,0,393,9,562,1,5,0,0,14,29,23,3,0,122,34,40,15,5,25,27,69,2,6,9,0,12,11,0,126,305,2,6,0,0,37,240,232,1,2,4,0,2,2,14,24,6,8,1,127,572,648,3,5,3,660,65,27,14,11,29,46,95,1,24,4,17,14,3,61,0,4,0,1,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null],[198,67,45,651,28,666,220,11,0,384,1354,1439,14,0,16,110,39,7,3,0,0,0,0,208,135,156,0,54,0,468,270,3439,297,67,67,232,111,214,0,27,106,420,598,781,451,3,0,0,173,1395,96,0,0,0,0,1,896,8,1261,124,4,0,0,1,23,0,102,136,6,2,0,0,2,313,2,2551,0,5,0,2,15,9,12,3,0,106,18,27,7,1,13,30,50,0,8,7,0,0,3,0,207,267,0,5,0,0,35,299,329,4,4,4,2,3,1,26,56,10,0,2,141,474,830,3,2,4,1084,149,22,16,7,20,29,35,0,17,1,3,7,4,55,0,2,0,0,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null],[248,156,51,1095,23,1355,317,65,0,1037,2555,2452,42,0,50,47,124,1,1,0,0,0,0,209,347,304,2,62,0,624,336,3236,561,131,131,284,170,144,0,62,189,720,992,842,738,4,0,0,308,2118,147,0,0,0,0,1,1427,100,1849,178,7,0,0,1,36,0,228,91,39,0,0,0,0,425,8,714,0,4,1,0,18,49,35,5,2,143,50,52,18,4,19,28,73,3,8,15,0,19,17,0,131,331,1,3,0,0,41,238,237,8,2,6,0,2,1,20,34,14,18,1,174,687,662,4,5,3,702,124,26,16,8,26,36,77,0,35,5,8,9,3,35,0,7,0,0,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null],[527,130,31,1324,33,737,501,64,0,343,1100,2604,17,0,81,60,124,2,2,0,0,0,0,324,304,160,3,276,0,611,463,3093,767,73,73,279,203,292,0,39,226,615,817,933,639,13,0,0,386,2257,111,0,0,0,0,2,2050,23,1034,211,8,2,0,1,9,0,143,370,18,5,0,0,0,433,10,4167,2,10,0,0,10,51,25,6,0,119,48,91,23,9,12,24,49,0,12,14,0,35,33,0,319,619,0,9,0,0,65,434,506,14,3,8,0,6,3,16,34,11,12,1,111,369,1579,10,8,8,2356,567,35,19,22,18,35,49,0,26,3,12,16,5,59,0,4,0,1,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null],[417,120,48,702,13,532,198,39,0,299,941,812,26,0,16,137,139,8,1,0,0,0,0,260,253,55,1,115,0,376,219,1566,464,46,46,148,138,166,0,24,146,209,288,144,401,1,0,0,187,864,63,0,0,0,0,1,1258,17,660,123,8,0,0,0,12,0,65,174,16,2,0,0,0,238,3,13987,7,44,2,0,54,60,32,19,3,365,232,158,67,7,30,50,59,2,57,39,0,4,5,0,959,2153,21,23,0,0,208,1290,1195,19,5,32,6,82,5,43,74,18,51,0,195,746,5819,38,23,3,4738,331,12,10,19,21,31,65,0,23,0,3,7,2,59,0,2,0,2,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null],[650,80,177,671,23,635,164,23,0,267,1188,1546,16,0,31,131,108,4,2,0,0,0,0,318,148,149,1,136,0,405,248,2150,319,44,44,232,109,182,0,15,135,305,353,517,357,1,0,0,256,1623,93,1,0,0,0,4,920,2,901,136,11,1,0,0,15,0,91,478,9,1,0,0,1,327,3,4184,1,13,0,0,21,18,28,3,2,110,36,47,25,6,18,18,45,0,3,8,0,2,1,0,272,409,3,7,0,1,50,292,313,5,0,4,0,3,2,29,62,15,3,1,156,501,1138,7,0,8,1840,179,14,6,14,28,17,41,0,17,4,5,12,1,47,0,1,0,3,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null],[671,158,442,995,20,1004,159,24,0,373,1511,2145,33,0,27,113,125,3,2,0,0,0,0,321,231,111,2,92,0,464,211,1486,560,40,40,231,198,130,0,16,178,356,338,302,379,0,0,0,287,2106,65,0,0,0,0,1,994,27,1103,121,4,0,0,0,39,2,119,372,16,1,0,0,0,414,6,7310,3,9,3,0,30,51,35,17,1,205,101,173,34,8,22,32,56,2,37,33,0,24,25,0,434,1036,8,12,0,0,107,579,630,17,3,10,1,30,1,23,41,16,30,2,151,549,3176,9,21,6,2986,323,14,4,17,30,28,52,0,21,2,5,9,3,38,0,2,0,3,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>M1A<\/th>\n      <th>M1B<\/th>\n      <th>A1A<\/th>\n      <th>A1B<\/th>\n      <th>V1A<\/th>\n      <th>V1B<\/th>\n      <th>M6A<\/th>\n      <th>M6B<\/th>\n      <th>A6A<\/th>\n      <th>A6B<\/th>\n      <th>V6A<\/th>\n      <th>V6B<\/th>\n      <th>M12A<\/th>\n      <th>M12B<\/th>\n      <th>A12A<\/th>\n      <th>A12B<\/th>\n      <th>V12A<\/th>\n      <th>V12B<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"scrollX":true,"autoWidth":true,"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]},{"orderable":false,"targets":0}],"order":[],"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

A data slice of RPKM table (`rpkmDFeByg.xls`) is shown here.

``` r
read.delim("results/rpkmDFeByg.xls", row.names = 1, check.names = FALSE)[1:4,
    1:4]
```

    ##                M1A      M1B      A1A      A1B
    ## AT1G01010 5179.326 5955.322 7558.353 4856.097
    ## AT1G01020 1792.088 2964.078 2746.372 3344.252
    ## AT1G01030 1925.599 2212.258 3072.697 1402.614
    ## AT1G01040 4452.871 4494.494 5772.681 4540.367

Note, for most statistical differential expression or abundance analysis
methods, such as `edgeR` or `DESeq2`, the raw count values should be used as input. The
usage of RPKM values should be restricted to specialty applications
required by some users, *e.g.* manually comparing the expression levels
among different genes or features.

#### Sample-wise correlation analysis

The following computes the sample-wise Spearman correlation coefficients from
the `rlog` transformed expression values generated with the `DESeq2` package. After
transformation to a distance matrix, hierarchical clustering is performed with
the `hclust` function and the result is plotted as a dendrogram
(also see file `sample_tree.pdf`).

``` r
appendStep(sal) <- LineWise(code = {
    library(DESeq2, quietly = TRUE)
    library(ape, warn.conflicts = FALSE)
    ## Extracting SummarizedExperiment object
    se <- SE(sal, "read_counting")
    dds <- DESeqDataSet(se, design = ~condition)
    d <- cor(assay(rlog(dds)), method = "spearman")
    hc <- hclust(dist(1 - d))
    pdf("results/sample_tree.pdf")
    plot.phylo(as.phylo(hc), type = "p", edge.col = "blue", edge.width = 2,
        show.node.label = TRUE, no.margin = TRUE)
    dev.off()
}, step_name = "sample_tree", dependency = "read_counting")
```

![](../results/sample_tree.png)

<div data-align="center">

Figure 2: Correlation dendrogram of samples

</div>

</br>

### Analysis of DEGs

The analysis of differentially expressed genes (DEGs) is performed with
the glm method of the `edgeR` package (Robinson, McCarthy, and Smyth 2010). The sample
comparisons used by this analysis are defined in the header lines of the
`targets.txt` file starting with `<CMP>`.

#### Run `edgeR`

``` r
appendStep(sal) <- LineWise(code = {
    library(edgeR)
    countDF <- read.delim("results/countDFeByg.xls", row.names = 1,
        check.names = FALSE)
    cmp <- readComp(stepsWF(sal)[["hisat2_mapping"]], format = "matrix",
        delim = "-")
    edgeDF <- run_edgeR(countDF = countDF, targets = targetsWF(sal)[["hisat2_mapping"]],
        cmp = cmp[[1]], independent = FALSE, mdsplot = "")
}, step_name = "run_edger", dependency = "read_counting")
```

#### Add gene descriptions

``` r
appendStep(sal) <- LineWise(code = {
    library("biomaRt")
    m <- useMart("plants_mart", dataset = "athaliana_eg_gene",
        host = "https://plants.ensembl.org")
    desc <- getBM(attributes = c("tair_locus", "description"),
        mart = m)
    desc <- desc[!duplicated(desc[, 1]), ]
    descv <- as.character(desc[, 2])
    names(descv) <- as.character(desc[, 1])
    edgeDF <- data.frame(edgeDF, Desc = descv[rownames(edgeDF)],
        check.names = FALSE)
    write.table(edgeDF, "./results/edgeRglm_allcomp.xls", quote = FALSE,
        sep = "\t", col.names = NA)
}, step_name = "custom_annot", dependency = "run_edger")
```

### Plot DEG results

Filter and plot DEG results for up and down regulated genes. The
definition of *up* and *down* is given in the corresponding help
file. To open it, type `?filterDEGs` in the R console.

``` r
appendStep(sal) <- LineWise(code = {
    edgeDF <- read.delim("results/edgeRglm_allcomp.xls", row.names = 1,
        check.names = FALSE)
    pdf("results/DEGcounts.pdf")
    DEG_list <- filterDEGs(degDF = edgeDF, filter = c(Fold = 2,
        FDR = 20))
    dev.off()
    write.table(DEG_list$Summary, "./results/DEGcounts.xls",
        quote = FALSE, sep = "\t", row.names = FALSE)
}, step_name = "filter_degs", dependency = "custom_annot")
```

![](../results/DEGcounts.png)

<div data-align="center">

Figure 3: Up and down regulated DEGs with FDR of 1%

</div>

</br>

### Venn diagrams of DEG sets

The `overLapper` function can compute Venn intersects for large numbers of sample
sets (up to 20 or more) and plots 2-5 way Venn diagrams. A useful
feature is the possibility to combine the counts from several Venn
comparisons with the same number of sample sets in a single Venn diagram
(here for 4 up and down DEG sets).

``` r
appendStep(sal) <- LineWise(code = {
    vennsetup <- overLapper(DEG_list$Up[6:9], type = "vennsets")
    vennsetdown <- overLapper(DEG_list$Down[6:9], type = "vennsets")
    pdf("results/vennplot.pdf")
    vennPlot(list(vennsetup, vennsetdown), mymain = "", mysub = "",
        colmode = 2, ccol = c("blue", "red"))
    dev.off()
}, step_name = "venn_diagram", dependency = "filter_degs")
```

![](../results/vennplot.png)

<div data-align="center">

Figure 4: Venn Diagram for 4 Up and Down DEG Sets

</div>

</br>

### GO term enrichment analysis

#### Obtain gene-to-GO mappings

The following shows how to obtain gene-to-GO mappings from `biomaRt` (here for *A.
thaliana*) and how to organize them for the downstream GO term
enrichment analysis. Alternatively, the gene-to-GO mappings can be
obtained for many organisms from Bioconductor’s `*.db` genome annotation
packages or GO annotation files provided by various genome databases.
For each annotation this relatively slow preprocessing step needs to be
performed only once. Subsequently, the preprocessed data can be loaded
with the `load` function as shown in the next subsection.

``` r
appendStep(sal) <- LineWise(code = {
    library("biomaRt")
    # listMarts() # To choose BioMart database
    # listMarts(host='plants.ensembl.org')
    m <- useMart("plants_mart", host = "https://plants.ensembl.org")
    m <- useMart("plants_mart", dataset = "athaliana_eg_gene",
        host = "https://plants.ensembl.org")
    go <- getBM(attributes = c("go_id", "tair_locus", "namespace_1003"),
        mart = m)
    go <- go[go[, 3] != "", ]
    go[, 3] <- as.character(go[, 3])
    go[go[, 3] == "molecular_function", 3] <- "F"
    go[go[, 3] == "biological_process", 3] <- "P"
    go[go[, 3] == "cellular_component", 3] <- "C"
    go[1:4, ]
    if (!dir.exists("./data/GO"))
        dir.create("./data/GO")
    write.table(go, "data/GO/GOannotationsBiomart_mod.txt", quote = FALSE,
        row.names = FALSE, col.names = FALSE, sep = "\t")
    catdb <- makeCATdb(myfile = "data/GO/GOannotationsBiomart_mod.txt",
        lib = NULL, org = "", colno = c(1, 2, 3), idconv = NULL)
    save(catdb, file = "data/GO/catdb.RData")
}, step_name = "get_go_annot", dependency = "filter_degs")
```

#### Batch GO term enrichment analysis

Apply the enrichment analysis to the DEG sets obtained the above differential
expression analysis. Note, in the following example the `FDR` filter is set
here to an unreasonably high value, simply because of the small size of the toy
data set used in this vignette. Batch enrichment analysis of many gene sets is
performed with the function. When `method=all`, it returns all GO terms passing
the p-value cutoff specified under the `cutoff` arguments. When `method=slim`,
it returns only the GO terms specified under the `myslimv` argument. The given
example shows how a GO slim vector for a specific organism can be obtained from
`BioMart`.

``` r
appendStep(sal) <- LineWise(code = {
    library("biomaRt")
    load("data/GO/catdb.RData")
    DEG_list <- filterDEGs(degDF = edgeDF, filter = c(Fold = 2,
        FDR = 50), plot = FALSE)
    up_down <- DEG_list$UporDown
    names(up_down) <- paste(names(up_down), "_up_down", sep = "")
    up <- DEG_list$Up
    names(up) <- paste(names(up), "_up", sep = "")
    down <- DEG_list$Down
    names(down) <- paste(names(down), "_down", sep = "")
    DEGlist <- c(up_down, up, down)
    DEGlist <- DEGlist[sapply(DEGlist, length) > 0]
    BatchResult <- GOCluster_Report(catdb = catdb, setlist = DEGlist,
        method = "all", id_type = "gene", CLSZ = 2, cutoff = 0.9,
        gocats = c("MF", "BP", "CC"), recordSpecGO = NULL)
    m <- useMart("plants_mart", dataset = "athaliana_eg_gene",
        host = "https://plants.ensembl.org")
    goslimvec <- as.character(getBM(attributes = c("goslim_goa_accession"),
        mart = m)[, 1])
    BatchResultslim <- GOCluster_Report(catdb = catdb, setlist = DEGlist,
        method = "slim", id_type = "gene", myslimv = goslimvec,
        CLSZ = 10, cutoff = 0.01, gocats = c("MF", "BP", "CC"),
        recordSpecGO = NULL)
    write.table(BatchResultslim, "results/GOBatchSlim.xls", row.names = FALSE,
        quote = FALSE, sep = "\t")
}, step_name = "go_enrich", dependency = "get_go_annot")
```

Shows GO term enrichment results from previous step. The last gene identifier column (10)
of this table has been excluded in this viewing instance to minimize the complexity of the
result.
To avoid slowdowns of the load time of this page, only 10 rows of the source table are shown
below.

``` r
BatchResult <- read.delim("results/GOBatchAll.xls")[1:10, ]
knitr::kable(BatchResult[, -10])
```

| CLID            | CLSZ | GOID       | NodeSize | SampleMatch |    Phyper |      Padj | Term                                                           | Ont |
| :-------------- | ---: | :--------- | -------: | ----------: | --------: | --------: | :------------------------------------------------------------- | :-- |
| M1-A1\_up\_down |   26 | GO:0050291 |        4 |           1 | 0.0039621 | 0.0396207 | sphingosine N-acyltransferase activity                         | MF  |
| M1-A1\_up\_down |   26 | GO:0004345 |        6 |           1 | 0.0059375 | 0.0593750 | glucose-6-phosphate dehydrogenase activity                     | MF  |
| M1-A1\_up\_down |   26 | GO:0050664 |       11 |           1 | 0.0108597 | 0.1085975 | oxidoreductase activity, acting on NAD(P)H, oxygen as acceptor | MF  |
| M1-A1\_up\_down |   26 | GO:0052593 |       11 |           1 | 0.0108597 | 0.1085975 | tryptamine:oxygen oxidoreductase (deaminating) activity        | MF  |
| M1-A1\_up\_down |   26 | GO:0052594 |       11 |           1 | 0.0108597 | 0.1085975 | aminoacetone:oxygen oxidoreductase(deaminating) activity       | MF  |
| M1-A1\_up\_down |   26 | GO:0052595 |       11 |           1 | 0.0108597 | 0.1085975 | aliphatic-amine oxidase activity                               | MF  |
| M1-A1\_up\_down |   26 | GO:0052596 |       11 |           1 | 0.0108597 | 0.1085975 | phenethylamine:oxygen oxidoreductase (deaminating) activity    | MF  |
| M1-A1\_up\_down |   26 | GO:0052793 |       12 |           1 | 0.0118414 | 0.1184141 | pectin acetylesterase activity                                 | MF  |
| M1-A1\_up\_down |   26 | GO:0008131 |       15 |           1 | 0.0147808 | 0.1478083 | primary amine oxidase activity                                 | MF  |
| M1-A1\_up\_down |   26 | GO:0016018 |       16 |           1 | 0.0157588 | 0.1575878 | cyclosporin A binding                                          | MF  |

#### Plot batch GO term results

The `data.frame` generated by `GOCluster` can be plotted with the `goBarplot` function. Because of the
variable size of the sample sets, it may not always be desirable to show
the results from different DEG sets in the same bar plot. Plotting
single sample sets is achieved by subsetting the input data frame as
shown in the first line of the following example.

``` r
appendStep(sal) <- LineWise(code = {
    gos <- BatchResultslim[grep("M6-V6_up_down", BatchResultslim$CLID),
        ]
    gos <- BatchResultslim
    png("results/GOslimbarplotMF.png")
    goBarplot(gos, gocat = "MF")
    dev.off()
    png("results/GOslimbarplotBP.png")
    goBarplot(gos, gocat = "BP")
    dev.off()
    png("results/GOslimbarplotCC.png")
    goBarplot(gos, gocat = "CC")
    dev.off()
}, step_name = "go_plot", dependency = "go_enrich")
```

![](../results/GOslimbarplotMF.png)

<div data-align="center">

Figure 5: GO Slim Barplot for MF Ontology

</div>

</br>

### Clustering and heat maps

The following example performs hierarchical clustering on the `rlog`
transformed expression matrix subsetted by the DEGs identified in the above
differential expression analysis. It uses a Pearson correlation-based distance
measure and complete linkage for cluster joining.

``` r
appendStep(sal) <- LineWise(code = {
    library(pheatmap)
    geneids <- unique(as.character(unlist(DEG_list[[1]])))
    y <- assay(rlog(dds))[geneids, ]
    pdf("results/heatmap1.pdf")
    pheatmap(y, scale = "row", clustering_distance_rows = "correlation",
        clustering_distance_cols = "correlation")
    dev.off()
}, step_name = "heatmap", dependency = "go_enrich")
```

![](../results/heatmap1.png)

<div data-align="center">

Figure 6: Heat Map with Hierarchical Clustering Dendrograms of DEGs

</div>

</br>

### Version Information

``` r
appendStep(sal) <- LineWise(code = {
    sessionInfo()
}, step_name = "sessionInfo", dependency = "heatmap")
```

## Running workflow

### Interactive job submissions in a single machine

For running the workflow, `runWF` function will execute all the steps store in
the `SYSargsList` workflow container. The execution will be on a single machine without
submitting to a queuing system of a computer cluster. Besides, `runWF` allows the
user to create a dedicated results folder for each workflow. This includes all the
output and log files for each step. When these options are used, the output
location will be updated by default and can be assigned to the same object.

``` r
sal <- runWF(sal)
```

### Parallelization on clusters

Alternatively, the computation can be greatly accelerated by processing many files
in parallel using several compute nodes of a cluster, where a scheduling/queuing
system is used for load balancing.

The `resources` list object provides the number of independent parallel cluster
processes defined under the `Njobs` element in the list. The following example
will run 18 processes in parallel using each 4 CPU cores.
If the resources available on a cluster allow running all 18 processes at the
same time, then the shown sample submission will utilize in a total of 72 CPU cores.

Note, `runWF` can be used with most queueing systems as it is based on utilities
from the `batchtools` package, which supports the use of template files (*`*.tmpl`*)
for defining the run parameters of different schedulers. To run the following
code, one needs to have both a `conffile` (see *`.batchtools.conf.R`* samples [here](https://mllg.github.io/batchtools/))
and a `template` file (see *`*.tmpl`* samples [here](https://github.com/mllg/batchtools/tree/master/inst/templates))
for the queueing available on a system. The following example uses the sample
`conffile` and `template` files for the Slurm scheduler provided by this package.

The resources can be appended when the step is generated, or it is possible to
add these resources later, as the following example using the `addResources`
function:

``` r
resources <- list(conffile=".batchtools.conf.R",
                  template="batchtools.slurm.tmpl", 
                  Njobs=18, 
                  walltime=120, ## minutes
                  ntasks=1,
                  ncpus=4, 
                  memory=1024, ## Mb
                  partition = "short"
                  )
sal <- addResources(sal, step = c("hisat2_mapping"), resources = resources)
sal <- runWF(sal)
```

### Visualize workflow

*`systemPipeR`* workflows instances can be visualized with the `plotWF` function.

``` r
plotWF(sal, out_format = "html", out_path = "plotWF.html")
```

<div id="htmlwidget-3" style="width:672px;height:480px;" class="plotwf html-widget"></div>
<script type="application/json" data-for="htmlwidget-3">{"x":{"dot":"digraph {\n    node[fontsize=20];\n    subgraph {\n        load_SPR -> hisat2_index -> hisat2_mapping -> create_db -> read_counting -> run_edger -> custom_annot -> filter_degs -> get_go_annot -> go_enrich -> heatmap -> sessionInfo\n   }\n    load_SPR -> preprocessing\n    preprocessing -> hisat2_mapping\n    hisat2_mapping -> align_stats\n    hisat2_mapping -> bam_urls\n    read_counting -> sample_tree\n    filter_degs -> venn_diagram\n    go_enrich -> go_plot\n    preprocessing -> fastq_report\n    load_SPR[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">load_SPR<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step load_SPR: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-25 13:50:56; End time: 2022-04-25 13:50:56; Duration: 00:00:00\"]\n    preprocessing[fillcolor=\"#d3d6eb\" style=\"filled, rounded\" label=<<b><font color=\"black\">preprocessing<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">36<\/font><\/b>; <font color=\"black\">0s<\/font>> , shape=\"box\"  tooltip=\"step preprocessing: 0 samples passed; 0 samples have warnings; 0 samples have errors; 36 samples in total; Start time: 2022-04-25 13:50:56; End time: 2022-04-25 13:50:56; Duration: 00:00:00\"]\n    fastq_report[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">fastq_report<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step fastq_report: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-25 13:50:56; End time: 2022-04-25 13:50:56; Duration: 00:00:00\"]\n    hisat2_index[fillcolor=\"#d3d6eb\" style=\"filled, rounded\" label=<<b><font color=\"black\">hisat2_index<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">8<\/font><\/b>; <font color=\"black\">0s<\/font>> , shape=\"box\"  tooltip=\"step hisat2_index: 0 samples passed; 0 samples have warnings; 0 samples have errors; 8 samples in total; Start time: 2022-04-25 13:50:56; End time: 2022-04-25 13:50:56; Duration: 00:00:00\"]\n    hisat2_mapping[fillcolor=\"#d3d6eb\" style=\"filled, rounded\" label=<<b><font color=\"black\">hisat2_mapping<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">72<\/font><\/b>; <font color=\"black\">0s<\/font>> , shape=\"box\"  tooltip=\"step hisat2_mapping: 0 samples passed; 0 samples have warnings; 0 samples have errors; 72 samples in total; Start time: 2022-04-25 13:50:56; End time: 2022-04-25 13:50:56; Duration: 00:00:00\"]\n    align_stats[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">align_stats<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step align_stats: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-25 13:50:56; End time: 2022-04-25 13:50:56; Duration: 00:00:00\"]\n    bam_urls[style=\"solid, \"label=<<b><font color=\"black\">bam_urls<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step bam_urls: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-25 13:50:56; End time: 2022-04-25 13:50:56; Duration: 00:00:00\"]\n    create_db[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">create_db<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step create_db: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-25 13:50:56; End time: 2022-04-25 13:50:56; Duration: 00:00:00\"]\n    read_counting[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">read_counting<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step read_counting: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-25 13:50:56; End time: 2022-04-25 13:50:56; Duration: 00:00:00\"]\n    sample_tree[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">sample_tree<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step sample_tree: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-25 13:50:56; End time: 2022-04-25 13:50:56; Duration: 00:00:00\"]\n    run_edger[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">run_edger<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step run_edger: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-25 13:50:56; End time: 2022-04-25 13:50:56; Duration: 00:00:00\"]\n    custom_annot[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">custom_annot<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step custom_annot: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-25 13:50:56; End time: 2022-04-25 13:50:56; Duration: 00:00:00\"]\n    filter_degs[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">filter_degs<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step filter_degs: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-25 13:50:56; End time: 2022-04-25 13:50:56; Duration: 00:00:00\"]\n    venn_diagram[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">venn_diagram<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step venn_diagram: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-25 13:50:56; End time: 2022-04-25 13:50:56; Duration: 00:00:00\"]\n    get_go_annot[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">get_go_annot<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step get_go_annot: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-25 13:50:56; End time: 2022-04-25 13:50:56; Duration: 00:00:00\"]\n    go_enrich[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">go_enrich<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step go_enrich: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-25 13:50:56; End time: 2022-04-25 13:50:56; Duration: 00:00:00\"]\n    go_plot[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">go_plot<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step go_plot: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-25 13:50:56; End time: 2022-04-25 13:50:56; Duration: 00:00:00\"]\n    heatmap[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">heatmap<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step heatmap: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-25 13:50:56; End time: 2022-04-25 13:50:56; Duration: 00:00:00\"]\n    sessionInfo[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">sessionInfo<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step sessionInfo: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-25 13:50:56; End time: 2022-04-25 13:50:56; Duration: 00:00:00\"]\n        subgraph cluster_legend {\n        rankdir=TB;\n        color=\"#eeeeee\";\n        style=filled;\n        ranksep =1;\n        label=\"Legends\";\n        fontsize = 30;\n        node [style=filled, fontsize=10];\n        legend_img-> step_state[color=\"#eeeeee\"];\n\n        legend_img[shape=none, image=\"plotwf_legend-src.png\", label = \" \", height=1, width=3, style=\"\"];\n\n        step_state[style=\"filled\", shape=\"box\" color=white, label =<\n            <table>\n            <tr><td><b>Step Colors<\/b><\/td><\/tr>\n            <tr><td><font color=\"black\">Pending steps<\/font>; <font color=\"#5cb85c\">Successful steps<\/font>; <font color=\"#d9534f\">Failed steps<\/font><\/td><\/tr>\n            <tr><td><b>Targets Files / Code Chunk <\/b><\/td><\/tr><tr><td><font color=\"#5cb85c\">0 (pass) <\/font> | <font color=\"#f0ad4e\">0 (warning) <\/font> | <font color=\"#d9534f\">0 (error) <\/font> | <font color=\"blue\">0 (total)<\/font>; Duration<\/td><\/tr><\/table>\n            >];\n\n    }\n\n}\n","plotid":"sprwf-43175286","responsive":true,"width":null,"height":null,"plot_method":"renderSVGElement","rmd":true,"msg":"","plot_ctr":true,"pan_zoom":false,"legend_uri":"data:image/svg+xml;base64,PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iVVRGLTgiPz4KPCEtLSBEbyBub3QgZWRpdCB0aGlzIGZpbGUgd2l0aCBlZGl0b3JzIG90aGVyIHRoYW4gZGlhZ3JhbXMubmV0IC0tPgo8IURPQ1RZUEUgc3ZnIFBVQkxJQyAiLS8vVzNDLy9EVEQgU1ZHIDEuMS8vRU4iICJodHRwOi8vd3d3LnczLm9yZy9HcmFwaGljcy9TVkcvMS4xL0RURC9zdmcxMS5kdGQiPgo8c3ZnIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyIgeG1sbnM6eGxpbms9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGxpbmsiIHZlcnNpb249IjEuMSIgd2lkdGg9IjQ5NnB4IiBoZWlnaHQ9IjI3OHB4IiB2aWV3Qm94PSItMC41IC0wLjUgNDk2IDI3OCIgY29udGVudD0iJmx0O214ZmlsZSBob3N0PSZxdW90O2FwcC5kaWFncmFtcy5uZXQmcXVvdDsgbW9kaWZpZWQ9JnF1b3Q7MjAyMS0xMS0yNFQyMDozOTo0NC45MjNaJnF1b3Q7IGFnZW50PSZxdW90OzUuMCAoWDExOyBMaW51eCB4ODZfNjQpIEFwcGxlV2ViS2l0LzUzNy4zNiAoS0hUTUwsIGxpa2UgR2Vja28pIENocm9tZS85Ni4wLjQ2NjQuNDUgU2FmYXJpLzUzNy4zNiZxdW90OyB2ZXJzaW9uPSZxdW90OzE1LjguMyZxdW90OyBldGFnPSZxdW90O1RfZEV4bkw3U0xYMmJ1WDhQYVgzJnF1b3Q7IHR5cGU9JnF1b3Q7Z29vZ2xlJnF1b3Q7Jmd0OyZsdDtkaWFncmFtIGlkPSZxdW90O1Z3OVZWaUdzeC1zTF9OaktlN0pQJnF1b3Q7Jmd0OzdWcmJjdUk0RVAwYUhwT3lKZDk0REFSbVhtWjNLMnpWUGl0WUdGVmt5V1BMZ2N6WHIyVEwrQ0lEM3NHUVRSVWtsZGd0cVMyZDAycDF0NW5BZWJ6L2xxSmsrNE9IbUU2QUZlNG44SGtDZ0EwQW1LaGZLL3pRRXN2eFMwbVVrbERMYXNHSy9NSlZSeTNOU1lpelZrZkJPUlVrYVF2WG5ERzhGaTBaU2xPK2EzZmJjTnArYW9JaWJBaFdhMFJONlQ4a0ZOdFNPcldzV3Y0ZGsyaGJQZG1yV21KVWRkYUNiSXRDdm11STRHSUM1eW5ub3J5SzkzTk1GWG9WTHVXNDVaSFd3OFJTek1TUUFacUpkMFJ6dlRZOUwvRlJMVFpLZVo1TTRFeitZeUZXNHl4NWQ1aTR1dEZLY0Nyd3ZvOEE5Rm9wczh3SjJvZGxTNFBCUE1ZaS9aQmRLa1dCSHFKdEJVQjl2NnVCQjQ1YnlyWk4wS0duQ2Rka1J3ZmROUjd5UWtQU0R3ODhEMDhibHQyV0NMeEswRnExN3FUNVM5bFd4Rkwvc3kwdk40VFNPYWM4TGNiQzBNVkI2Q2cwUmNyZmNOWENPSlBEWnhGRldkWUhkL2FHeFhwN0R2c214cUFmWTQycEE0ZEI2b3lBcUhOZFJKZkxSZkM4T0lib1ZWRU1wZ05SOUM5SDBiMHVpb3ZsRWl6bnQwWHhMR3p3Y3RpOHNXSGpUT2pUS1RCUjNManFSOG9SSlJHVE1vbzM0aGlvU2xWamJQbVJjaTRmVG9ReXNXQkVoSDNyU2c3VE54RE9PSldOWFpqbDlFVWJ5d3FrdFZ3R2xpak0xQ0tKUEhPZmRFTk13bEFObjZVNEk3LzBrYUl3U1RoaG9waXpPNXU0ejBwWExuaFc4bUlmQmJ4QlZvT0RKcDN5ZG9saVFoWDZmNU5ZUmh6QStnUHY1TjhYSGlNMmxJN1RCeHp3L1VlM2ZjUUIwNVBBSHNLcXcvUVN2Z0tEcnhCbFczd243RGhoMEhVTndnS0RNT2RLaEUwTnduNGdKdVBOV00zK3E1Rm0yN2ZiWm03dzZIZG9jODFRc3M4eGprRmJsYmswZUp2ek9Na0Z2cE4yWXFzNXhsYnI0Y3k3Rm1kMlQ3amcwZUlJVCtRU202eDVQM09Wc2hXNFBKVDRQc2tPUWJLdjIrUlZwUDdMcDhoRUZWZTY1RFJLZFdXcllROVY5LzhVeElWNGczSXFicmU5dktETGxPOFpUTG4yWTA5a0IwWUlpRzB6a1gyNTc2d1RmQVYyaHk5bzlXU0IxOXBaWmw0dHZhRmNYUGhBQ2J1N3hGTXVFWFpkSXJSc2c3amdXc1NaNlh1cTZqc3NvaVpyZFNwbFgrSzhQaUhBQTBOOEdlajFaZTRJSUEvSTdqRUxuMVFoc3piSEJwYWRTbDFoejFXZEVoNHd3cUZSNUR5TFVIUDlQU1pXeVZKTWtTRHZiZlY5a09nbi9LWDJaUE13bVQ2Q0RnUGRxa25HODNTTjlVRFFLSEIyZEVGZ0czRmZWNWRBYVlTRm9hdmc2YkQ0WWRRTnFEQmNSSjNPenNyTy8zY2VaYWJVQWQ1N0RLeHAvWEYvajFQSE1sTG1NNXBIWk5pc2NNZ0VMRVNDNjAxK1A3ZU9CSWhkcHdvZDA2bGVMWlEzNnh4L0pvSndodWlkdEZNbm9UV0F0RDRQTWdwcFpxMWpsUHhMd1RzMCticm83UTRNUGZ4NnF1QjdpK3BpKy9VWkJHYVU3L2k5a2N3SURBS3o2akVLZzR6TGxsdlJ1Q3crbjB1akdaTGVtRWl6RlBLU00wWllORkhxU3hKZTB5R0UycFprdEtDeVErb0taNW4weUYvSUgrczVkdDRFMmRiTnpNSU9nbllZQkwzS0FzNzRaemlHVVpnMWw4b1VRdkora1MyYzlSS3ZhUDBXRmJ2NllWMHlvL1NKRkxHc2dtK21mTVRBNXgyTStZZzdPWWlMaGVreCtHZE8wcTlXZVA5OG0vWGJOZ3VCR1ZIMFZlRkhzVml6N3JRU09CblBoODJMcjJ6Y3JXRzROVXhoMjRNRlpqR3I3MVhhYjFpRHZLMi9WMVZtZnZYWDArRGlYdz09Jmx0Oy9kaWFncmFtJmd0OyZsdDsvbXhmaWxlJmd0OyI+PGRlZnMvPjxnPjxyZWN0IHg9IjQiIHk9IjkwIiB3aWR0aD0iNDkwIiBoZWlnaHQ9IjkyIiBmaWxsPSIjZDVlOGQ0IiBzdHJva2U9Im5vbmUiIHBvaW50ZXItZXZlbnRzPSJub25lIi8+PHJlY3QgeD0iNCIgeT0iMTgyIiB3aWR0aD0iNDkwIiBoZWlnaHQ9Ijk0IiBmaWxsPSIjZmZlOGRlIiBzdHJva2U9Im5vbmUiIHBvaW50ZXItZXZlbnRzPSJub25lIi8+PHJlY3QgeD0iNCIgeT0iNCIgd2lkdGg9IjQ5MCIgaGVpZ2h0PSI4NiIgZmlsbD0iI2VmZjJmYyIgc3Ryb2tlPSJub25lIiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxyZWN0IHg9IjQiIHk9IjQiIHdpZHRoPSIxNDAiIGhlaWdodD0iMjcyIiBmaWxsLW9wYWNpdHk9IjAuOCIgZmlsbD0iI2Y1ZjVmNSIgc3Ryb2tlPSJub25lIiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiAxMXB4OyBtYXJnaW4tbGVmdDogMTE1cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogOHB4OyBmb250LWZhbWlseTogJnF1b3Q7VGltZXMgTmV3IFJvbWFuJnF1b3Q7OyBjb2xvcjogcmdiKDAsIDAsIDApOyBsaW5lLWhlaWdodDogMS4yOyBwb2ludGVyLWV2ZW50czogbm9uZTsgd2hpdGUtc3BhY2U6IG5vd3JhcDsiPnNvbGlkPC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjExNSIgeT0iMTMiIGZpbGw9InJnYmEoMCwgMCwgMCwgMSkiIGZvbnQtZmFtaWx5PSJUaW1lcyBOZXcgUm9tYW4iIGZvbnQtc2l6ZT0iOHB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5zb2xpZDwvdGV4dD48L3N3aXRjaD48L2c+PGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTAuNSAtMC41KXNjYWxlKDIpIj48c3dpdGNoPjxmb3JlaWduT2JqZWN0IHBvaW50ZXItZXZlbnRzPSJub25lIiB3aWR0aD0iMTAwJSIgaGVpZ2h0PSIxMDAlIiByZXF1aXJlZEZlYXR1cmVzPSJodHRwOi8vd3d3LnczLm9yZy9UUi9TVkcxMS9mZWF0dXJlI0V4dGVuc2liaWxpdHkiIHN0eWxlPSJvdmVyZmxvdzogdmlzaWJsZTsgdGV4dC1hbGlnbjogbGVmdDsiPjxkaXYgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGh0bWwiIHN0eWxlPSJkaXNwbGF5OiBmbGV4OyBhbGlnbi1pdGVtczogdW5zYWZlIGNlbnRlcjsganVzdGlmeS1jb250ZW50OiB1bnNhZmUgY2VudGVyOyB3aWR0aDogMXB4OyBoZWlnaHQ6IDFweDsgcGFkZGluZy10b3A6IDEwcHg7IG1hcmdpbi1sZWZ0OiAxOThweDsiPjxkaXYgZGF0YS1kcmF3aW8tY29sb3JzPSJjb2xvcjogcmdiYSgwLCAwLCAwLCAxKTsgIiBzdHlsZT0iYm94LXNpemluZzogYm9yZGVyLWJveDsgZm9udC1zaXplOiAwcHg7IHRleHQtYWxpZ246IGNlbnRlcjsiPjxkaXYgc3R5bGU9ImRpc3BsYXk6IGlubGluZS1ibG9jazsgZm9udC1zaXplOiA4cHg7IGZvbnQtZmFtaWx5OiAmcXVvdDtUaW1lcyBOZXcgUm9tYW4mcXVvdDs7IGNvbG9yOiByZ2IoMCwgMCwgMCk7IGxpbmUtaGVpZ2h0OiAxLjI7IHBvaW50ZXItZXZlbnRzOiBub25lOyB3aGl0ZS1zcGFjZTogbm93cmFwOyI+ZGFzaGVkPC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjE5OCIgeT0iMTIiIGZpbGw9InJnYmEoMCwgMCwgMCwgMSkiIGZvbnQtZmFtaWx5PSJUaW1lcyBOZXcgUm9tYW4iIGZvbnQtc2l6ZT0iOHB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5kYXNoZWQ8L3RleHQ+PC9zd2l0Y2g+PC9nPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiAzMnB4OyBtYXJnaW4tbGVmdDogMTE2cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTFweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3dyYXA7Ij5NYW5hZ2VtZW50PC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjExNiIgeT0iMzUiIGZpbGw9InJnYmEoMCwgMCwgMCwgMSkiIGZvbnQtZmFtaWx5PSJUaW1lcyBOZXcgUm9tYW4iIGZvbnQtc2l6ZT0iMTFweCIgdGV4dC1hbmNob3I9Im1pZGRsZSI+TWFuYWdlbWVudDwvdGV4dD48L3N3aXRjaD48L2c+PGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTAuNSAtMC41KXNjYWxlKDIpIj48c3dpdGNoPjxmb3JlaWduT2JqZWN0IHBvaW50ZXItZXZlbnRzPSJub25lIiB3aWR0aD0iMTAwJSIgaGVpZ2h0PSIxMDAlIiByZXF1aXJlZEZlYXR1cmVzPSJodHRwOi8vd3d3LnczLm9yZy9UUi9TVkcxMS9mZWF0dXJlI0V4dGVuc2liaWxpdHkiIHN0eWxlPSJvdmVyZmxvdzogdmlzaWJsZTsgdGV4dC1hbGlnbjogbGVmdDsiPjxkaXYgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGh0bWwiIHN0eWxlPSJkaXNwbGF5OiBmbGV4OyBhbGlnbi1pdGVtczogdW5zYWZlIGNlbnRlcjsganVzdGlmeS1jb250ZW50OiB1bnNhZmUgY2VudGVyOyB3aWR0aDogMXB4OyBoZWlnaHQ6IDFweDsgcGFkZGluZy10b3A6IDMycHg7IG1hcmdpbi1sZWZ0OiAxOThweDsiPjxkaXYgZGF0YS1kcmF3aW8tY29sb3JzPSJjb2xvcjogcmdiYSgwLCAwLCAwLCAxKTsgIiBzdHlsZT0iYm94LXNpemluZzogYm9yZGVyLWJveDsgZm9udC1zaXplOiAwcHg7IHRleHQtYWxpZ246IGNlbnRlcjsiPjxkaXYgc3R5bGU9ImRpc3BsYXk6IGlubGluZS1ibG9jazsgZm9udC1zaXplOiAxMXB4OyBmb250LWZhbWlseTogJnF1b3Q7VGltZXMgTmV3IFJvbWFuJnF1b3Q7OyBjb2xvcjogcmdiKDAsIDAsIDApOyBsaW5lLWhlaWdodDogMS4yOyBwb2ludGVyLWV2ZW50czogbm9uZTsgd2hpdGUtc3BhY2U6IG5vd3JhcDsiPkNvbXB1dGU8L2Rpdj48L2Rpdj48L2Rpdj48L2ZvcmVpZ25PYmplY3Q+PHRleHQgeD0iMTk4IiB5PSIzNSIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMXB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5Db21wdXRlPC90ZXh0Pjwvc3dpdGNoPjwvZz48ZWxsaXBzZSBjeD0iMjMyLjUiIGN5PSIxMjMiIHJ4PSI1MS41IiByeT0iMjciIGZpbGw9InJnYmEoMjU1LCAyNTUsIDI1NSwgMSkiIHN0cm9rZT0icmdiYSgwLCAwLCAwLCAxKSIgc3Ryb2tlLXdpZHRoPSIyIiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDUwcHg7IGhlaWdodDogMXB4OyBwYWRkaW5nLXRvcDogNjJweDsgbWFyZ2luLWxlZnQ6IDkycHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTJweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3JtYWw7IG92ZXJmbG93LXdyYXA6IG5vcm1hbDsiPjxzcGFuIHN0eWxlPSJmb250LXNpemU6IDhweCI+ZWxsaXBzZTwvc3Bhbj48L2Rpdj48L2Rpdj48L2Rpdj48L2ZvcmVpZ25PYmplY3Q+PHRleHQgeD0iMTE2IiB5PSI2NSIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMnB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5lbGxpcHNlPC90ZXh0Pjwvc3dpdGNoPjwvZz48ZyB0cmFuc2Zvcm09InRyYW5zbGF0ZSgtMC41IC0wLjUpc2NhbGUoMikiPjxzd2l0Y2g+PGZvcmVpZ25PYmplY3QgcG9pbnRlci1ldmVudHM9Im5vbmUiIHdpZHRoPSIxMDAlIiBoZWlnaHQ9IjEwMCUiIHJlcXVpcmVkRmVhdHVyZXM9Imh0dHA6Ly93d3cudzMub3JnL1RSL1NWRzExL2ZlYXR1cmUjRXh0ZW5zaWJpbGl0eSIgc3R5bGU9Im92ZXJmbG93OiB2aXNpYmxlOyB0ZXh0LWFsaWduOiBsZWZ0OyI+PGRpdiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMTk5OS94aHRtbCIgc3R5bGU9ImRpc3BsYXk6IGZsZXg7IGFsaWduLWl0ZW1zOiB1bnNhZmUgY2VudGVyOyBqdXN0aWZ5LWNvbnRlbnQ6IHVuc2FmZSBjZW50ZXI7IHdpZHRoOiAxcHg7IGhlaWdodDogMXB4OyBwYWRkaW5nLXRvcDogODVweDsgbWFyZ2luLWxlZnQ6IDExNHB4OyI+PGRpdiBkYXRhLWRyYXdpby1jb2xvcnM9ImNvbG9yOiByZ2JhKDAsIDAsIDAsIDEpOyAiIHN0eWxlPSJib3gtc2l6aW5nOiBib3JkZXItYm94OyBmb250LXNpemU6IDBweDsgdGV4dC1hbGlnbjogY2VudGVyOyI+PGRpdiBzdHlsZT0iZGlzcGxheTogaW5saW5lLWJsb2NrOyBmb250LXNpemU6IDExcHg7IGZvbnQtZmFtaWx5OiAmcXVvdDtUaW1lcyBOZXcgUm9tYW4mcXVvdDs7IGNvbG9yOiByZ2IoMCwgMCwgMCk7IGxpbmUtaGVpZ2h0OiAxLjI7IHBvaW50ZXItZXZlbnRzOiBub25lOyB3aGl0ZS1zcGFjZTogbm93cmFwOyI+UjwvZGl2PjwvZGl2PjwvZGl2PjwvZm9yZWlnbk9iamVjdD48dGV4dCB4PSIxMTQiIHk9Ijg4IiBmaWxsPSJyZ2JhKDAsIDAsIDAsIDEpIiBmb250LWZhbWlseT0iVGltZXMgTmV3IFJvbWFuIiBmb250LXNpemU9IjExcHgiIHRleHQtYW5jaG9yPSJtaWRkbGUiPlI8L3RleHQ+PC9zd2l0Y2g+PC9nPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiA4M3B4OyBtYXJnaW4tbGVmdDogMTk4cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTFweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3dyYXA7Ij5Db21tYW5kLWxpbmU8L2Rpdj48L2Rpdj48L2Rpdj48L2ZvcmVpZ25PYmplY3Q+PHRleHQgeD0iMTk4IiB5PSI4NiIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMXB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5Db21tYW5kLWxpbmU8L3RleHQ+PC9zd2l0Y2g+PC9nPjxyZWN0IHg9IjM0OSIgeT0iOTYiIHdpZHRoPSIxMDUiIGhlaWdodD0iNTAiIHJ4PSI3LjUiIHJ5PSI3LjUiIGZpbGw9InJnYmEoMjU1LCAyNTUsIDI1NSwgMSkiIHN0cm9rZT0icmdiYSgwLCAwLCAwLCAxKSIgc3Ryb2tlLXdpZHRoPSIyIiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDUxcHg7IGhlaWdodDogMXB4OyBwYWRkaW5nLXRvcDogNjFweDsgbWFyZ2luLWxlZnQ6IDE3NnB4OyI+PGRpdiBkYXRhLWRyYXdpby1jb2xvcnM9ImNvbG9yOiByZ2JhKDAsIDAsIDAsIDEpOyAiIHN0eWxlPSJib3gtc2l6aW5nOiBib3JkZXItYm94OyBmb250LXNpemU6IDBweDsgdGV4dC1hbGlnbjogY2VudGVyOyI+PGRpdiBzdHlsZT0iZGlzcGxheTogaW5saW5lLWJsb2NrOyBmb250LXNpemU6IDhweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3JtYWw7IG92ZXJmbG93LXdyYXA6IG5vcm1hbDsiPnJlY3RhbmdsZTwvZGl2PjwvZGl2PjwvZGl2PjwvZm9yZWlnbk9iamVjdD48dGV4dCB4PSIyMDEiIHk9IjYzIiBmaWxsPSJyZ2JhKDAsIDAsIDAsIDEpIiBmb250LWZhbWlseT0iVGltZXMgTmV3IFJvbWFuIiBmb250LXNpemU9IjhweCIgdGV4dC1hbmNob3I9Im1pZGRsZSI+cmVjdGFuZ2xlPC90ZXh0Pjwvc3dpdGNoPjwvZz48cGF0aCBkPSJNIDE4Mi41IDM4IEwgMjg3LjUgMzgiIGZpbGw9Im5vbmUiIHN0cm9rZT0icmdiYSgwLCAwLCAwLCAxKSIgc3Ryb2tlLXdpZHRoPSI2IiBzdHJva2UtbWl0ZXJsaW1pdD0iMTAiIHBvaW50ZXItZXZlbnRzPSJub25lIi8+PHBhdGggZD0iTSAzNTQgMzcuNjIgTCA0NTkgMzcuNjIiIGZpbGw9Im5vbmUiIHN0cm9rZT0icmdiYSgwLCAwLCAwLCAxKSIgc3Ryb2tlLXdpZHRoPSI2IiBzdHJva2UtbWl0ZXJsaW1pdD0iMTAiIHN0cm9rZS1kYXNoYXJyYXk9IjE4IDE4IiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiAxMjhweDsgbWFyZ2luLWxlZnQ6IDExNXB4OyI+PGRpdiBkYXRhLWRyYXdpby1jb2xvcnM9ImNvbG9yOiByZ2JhKDAsIDAsIDAsIDEpOyAiIHN0eWxlPSJib3gtc2l6aW5nOiBib3JkZXItYm94OyBmb250LXNpemU6IDBweDsgdGV4dC1hbGlnbjogY2VudGVyOyI+PGRpdiBzdHlsZT0iZGlzcGxheTogaW5saW5lLWJsb2NrOyBmb250LXNpemU6IDExcHg7IGZvbnQtZmFtaWx5OiAmcXVvdDtUaW1lcyBOZXcgUm9tYW4mcXVvdDs7IGNvbG9yOiByZ2IoMCwgMCwgMCk7IGxpbmUtaGVpZ2h0OiAxLjI7IHBvaW50ZXItZXZlbnRzOiBub25lOyB3aGl0ZS1zcGFjZTogbm93cmFwOyI+TWFuZGF0b3J5PC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjExNSIgeT0iMTMxIiBmaWxsPSJyZ2JhKDAsIDAsIDAsIDEpIiBmb250LWZhbWlseT0iVGltZXMgTmV3IFJvbWFuIiBmb250LXNpemU9IjExcHgiIHRleHQtYW5jaG9yPSJtaWRkbGUiPk1hbmRhdG9yeTwvdGV4dD48L3N3aXRjaD48L2c+PGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTAuNSAtMC41KXNjYWxlKDIpIj48c3dpdGNoPjxmb3JlaWduT2JqZWN0IHBvaW50ZXItZXZlbnRzPSJub25lIiB3aWR0aD0iMTAwJSIgaGVpZ2h0PSIxMDAlIiByZXF1aXJlZEZlYXR1cmVzPSJodHRwOi8vd3d3LnczLm9yZy9UUi9TVkcxMS9mZWF0dXJlI0V4dGVuc2liaWxpdHkiIHN0eWxlPSJvdmVyZmxvdzogdmlzaWJsZTsgdGV4dC1hbGlnbjogbGVmdDsiPjxkaXYgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGh0bWwiIHN0eWxlPSJkaXNwbGF5OiBmbGV4OyBhbGlnbi1pdGVtczogdW5zYWZlIGNlbnRlcjsganVzdGlmeS1jb250ZW50OiB1bnNhZmUgY2VudGVyOyB3aWR0aDogMXB4OyBoZWlnaHQ6IDFweDsgcGFkZGluZy10b3A6IDEyOHB4OyBtYXJnaW4tbGVmdDogMTk4cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTFweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3dyYXA7Ij5PcHRpb25hbDwvZGl2PjwvZGl2PjwvZGl2PjwvZm9yZWlnbk9iamVjdD48dGV4dCB4PSIxOTgiIHk9IjEzMSIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMXB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5PcHRpb25hbDwvdGV4dD48L3N3aXRjaD48L2c+PHJlY3QgeD0iMTg0IiB5PSIxOTAiIHdpZHRoPSI5NSIgaGVpZ2h0PSI0MCIgZmlsbD0iI2QzZDZlYiIgc3Ryb2tlPSJub25lIiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDQ2cHg7IGhlaWdodDogMXB4OyBwYWRkaW5nLXRvcDogMTA1cHg7IG1hcmdpbi1sZWZ0OiA5M3B4OyI+PGRpdiBkYXRhLWRyYXdpby1jb2xvcnM9ImNvbG9yOiByZ2JhKDAsIDAsIDAsIDEpOyAiIHN0eWxlPSJib3gtc2l6aW5nOiBib3JkZXItYm94OyBmb250LXNpemU6IDBweDsgdGV4dC1hbGlnbjogY2VudGVyOyI+PGRpdiBzdHlsZT0iZGlzcGxheTogaW5saW5lLWJsb2NrOyBmb250LXNpemU6IDEycHg7IGZvbnQtZmFtaWx5OiAmcXVvdDtUaW1lcyBOZXcgUm9tYW4mcXVvdDs7IGNvbG9yOiByZ2IoMCwgMCwgMCk7IGxpbmUtaGVpZ2h0OiAxLjI7IHBvaW50ZXItZXZlbnRzOiBub25lOyB3aGl0ZS1zcGFjZTogbm9ybWFsOyBvdmVyZmxvdy13cmFwOiBub3JtYWw7Ij48c3BhbiBzdHlsZT0iZm9udC1zaXplOiA4cHgiPmZpbGw8L3NwYW4+PC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjExNiIgeT0iMTA5IiBmaWxsPSJyZ2JhKDAsIDAsIDAsIDEpIiBmb250LWZhbWlseT0iVGltZXMgTmV3IFJvbWFuIiBmb250LXNpemU9IjEycHgiIHRleHQtYW5jaG9yPSJtaWRkbGUiPmZpbGw8L3RleHQ+PC9zd2l0Y2g+PC9nPjxyZWN0IHg9IjM0OSIgeT0iMTkwIiB3aWR0aD0iOTUiIGhlaWdodD0iNDAiIGZpbGw9IiNmZmZmZmYiIHN0cm9rZT0ibm9uZSIgcG9pbnRlci1ldmVudHM9Im5vbmUiLz48ZyB0cmFuc2Zvcm09InRyYW5zbGF0ZSgtMC41IC0wLjUpc2NhbGUoMikiPjxzd2l0Y2g+PGZvcmVpZ25PYmplY3QgcG9pbnRlci1ldmVudHM9Im5vbmUiIHdpZHRoPSIxMDAlIiBoZWlnaHQ9IjEwMCUiIHJlcXVpcmVkRmVhdHVyZXM9Imh0dHA6Ly93d3cudzMub3JnL1RSL1NWRzExL2ZlYXR1cmUjRXh0ZW5zaWJpbGl0eSIgc3R5bGU9Im92ZXJmbG93OiB2aXNpYmxlOyB0ZXh0LWFsaWduOiBsZWZ0OyI+PGRpdiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMTk5OS94aHRtbCIgc3R5bGU9ImRpc3BsYXk6IGZsZXg7IGFsaWduLWl0ZW1zOiB1bnNhZmUgY2VudGVyOyBqdXN0aWZ5LWNvbnRlbnQ6IHVuc2FmZSBjZW50ZXI7IHdpZHRoOiA0NnB4OyBoZWlnaHQ6IDFweDsgcGFkZGluZy10b3A6IDEwNXB4OyBtYXJnaW4tbGVmdDogMTc2cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTJweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3JtYWw7IG92ZXJmbG93LXdyYXA6IG5vcm1hbDsiPjxzcGFuIHN0eWxlPSJmb250LXNpemU6IDhweCI+bm8gZmlsbDwvc3Bhbj48L2Rpdj48L2Rpdj48L2Rpdj48L2ZvcmVpZ25PYmplY3Q+PHRleHQgeD0iMTk4IiB5PSIxMDkiIGZpbGw9InJnYmEoMCwgMCwgMCwgMSkiIGZvbnQtZmFtaWx5PSJUaW1lcyBOZXcgUm9tYW4iIGZvbnQtc2l6ZT0iMTJweCIgdGV4dC1hbmNob3I9Im1pZGRsZSI+bm8gZmlsbDwvdGV4dD48L3N3aXRjaD48L2c+PGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTAuNSAtMC41KXNjYWxlKDIpIj48c3dpdGNoPjxmb3JlaWduT2JqZWN0IHBvaW50ZXItZXZlbnRzPSJub25lIiB3aWR0aD0iMTAwJSIgaGVpZ2h0PSIxMDAlIiByZXF1aXJlZEZlYXR1cmVzPSJodHRwOi8vd3d3LnczLm9yZy9UUi9TVkcxMS9mZWF0dXJlI0V4dGVuc2liaWxpdHkiIHN0eWxlPSJvdmVyZmxvdzogdmlzaWJsZTsgdGV4dC1hbGlnbjogbGVmdDsiPjxkaXYgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGh0bWwiIHN0eWxlPSJkaXNwbGF5OiBmbGV4OyBhbGlnbi1pdGVtczogdW5zYWZlIGNlbnRlcjsganVzdGlmeS1jb250ZW50OiB1bnNhZmUgY2VudGVyOyB3aWR0aDogMXB4OyBoZWlnaHQ6IDFweDsgcGFkZGluZy10b3A6IDI0cHg7IG1hcmdpbi1sZWZ0OiAzNXB4OyI+PGRpdiBkYXRhLWRyYXdpby1jb2xvcnM9ImNvbG9yOiByZ2JhKDAsIDAsIDAsIDEpOyAiIHN0eWxlPSJib3gtc2l6aW5nOiBib3JkZXItYm94OyBmb250LXNpemU6IDBweDsgdGV4dC1hbGlnbjogY2VudGVyOyI+PGRpdiBzdHlsZT0iZGlzcGxheTogaW5saW5lLWJsb2NrOyBmb250LXNpemU6IDEwcHg7IGZvbnQtZmFtaWx5OiAmcXVvdDtUaW1lcyBOZXcgUm9tYW4mcXVvdDs7IGNvbG9yOiByZ2IoMCwgMCwgMCk7IGxpbmUtaGVpZ2h0OiAxLjI7IHBvaW50ZXItZXZlbnRzOiBub25lOyBmb250LXdlaWdodDogYm9sZDsgd2hpdGUtc3BhY2U6IG5vd3JhcDsiPlJ1bm5pbmcgPGJyIHN0eWxlPSJmb250LXNpemU6IDEwcHgiIC8+U2Vzc2lvbjwvZGl2PjwvZGl2PjwvZGl2PjwvZm9yZWlnbk9iamVjdD48dGV4dCB4PSIzNSIgeT0iMjciIGZpbGw9InJnYmEoMCwgMCwgMCwgMSkiIGZvbnQtZmFtaWx5PSJUaW1lcyBOZXcgUm9tYW4iIGZvbnQtc2l6ZT0iMTBweCIgdGV4dC1hbmNob3I9Im1pZGRsZSIgZm9udC13ZWlnaHQ9ImJvbGQiPlJ1bm5pbmcuLi48L3RleHQ+PC9zd2l0Y2g+PC9nPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiAxMTNweDsgbWFyZ2luLWxlZnQ6IDM1cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTBweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IGZvbnQtd2VpZ2h0OiBib2xkOyB3aGl0ZS1zcGFjZTogbm93cmFwOyI+PGRpdiBzdHlsZT0iZm9udC1zaXplOiAxMHB4Ij48c3BhbiBzdHlsZT0iYmFja2dyb3VuZC1jb2xvcjogdHJhbnNwYXJlbnQgOyBmb250LXNpemU6IDEwcHgiPlJ1bm5pbmc8L3NwYW4+PC9kaXY+UmVxdWlyZW1lbnQ8L2Rpdj48L2Rpdj48L2Rpdj48L2ZvcmVpZ25PYmplY3Q+PHRleHQgeD0iMzUiIHk9IjExNiIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMHB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIiBmb250LXdlaWdodD0iYm9sZCI+UnVubmluZ1JlcXVpcmUuLi48L3RleHQ+PC9zd2l0Y2g+PC9nPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiA2OHB4OyBtYXJnaW4tbGVmdDogMzVweDsiPjxkaXYgZGF0YS1kcmF3aW8tY29sb3JzPSJjb2xvcjogcmdiYSgwLCAwLCAwLCAxKTsgIiBzdHlsZT0iYm94LXNpemluZzogYm9yZGVyLWJveDsgZm9udC1zaXplOiAwcHg7IHRleHQtYWxpZ246IGNlbnRlcjsiPjxkaXYgc3R5bGU9ImRpc3BsYXk6IGlubGluZS1ibG9jazsgZm9udC1zaXplOiAxMHB4OyBmb250LWZhbWlseTogJnF1b3Q7VGltZXMgTmV3IFJvbWFuJnF1b3Q7OyBjb2xvcjogcmdiKDAsIDAsIDApOyBsaW5lLWhlaWdodDogMS4yOyBwb2ludGVyLWV2ZW50czogbm9uZTsgZm9udC13ZWlnaHQ6IGJvbGQ7IHdoaXRlLXNwYWNlOiBub3dyYXA7Ij5TdGVwIDxiciBzdHlsZT0iZm9udC1zaXplOiAxMHB4IiAvPkNsYXNzPC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjM1IiB5PSI3MSIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMHB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIiBmb250LXdlaWdodD0iYm9sZCI+U3RlcC4uLjwvdGV4dD48L3N3aXRjaD48L2c+PC9nPjxzd2l0Y2g+PGcgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5Ii8+PGEgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoMCwtNSkiIHhsaW5rOmhyZWY9Imh0dHBzOi8vd3d3LmRpYWdyYW1zLm5ldC9kb2MvZmFxL3N2Zy1leHBvcnQtdGV4dC1wcm9ibGVtcyIgdGFyZ2V0PSJfYmxhbmsiPjx0ZXh0IHRleHQtYW5jaG9yPSJtaWRkbGUiIGZvbnQtc2l6ZT0iMTBweCIgeD0iNTAlIiB5PSIxMDAlIj5WaWV3ZXIgZG9lcyBub3Qgc3VwcG9ydCBmdWxsIFNWRyAxLjE8L3RleHQ+PC9hPjwvc3dpdGNoPjwvc3ZnPg=="},"evals":[],"jsHooks":[]}</script>

## Checking workflow status

To check the summary of the workflow, we can use:

``` r
sal
statusWF(sal)
```

## Technical report

*`systemPipeR`* compiles all the workflow execution logs in one central location,
making it easier to check any standard output (`stdout`) or standard error
(`stderr`) for any command-line tools used on the workflow or the R code stdout.

``` r
sal <- renderLogs(sal)
```

## Scientific report

*`systemPipeR`* auto-generates scientific analysis reports in HTML format.

``` r
sal <- renderReport(sal)
```

## Funding

This project is funded by NSF award [ABI-1661152](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1661152).

## References

<div id="refs" class="references">

<div id="ref-H_Backman2016-bt">

H Backman, Tyler W, and Thomas Girke. 2016. “systemPipeR: NGS workflow and report generation environment.” *BMC Bioinformatics* 17 (1): 388. <https://doi.org/10.1186/s12859-016-1241-0>.

</div>

<div id="ref-Howard2013-fq">

Howard, Brian E, Qiwen Hu, Ahmet Can Babaoglu, Manan Chandra, Monica Borghi, Xiaoping Tan, Luyan He, et al. 2013. “High-Throughput RNA Sequencing of Pseudomonas-Infected Arabidopsis Reveals Hidden Transcriptome Complexity and Novel Splice Variants.” *PLoS One* 8 (10): e74183. <https://doi.org/10.1371/journal.pone.0074183>.

</div>

<div id="ref-Kim2015-ve">

Kim, Daehwan, Ben Langmead, and Steven L Salzberg. 2015. “HISAT: A Fast Spliced Aligner with Low Memory Requirements.” *Nat. Methods* 12 (4): 357–60.

</div>

<div id="ref-Lawrence2013-kt">

Lawrence, Michael, Wolfgang Huber, Hervé Pagès, Patrick Aboyoun, Marc Carlson, Robert Gentleman, Martin T Morgan, and Vincent J Carey. 2013. “Software for Computing and Annotating Genomic Ranges.” *PLoS Comput. Biol.* 9 (8): e1003118. <https://doi.org/10.1371/journal.pcbi.1003118>.

</div>

<div id="ref-Robinson2010-uk">

Robinson, M D, D J McCarthy, and G K Smyth. 2010. “EdgeR: A Bioconductor Package for Differential Expression Analysis of Digital Gene Expression Data.” *Bioinformatics* 26 (1): 139–40. <https://doi.org/10.1093/bioinformatics/btp616>.

</div>

</div>
