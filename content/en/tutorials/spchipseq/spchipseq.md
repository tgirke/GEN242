---
title: "ChIP-Seq Workflow Template" 
author: "Author: First Last"
date: "Last update: 05 June, 2023" 
output:
  BiocStyle::html_document:
    toc_float: true
    code_folding: show
package: systemPipeR
vignette: |
  %\VignetteEncoding{UTF-8}
  %\VignetteIndexEntry{WF: ChIP-Seq Workflow Template}
  %\VignetteEngine{knitr::rmarkdown}
fontsize: 14pt
bibliography: bibtex.bib
weight: 9
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

<link id="plotwf_legend-1-attachment" rel="attachment" href="spchipseq_files/plotwf_legend/plotwf_legend.svg"/>

<script src="/rmarkdown-libs/plotwf-binding/plotwf.js"></script>

<!--
- Compile from command-line
Rscript -e "rmarkdown::render('spchipseq.Rmd', c('BiocStyle::html_document'), clean=F); knitr::knit('spchipseq.Rmd', tangle=TRUE)"
-->

<style type="text/css">
pre code {
white-space: pre !important;
overflow-x: scroll !important;
word-break: keep-all !important;
word-wrap: initial !important;
}
</style>

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

Source code downloads:    
\[ [.Rmd](https://raw.githubusercontent.com/tgirke/GEN242/main/static/custom/spWFtemplates/spchipseq.Rmd) \]    
\[ [.R](https://raw.githubusercontent.com/tgirke/GEN242/main/content/en/tutorials/spchipseq/spchipseq.R) \]

</div>

## Introduction

The following analyzes the ChIP-Seq data from Kaufman et al. (2010) using
for peak calling MACS2 where the uninduced sample serves as input (reference). Prior to running
this analysis the corresponding FASTQ files need to be downloaded following the instructions
[here](https://girke.bioinformatics.ucr.edu/GEN242/assignments/projects/project_data/).

For learning purposes one can use the much smaller toy data set and ChIP-Seq `Rmd`
workflow instance provided by the `systemPipeRdata` package. This toy workflow instance
can be conveniently obtained by running `genWorkenvir(workflow='chipseq')` (see
below). The FASTQ data used for this toy instance are the same as for the
RNA-Seq workflow.

For the analysis of the Kaufman et al. (2010) data set (see download [here](https://girke.bioinformatics.ucr.edu/GEN242/assignments/projects/project_data/))
users want to use the `Rmd` instance linked from the top right corner of this page.
Additional detail about this is provided [here](https://girke.bioinformatics.ucr.edu/GEN242/assignments/projects/project_data/#workflow-rmd-file).

### Experimental design

Typically, users want to specify here all information relevant for the
analysis of their NGS study. This includes detailed descriptions of
FASTQ files, experimental design, reference genome, gene annotations,
etc.

### Workflow environment

<font color="red">NOTE: this section</font> describes how to set up the proper
environment (directory structure) for running `systemPipeR` workflows. After
mastering this task the workflow run instructions <font color="red">can be deleted</font>
since they are not expected to be included in a final HTML/PDF report of a workflow.

1.  If a remote system or cluster is used, then users need to log in to the
    remote system first. The following applies to an HPC cluster (*e.g.* HPCC
    cluster).
    
    A terminal application needs to be used to log in to a user’s cluster account. Next, one
    can open an interactive session on a computer node with `srun --x11`. More details about
    argument settings for `srun` are available in this [HPCC
    manual](http://hpcc.ucr.edu/manuals_linux-cluster_jobs.html#partitions) or
    the HPCC section of this website
    [here](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/linux/linux/#job-submission-with-sbatch).
    Next, load the R version required for running the workflow with `module load`. Sometimes it may be necessary to
    first unload an active software version before loading another version, *e.g.* `module unload R`.

<!-- end list -->

``` bash
srun --x11 --partition=gen242 --mem=20gb --cpus-per-task 8 --ntasks 1 --time 20:00:00 --pty bash -l
module unload R; module load R/4.2.2
```

If the above `srun` command denies access to the gen242 parition, please try adding `--account gen242`.

2.  Load a workflow template with the `genWorkenvir` function. This can be done from the command-line or from within R.
    However, only one of the two options needs to be used.

From command-line

``` bash
$ Rscript -e "systemPipeRdata::genWorkenvir(workflow='chipseq')"
$ cd chipseq
```

From R

``` r
library(systemPipeRdata)
genWorkenvir(workflow = "chipseq")
setwd("chipseq")
```

3.  Optional: if the user wishes to use another `Rmd` file than the template
    instance provided by the `genWorkenvir` function, then it can be copied or downloaded
    into the root directory of the workflow environment (*e.g.* with `cp` or `wget`).

4.  Now one can open from the root directory of the workflow the corresponding
    R Markdown script (*e.g.* systemPipeChIPseq.Rmd) using an R IDE, such as *nvim-r*,
    *ESS* or RStudio. Subsequently, the workflow can be run as outlined below. For
    learning purposes it is recommended to run workflows for the first time interactively.
    Once all workflow steps are understood and possibly modified to custom needs,
    one can run the workflow from start to finish with a single command using `runWF()`.

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
[here](http://www.bioconductor.org/packages/release/bioc/vignettes/systemPipeR/inst/doc/systemPipeR.html#25_structure_of_targets_file). More details about the new parameter files from `systemPipeR` can be found [here](http://www.bioconductor.org/packages/release/bioc/vignettes/systemPipeR/inst/doc/systemPipeR.html#26_structure_of_the_new_param_files_and_construct_sysargs2_container).

### Import custom functions

Custom functions for the challenge projects can be imported with the source
command from a local R script (here [challengeProject\_Fct.R](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/spchipseq/challengeProject_Fct.R)). Skip this step if such a script is not available. Alternatively,
these functions can be loaded from a custom R package.

``` r
source("challengeProject_Fct.R")
```

### Experiment definition provided by `targets` file

The `targets` file defines all FASTQ files and sample comparisons of the analysis workflow.
If needed the tab separated (TSV) version of this file can be downloaded from [here](https://github.com/tgirke/GEN242/tree/main/content/en/assignments/Projects/targets_files)
and the corresponding Google Sheet is [here](https://docs.google.com/spreadsheets/d/1w9V3JDOsXR8qW_qNqoXO0E0UR_Oev4k7-o2y6NMDTt8/edit#gid=472150521).

``` r
targetspath <- "targets_chipseq.txt"
targets <- read.delim(targetspath, comment.char = "#")
knitr::kable(targets)
```

| FileName                     | SampleName | Factor | SampleLong       | Experiment | Date      | SampleReference |
| :--------------------------- | :--------- | :----- | :--------------- | ---------: | :-------- | :-------------- |
| ./data/SRR038845\_1.fastq.gz | AP1\_1     | AP1    | APETALA1 Induced |          1 | 23-Mar-12 |                 |
| ./data/SRR038846\_1.fastq.gz | AP1\_2A    | AP1    | APETALA1 Induced |          1 | 23-Mar-12 |                 |
| ./data/SRR038847\_1.fastq.gz | AP1\_2B    | AP1    | APETALA1 Induced |          1 | 23-Mar-12 |                 |
| ./data/SRR038848\_1.fastq.gz | C\_1A      | C      | Control Mock     |          1 | 23-Mar-12 | AP1\_1          |
| ./data/SRR038849\_1.fastq.gz | C\_1B      | C      | Control Mock     |          1 | 23-Mar-12 | AP1\_1          |
| ./data/SRR038850\_1.fastq.gz | C\_2A      | C      | Control Mock     |          1 | 23-Mar-12 | AP1\_2A         |
| ./data/SRR038851\_1.fastq.gz | C\_2B      | C      | Control Mock     |          1 | 23-Mar-12 | AP1\_2B         |

## Workflow steps

This tutorial will demonstrate how to either build and run the workflow automatically,
or in an interactive mode by appending each step with the `appendStep` method. In both
cases the `SYSargsList` object will be populated with the instructions for running each
workflow step, while supporting both command-line steps as well as line-wise R commands
defined in the corresponding code chunks of this or any `Rmd` file that has been properly
formatted.

To create a workflow within *`systemPipeR`*, we can start by defining an empty
`SYSargsList` container. When restarting an existing workflow one
can set `resume=TRUE` under the `SPRproject()` function call.

``` r
library(systemPipeR)
sal <- SPRproject()
```

    ## Creating directory:  /home/tgirke/tmp/GEN242/content/en/tutorials/spchipseq/data 
    ## Creating directory '/home/tgirke/tmp/GEN242/content/en/tutorials/spchipseq/.SPRproject'
    ## Creating file '/home/tgirke/tmp/GEN242/content/en/tutorials/spchipseq/.SPRproject/SYSargsList.yml'

``` r
sal
```

    ## Instance of 'SYSargsList': 
    ##  No workflow steps added

Next, the `importWF` function will load the entire workflow into the `SYSargsList` object
(here `sal`). Subsequently, the `runWF()` function will run the workflow from start to finish.
If needed, specific workflow steps can be executed by assigning their corresponding
position numbers within the workflow to the `steps` argument (see `?runWF`). After completion of the workflow
one can render a scientific analysis report in HTML format with the `renderReport()`
function that uses R Markdown internally.

``` r
sal <- importWF(sal, file_path = "systemPipeChIPseq.Rmd")  ## Import all the Workflow steps
sal
sal <- runWF(sal)  # Runs workflow
sal <- renderReport(sal)  # Renders report
rmarkdown::render("systemPipeChIPseq.Rmd", clean = TRUE, output_format = "BiocStyle::html_document")  # Alternative report rendering
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
appendStep(sal) <- SYSargsList(step_name = "preprocessing", targets = "targets_chipseq.txt",
    dir = TRUE, wf_file = "preprocessReads/preprocessReads-se.cwl",
    input_file = "preprocessReads/preprocessReads-se.yml", dir_path = "param/cwl",
    inputvars = c(FileName = "_FASTQ_PATH1_", SampleName = "_SampleName_"),
    dependency = c("load_SPR"))
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

Figure 1: FASTQ quality report for 7 samples.

</div>

</br>

### Alignments

#### Read mapping with `Bowtie2`

The NGS reads of this project will be aligned with `Bowtie2` against the
reference genome sequence (Langmead and Salzberg 2012). The parameter settings of the
Bowtie2 index are defined in the `bowtie2-index.cwl` and `bowtie2-index.yml` files.

Building the index:

``` r
appendStep(sal) <- SYSargsList(step_name = "bowtie2_index", dir = FALSE,
    targets = NULL, wf_file = "bowtie2/bowtie2-index.cwl", input_file = "bowtie2/bowtie2-index.yml",
    dir_path = "param/cwl", inputvars = NULL, dependency = c("preprocessing"))
```

The parameter settings of the aligner are defined in the `workflow_bowtie2-se.cwl`
and `workflow_bowtie2-se.yml` files. The following shows how to construct the
corresponding *SYSargsList* object.

In ChIP-Seq experiments it is usually more appropriate to eliminate reads mapping
to multiple locations. To achieve this, users want to remove the argument setting
`-k 50 non-deterministic` in the configuration files.

``` r
appendStep(sal) <- SYSargsList(step_name = "bowtie2_alignment",
    dir = TRUE, targets = "targets_chipseq.txt", wf_file = "workflow-bowtie2/workflow_bowtie2-se.cwl",
    input_file = "workflow-bowtie2/workflow_bowtie2-se.yml",
    dir_path = "param/cwl", inputvars = c(FileName = "_FASTQ_PATH1_",
        SampleName = "_SampleName_"), dependency = c("bowtie2_index"))
```

To double-check the command line for each sample, please use the following:

``` r
cmdlist(sal, step = "bowtie2_alignment", targets = 1)
```

    ## $bowtie2_alignment
    ## $bowtie2_alignment$AP1_1
    ## $bowtie2_alignment$AP1_1$bowtie2
    ## [1] "bowtie2 -S ./results/AP1_1.sam  -x ./data/tair10.fasta  -k 50  --non-deterministic  -U ./data/SRR038845_1.fastq.gz -p 4"
    ## 
    ## $bowtie2_alignment$AP1_1$`samtools-view`
    ## [1] "samtools view -bS -o ./results/AP1_1.bam  ./results/AP1_1.sam "
    ## 
    ## $bowtie2_alignment$AP1_1$`samtools-sort`
    ## [1] "samtools sort -o ./results/AP1_1.sorted.bam  ./results/AP1_1.bam  -@ 4"
    ## 
    ## $bowtie2_alignment$AP1_1$`samtools-index`
    ## [1] "samtools index -b results/AP1_1.sorted.bam  results/AP1_1.sorted.bam.bai  ./results/AP1_1.sorted.bam "

#### Read and alignment stats

The following provides an overview of the number of reads in each sample
and how many of them aligned to the reference.

``` r
appendStep(sal) <- LineWise(code = {
    fqpaths <- getColumn(sal, step = "bowtie2_alignment", "targetsWF",
        column = "FileName")
    bampaths <- getColumn(sal, step = "bowtie2_alignment", "outfiles",
        column = "samtools_sort_bam")
    read_statsDF <- alignStats(args = bampaths, fqpaths = fqpaths,
        pairEnd = TRUE)
    write.table(read_statsDF, "results/alignStats.xls", row.names = FALSE,
        quote = FALSE, sep = "\t")
}, step_name = "align_stats", dependency = "bowtie2_alignment")
```

### Create symbolic links for viewing BAM files in IGV

The `symLink2bam` function creates symbolic links to view the BAM alignment
files in a genome browser such as IGV without moving these large files to a
local system. The corresponding URLs are written to a file with a path
specified under `urlfile`, here `IGVurl.txt`. Please replace the directory
and the user name. The following parameter settings will create a
subdirectory under `~/.html` called `somedir` of the user account. The user
name under `urlbase`, here `<username>`, needs to be changed to the corresponding
user name of the person running this function.

``` r
appendStep(sal) <- LineWise(code = {
    bampaths <- getColumn(sal, step = "bowtie2_alignment", "outfiles",
        column = "samtools_sort_bam")
    bampaths <- setNames(normalizePath(bampaths), names(bampaths))
    symLink2bam(sysargs = bampaths, htmldir = c("~/.html/", "somedir/"),
        urlbase = "http://cluster.hpcc.ucr.edu/~<username>/",
        urlfile = "./results/IGVurl.txt")
}, step_name = "bam_IGV", dependency = "bowtie2_alignment", run_step = "optional")
```

### Peak calling with MACS2

### Merge BAM files of replicates prior to peak calling

Merging BAM files of technical and/or biological replicates can improve
the sensitivity of the peak calling by increasing the depth of read
coverage. The `mergeBamByFactor` function merges BAM files based on grouping information
specified by a `factor`, here the `Factor` column of the imported targets file.
It also returns an updated `targets` object containing the paths to the
merged BAM files as well as to any unmerged files without replicates.
The updated `targets` object can be used to update the `SYSargsList` object.

This step can be skipped if merging of BAM files is not desired.

``` r
appendStep(sal) <- LineWise(code = {
    bampaths <- getColumn(sal, step = "bowtie2_alignment", "outfiles",
        column = "samtools_sort_bam")
    merge_bams <- mergeBamByFactor(args = bampaths, targetsDF = targetsWF(sal)[["bowtie2_alignment"]],
        overwrite = TRUE)
    updateColumn(sal, step = "merge_bams", position = "targetsWF") <- merge_bams
    writeTargets(sal, step = "merge_bams", file = "targets_merge_bams.txt",
        overwrite = TRUE)
}, step_name = "merge_bams", dependency = "bowtie2_alignment")
```

#### Peak calling with input/reference sample

MACS2 can perform peak calling on ChIP-Seq data with and without input
samples (Zhang et al. 2008).

The following performs peak calling with input sample. The input sample
can be most conveniently specified in the `SampleReference` column of the
initial `targets` file. The `writeTargetsRef` function uses this
information to create a `targets` file intermediate for running MACS2
with the corresponding input sample(s).

``` r
appendStep(sal) <- LineWise(code = {
    writeTargetsRef(infile = "targets_merge_bams.txt", outfile = "targets_bam_ref.txt",
        silent = FALSE, overwrite = TRUE)
}, step_name = "writeTargetsRef", dependency = "merge_bams")
```

``` r
appendStep(sal) <- SYSargsList(step_name = "call_peaks_macs_withref",
    targets = "targets_bam_ref.txt", wf_file = "MACS2/macs2-input.cwl",
    input_file = "MACS2/macs2-input.yml", dir_path = "param/cwl",
    inputvars = c(FileName1 = "_FASTQ_PATH1_", FileName2 = "_FASTQ_PATH2_",
        SampleReference = "_SampleName_"), id = "SampleReference",
    dependency = c("writeTargetsRef"))
```

The peak calling results from MACS2 are written for each sample to
separate files in the `results/call_peaks_macs_withref` directory. They are
named after the corresponding files with extensions used by MACS2.

### Annotate peaks with genomic context

#### Annotation with `ChIPseeker` package

The following annotates the identified peaks with genomic context information
using the `ChIPseeker` package (Yu, Wang, and He 2015).

``` r
appendStep(sal) <- LineWise(code = {
    library(ChIPseeker)
    library(GenomicFeatures)
    peaks_files <- getColumn(sal, step = "call_peaks_macs_withref",
        "outfiles", column = "peaks_xls")
    txdb <- suppressWarnings(makeTxDbFromGFF(file = "data/tair10.gff",
        format = "gff", dataSource = "TAIR", organism = "Arabidopsis thaliana"))
    for (i in seq(along = peaks_files)) {
        peakAnno <- annotatePeak(peaks_files[i], TxDb = txdb,
            verbose = FALSE)
        df <- as.data.frame(peakAnno)
        outpaths <- paste0("./results/", names(peaks_files),
            "_ChIPseeker_annotated.xls")
        names(outpaths) <- names(peaks_files)
        write.table(df, outpaths[i], quote = FALSE, row.names = FALSE,
            sep = "\t")
    }
    updateColumn(sal, step = "annotation_ChIPseeker", position = "outfiles") <- data.frame(outpaths)
}, step_name = "annotation_ChIPseeker", dependency = "call_peaks_macs_withref")
```

The peak annotation results are written to the `results` directory.

Summary plots provided by the `ChIPseeker` package. Here applied only to one sample
for demonstration purposes.

``` r
appendStep(sal) <- LineWise(code = {
    peaks_files <- getColumn(sal, step = "call_peaks_macs_withref",
        "outfiles", column = "peaks_xls")
    peak <- readPeakFile(peaks_files[1])
    pdf("results/peakscoverage.pdf")
    covplot(peak, weightCol = "X.log10.pvalue.")
    dev.off()
    pdf("results/peaksHeatmap.pdf")
    peakHeatmap(peaks_files[1], TxDb = txdb, upstream = 1000,
        downstream = 1000, color = "red")
    dev.off()
    pdf("results/peaksProfile.pdf")
    plotAvgProf2(peaks_files[1], TxDb = txdb, upstream = 1000,
        downstream = 1000, xlab = "Genomic Region (5'->3')",
        ylab = "Read Count Frequency", conf = 0.05)
    dev.off()
}, step_name = "ChIPseeker_plots", dependency = "annotation_ChIPseeker")
```

### Count reads overlapping peaks

The `countRangeset` function is a convenience wrapper to perform read counting
iteratively over serveral range sets, here peak range sets. Internally,
the read counting is performed with the `summarizeOverlaps` function from the
`GenomicAlignments` package. The resulting count tables are directly saved to
files, one for each peak set.

``` r
appendStep(sal) <- LineWise(code = {
    library(GenomicRanges)
    bam_files <- getColumn(sal, step = "bowtie2_alignment", "outfiles",
        column = "samtools_sort_bam")
    args <- getColumn(sal, step = "call_peaks_macs_withref",
        "outfiles", column = "peaks_xls")
    outfiles <- paste0("./results/", names(args), "_countDF.xls")
    bfl <- BamFileList(bam_files, yieldSize = 50000, index = character())
    countDFnames <- countRangeset(bfl, args, outfiles, mode = "Union",
        ignore.strand = TRUE)
    updateColumn(sal, step = "count_peak_ranges", position = "outfiles") <- data.frame(countDFnames)
}, step_name = "count_peak_ranges", dependency = "call_peaks_macs_withref",
    )
```

Shows count table generated in previous step (`results/AP1_1_countDF.xls`).
To avoid slowdowns of the load time of this page, ony 200 rows of the source
table are imported into the below `datatable` view .

``` r
countDF <- read.delim("results/AP1_1_countDF.xls")[1:200, ]
colnames(countDF)[1] <- "PeakIDs"
library(DT)
datatable(countDF)
```

<div id="htmlwidget-1" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-1">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120","121","122","123","124","125","126","127","128","129","130","131","132","133","134","135","136","137","138","139","140","141","142","143","144","145","146","147","148","149","150","151","152","153","154","155","156","157","158","159","160","161","162","163","164","165","166","167","168","169","170","171","172","173","174","175","176","177","178","179","180","181","182","183","184","185","186","187","188","189","190","191","192","193","194","195","196","197","198","199","200"],["Chr1_10000123-10001504","Chr1_10015050-10015329","Chr1_10044084-10044292","Chr1_10053060-10053403","Chr1_10065376-10065801","Chr1_10077918-10078299","Chr1_10111557-10111830","Chr1_10137133-10137852","Chr1_10139865-10140426","Chr1_10141346-10141599","Chr1_10142995-10143536","Chr1_10143619-10144022","Chr1_10145791-10146202","Chr1_10146333-10147178","Chr1_10186868-10187102","Chr1_10187659-10187893","Chr1_10203888-10204202","Chr1_10219688-10219985","Chr1_10228577-10229094","Chr1_10229262-10229710","Chr1_10242962-10243288","Chr1_10244590-10244954","Chr1_10247304-10247528","Chr1_10251234-10251447","Chr1_10263957-10264248","Chr1_102694-103029","Chr1_10294478-10294800","Chr1_10295972-10296231","Chr1_10353709-10354318","Chr1_10364919-10365439","Chr1_10369960-10370330","Chr1_10375634-10375976","Chr1_10376847-10377070","Chr1_10489119-10489357","Chr1_10489629-10489865","Chr1_1049818-1050057","Chr1_10504608-10505267","Chr1_10505986-10506492","Chr1_10519205-10519487","Chr1_10520174-10520649","Chr1_10520971-10521455","Chr1_10530766-10530996","Chr1_10535179-10535644","Chr1_10537198-10537646","Chr1_10538872-10539144","Chr1_10541237-10541831","Chr1_10542202-10542449","Chr1_10553478-10553874","Chr1_10562570-10563068","Chr1_10563816-10565000","Chr1_10565793-10566457","Chr1_10567369-10567945","Chr1_10569384-10569620","Chr1_105715-106037","Chr1_10581383-10581685","Chr1_10622431-10622757","Chr1_10631994-10632273","Chr1_10649137-10649367","Chr1_10668660-10668958","Chr1_10669009-10669489","Chr1_10676832-10677414","Chr1_10678636-10679743","Chr1_10679884-10680612","Chr1_10683747-10684382","Chr1_10684727-10685161","Chr1_10685376-10685989","Chr1_10691341-10691562","Chr1_10691767-10692130","Chr1_10692583-10693433","Chr1_10693954-10694617","Chr1_10695463-10696565","Chr1_10701975-10702395","Chr1_10704093-10704316","Chr1_10713567-10714255","Chr1_10715301-10715933","Chr1_1072348-1072586","Chr1_10726169-10726568","Chr1_1072802-1073152","Chr1_10737905-10738222","Chr1_10793236-10793520","Chr1_10794884-10795869","Chr1_10799998-10800389","Chr1_10800942-10801268","Chr1_10801874-10802349","Chr1_10803701-10803926","Chr1_1084359-1084961","Chr1_10845312-10845874","Chr1_10864966-10865601","Chr1_10871974-10872308","Chr1_10874129-10874328","Chr1_10874476-10875185","Chr1_10891639-10891997","Chr1_10909961-10910163","Chr1_10953463-10953804","Chr1_10954240-10954628","Chr1_10956271-10956623","Chr1_10958841-10959244","Chr1_10959439-10960076","Chr1_10972812-10973048","Chr1_1103200-1103416","Chr1_11033716-11034097","Chr1_11059998-11060940","Chr1_11061272-11061800","Chr1_11062322-11062750","Chr1_11067587-11068014","Chr1_11071124-11071500","Chr1_11084496-11084695","Chr1_11128-11443","Chr1_11155158-11155410","Chr1_11196849-11197415","Chr1_11200138-11200385","Chr1_11210773-11211568","Chr1_11212068-11212348","Chr1_11212560-11212978","Chr1_1121262-1121569","Chr1_11222722-11222952","Chr1_11227689-11228070","Chr1_11231069-11231415","Chr1_1129384-1129772","Chr1_1134217-1134971","Chr1_11343316-11343530","Chr1_11350683-11350893","Chr1_11374190-11375132","Chr1_11377612-11378008","Chr1_11385968-11386241","Chr1_11386635-11386872","Chr1_11388323-11388599","Chr1_11421622-11422073","Chr1_11428546-11428802","Chr1_11429167-11429366","Chr1_11429941-11430188","Chr1_11430639-11430926","Chr1_11441925-11442147","Chr1_11446002-11446822","Chr1_11447290-11447887","Chr1_115641-116419","Chr1_1159091-1159510","Chr1_11594952-11595190","Chr1_11598295-11598860","Chr1_11600614-11600983","Chr1_11601139-11601823","Chr1_11605622-11605926","Chr1_11620054-11620600","Chr1_11621059-11621447","Chr1_11621523-11623089","Chr1_11624914-11625625","Chr1_11626261-11626946","Chr1_11628379-11628622","Chr1_11628853-11629244","Chr1_11630475-11631803","Chr1_11662337-11662575","Chr1_1172920-1173129","Chr1_1173195-1173526","Chr1_11763531-11763800","Chr1_11770202-11770807","Chr1_11771164-11771564","Chr1_11798751-11798970","Chr1_11801046-11801322","Chr1_11804660-11805026","Chr1_11807363-11807831","Chr1_11808375-11808791","Chr1_11809718-11810175","Chr1_11829516-11829982","Chr1_11830901-11831419","Chr1_11874425-11874646","Chr1_1188104-1188901","Chr1_11925063-11925373","Chr1_11927715-11928115","Chr1_11928466-11928909","Chr1_11978528-11979276","Chr1_1197874-1198387","Chr1_12004740-12004982","Chr1_12005493-12006453","Chr1_12026542-12027136","Chr1_12050647-12051422","Chr1_12051618-12051822","Chr1_12054702-12055258","Chr1_12055957-12056342","Chr1_12056485-12057547","Chr1_12059598-12059861","Chr1_12064129-12064328","Chr1_12070355-12070698","Chr1_12123473-12123764","Chr1_12132602-12132955","Chr1_12177510-12177803","Chr1_12180624-12180825","Chr1_12254744-12255260","Chr1_1233120-1233354","Chr1_1235877-1236728","Chr1_12458546-12458848","Chr1_12474908-12475186","Chr1_128465-128823","Chr1_12916443-12916688","Chr1_12986436-12986729","Chr1_12986791-12987139","Chr1_13029530-13030312","Chr1_1303036-1303264","Chr1_130507-131207","Chr1_1306355-1306663","Chr1_1307140-1307396"],[384,19,54,31,89,82,36,95,135,29,44,44,59,134,27,14,50,49,193,75,36,29,22,24,33,21,19,16,108,125,38,17,20,20,28,36,135,53,23,168,101,19,67,26,31,55,40,39,40,290,229,40,21,39,28,13,36,16,30,40,96,271,81,63,48,106,18,40,207,86,121,85,19,158,100,38,72,40,30,16,277,51,54,97,22,88,120,77,19,18,137,64,23,44,53,47,33,150,22,29,39,107,161,49,94,23,8,60,53,91,27,58,32,69,33,26,26,47,49,222,15,32,191,35,17,26,43,71,17,16,13,39,6,120,80,93,116,17,86,41,178,36,125,48,178,59,92,23,41,284,80,17,38,23,104,104,12,13,35,72,80,82,49,98,15,100,16,36,32,251,72,25,143,127,177,16,94,37,219,22,12,80,55,37,30,19,102,9,268,26,13,74,31,23,37,304,30,93,40,16],[571,50,31,71,122,115,59,201,191,38,102,72,120,296,32,44,59,56,357,129,54,62,36,35,47,64,66,43,213,149,69,48,33,45,37,19,281,95,63,301,180,39,109,83,58,86,50,80,100,625,371,100,35,69,46,78,46,44,68,106,148,447,152,152,123,205,38,71,347,159,264,96,31,196,144,37,150,42,47,37,434,93,79,135,37,132,186,137,71,38,236,85,33,54,86,63,76,243,53,29,114,319,248,135,277,73,30,92,37,165,62,157,91,138,61,34,81,82,63,330,31,32,295,78,58,48,52,131,46,26,43,54,36,167,112,154,128,31,89,49,199,60,165,90,399,162,177,29,84,485,141,32,70,41,152,139,39,46,75,127,120,169,157,141,35,199,74,82,114,343,164,44,237,228,284,28,164,61,433,37,43,99,44,75,69,33,160,34,285,59,62,153,37,52,45,481,32,281,53,38],[282,26,19,32,89,46,21,90,95,26,58,30,58,149,25,20,35,24,158,80,26,22,22,21,23,26,49,17,96,83,24,28,13,18,16,15,107,40,32,152,74,21,49,36,27,36,19,31,45,348,151,54,14,31,15,32,31,18,21,41,83,222,82,62,40,102,10,30,143,98,119,38,12,92,66,24,94,21,20,27,189,41,52,70,19,41,94,59,33,19,132,26,11,24,29,25,39,119,13,14,53,155,110,55,110,50,14,42,15,77,21,91,27,58,37,24,42,44,31,151,15,25,142,32,22,20,18,67,28,15,22,21,19,75,46,58,43,14,45,27,90,25,51,44,191,76,83,19,48,248,67,21,35,16,68,81,14,19,42,57,78,67,74,58,9,76,24,23,45,169,72,28,117,96,141,14,67,29,195,20,17,57,25,38,29,13,78,26,118,16,25,71,15,28,16,200,16,129,27,23],[56,1,4,6,1,12,3,15,7,5,10,34,1,1,2,0,4,13,1,17,1,0,0,0,4,3,1,1,15,15,10,9,2,5,3,1,6,19,0,3,13,6,11,0,8,2,18,4,5,9,7,7,2,10,0,4,5,0,10,7,0,11,3,5,1,17,0,10,12,0,25,3,0,4,19,10,4,15,5,0,12,2,1,2,2,60,16,2,5,14,2,3,2,12,12,11,3,6,4,1,2,4,9,0,24,3,14,1,7,0,6,8,2,9,4,20,1,10,6,1,1,7,10,34,7,4,3,4,6,2,0,22,1,6,13,1,23,0,33,9,1,17,61,2,21,1,19,0,9,3,5,5,10,5,0,16,3,34,6,14,0,0,5,25,0,8,11,1,1,43,9,0,10,15,10,0,22,16,14,4,0,1,14,2,19,21,3,1,34,5,11,67,29,5,15,8,1,9,1,3],[43,3,3,6,6,4,8,8,8,8,10,12,1,4,3,1,9,5,2,5,7,1,6,1,6,3,4,1,9,10,8,6,2,8,1,0,6,22,2,3,18,4,10,3,8,7,9,8,5,32,17,7,7,10,0,3,4,0,8,3,0,21,23,9,1,20,1,11,12,1,33,7,2,9,13,8,9,16,3,2,9,4,8,8,2,19,14,11,5,5,17,1,0,7,6,5,5,9,2,4,6,17,12,3,16,1,8,8,5,1,7,2,3,3,9,16,0,7,7,19,9,2,23,19,4,0,7,14,2,3,0,15,0,6,24,4,19,2,10,7,8,12,17,0,39,14,11,4,3,6,22,4,12,7,1,11,5,10,11,15,0,1,9,17,2,13,11,0,4,26,12,0,21,18,17,0,15,5,19,1,1,5,6,8,11,6,1,1,22,9,5,27,7,1,9,14,0,6,1,7],[104,12,14,25,18,10,17,20,20,5,21,11,11,15,8,10,13,16,32,16,12,10,9,7,18,11,7,4,21,27,12,8,2,6,10,4,22,16,13,30,11,17,13,15,13,16,6,14,14,33,36,11,7,2,6,11,12,12,16,19,13,41,21,20,10,25,6,22,30,15,33,15,14,17,19,9,35,10,14,9,32,22,8,8,8,15,18,17,20,14,22,15,5,7,13,5,12,11,6,7,14,41,22,20,46,8,11,14,7,16,5,9,12,10,9,5,12,19,9,31,8,12,48,15,11,4,5,36,4,4,12,5,13,20,20,18,28,6,13,15,31,10,12,7,41,23,26,5,22,36,43,5,13,2,28,29,7,5,16,18,18,18,24,12,10,30,9,9,13,55,26,3,46,37,30,6,29,5,21,11,7,21,7,13,9,3,24,5,36,23,16,25,4,6,11,31,12,43,18,4],[55,5,3,4,9,4,9,15,17,6,4,8,10,14,4,4,7,8,19,5,4,12,5,5,2,7,5,4,12,17,6,6,4,1,3,5,16,5,2,19,10,4,7,6,2,5,8,7,9,14,24,3,7,3,3,9,7,5,7,11,9,22,11,2,8,16,4,10,14,7,13,10,2,10,11,4,30,5,8,3,20,3,6,7,5,9,5,8,10,4,12,9,10,0,9,2,7,18,2,6,11,26,11,15,23,9,4,10,6,9,6,9,4,9,4,1,3,8,7,13,3,8,20,10,5,5,6,12,3,1,5,9,3,8,11,10,17,1,4,5,18,4,12,6,25,13,11,0,10,19,22,5,10,4,19,17,2,8,6,18,7,9,7,7,4,8,3,2,4,23,11,1,24,13,14,4,16,7,12,5,2,9,3,11,7,4,20,6,13,6,14,14,2,4,9,17,4,32,3,4]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>PeakIDs<\/th>\n      <th>AP1_1<\/th>\n      <th>AP1_2A<\/th>\n      <th>AP1_2B<\/th>\n      <th>C_1A<\/th>\n      <th>C_1B<\/th>\n      <th>C_2A<\/th>\n      <th>C_2B<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7,8]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

### Differential binding analysis

The `runDiff` function performs differential binding analysis in batch mode for
several count tables using `edgeR` or `DESeq2` (Robinson, McCarthy, and Smyth 2010; Love, Huber, and Anders 2014).
Internally, it calls the functions `run_edgeR` and `run_DESeq2`. It also returns
the filtering results and plots from the downstream `filterDEGs` function using
the fold change and FDR cutoffs provided under the `dbrfilter` argument.

``` r
appendStep(sal) <- LineWise(code = {
    countDF_files <- getColumn(sal, step = "count_peak_ranges",
        "outfiles")
    outfiles <- paste0("./results/", names(countDF_files), "_peaks_edgeR.xls")
    names(outfiles) <- names(countDF_files)
    cmp <- readComp(file = stepsWF(sal)[["bowtie2_alignment"]],
        format = "matrix")
    dbrlist <- runDiff(args = countDF_files, outfiles = outfiles,
        diffFct = run_edgeR, targets = targetsWF(sal)[["bowtie2_alignment"]],
        cmp = cmp[[1]], independent = TRUE, dbrfilter = c(Fold = 2,
            FDR = 1))
}, step_name = "diff_bind_analysis", dependency = "count_peak_ranges",
    )
```

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

#### GO term enrichment analysis

The following performs GO term enrichment analysis for each annotated peak
set. Note: the following assumes that the GO annotation data exists under
`data/GO/catdb.RData`.

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
}, step_name = "get_go_annot", dependency = "annotation_ChIPseeker",
    )
```

#### GO term enrichment analysis

The following performs GO term enrichment analysis for each annotated peak
set. Note: the following assumes that the GO annotation data exists under
`data/GO/catdb.RData`.

``` r
appendStep(sal) <- LineWise(code = {
    annofiles <- getColumn(sal, step = "annotation_ChIPseeker",
        "outfiles")
    gene_ids <- sapply(annofiles, function(x) unique(as.character(read.delim(x)[,
        "geneId"])), simplify = FALSE)
    load("data/GO/catdb.RData")
    BatchResult <- GOCluster_Report(catdb = catdb, setlist = gene_ids,
        method = "all", id_type = "gene", CLSZ = 2, cutoff = 0.9,
        gocats = c("MF", "BP", "CC"), recordSpecGO = NULL)
    write.table(BatchResult, "results/GOBatchAll.xls", quote = FALSE,
        row.names = FALSE, sep = "\t")
}, step_name = "go_enrich", dependency = "annotation_ChIPseeker",
    )
```

Shows GO term enrichment results from previous step. The last gene identifier column (10)
of this table has been excluded in this viewing instance to minimize the complexity of the
result.
To avoid slowdowns of the load time of this page, only 50 rows of the source
table are imported into the below `datatable` view .

``` r
BatchResult <- read.delim("results/GOBatchAll.xls")[1:50, ]
library(DT)
datatable(BatchResult[, -10], options = list(scrollX = TRUE,
    autoWidth = TRUE))
```

<div id="htmlwidget-2" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-2">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50"],["C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A"],[4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903],["GO:0003700","GO:0140110","GO:0098772","GO:0043565","GO:0000976","GO:0001067","GO:1990837","GO:0003677","GO:0003690","GO:0042802","GO:0000981","GO:0046983","GO:0000987","GO:0003676","GO:0000978","GO:1901363","GO:0033612","GO:0097159","GO:0042562","GO:0005515","GO:0000977","GO:0060089","GO:0004672","GO:0019840","GO:0016773","GO:0043178","GO:0004860","GO:0019210","GO:0016705","GO:0015291","GO:0038023","GO:0001653","GO:0016301","GO:0022857","GO:0004861","GO:0030291","GO:0022804","GO:0042803","GO:0000156","GO:0005215","GO:0004674","GO:0019900","GO:0015293","GO:0005102","GO:0004805","GO:0005516","GO:0010427","GO:0080161","GO:0008509","GO:0019207"],[1667,1837,2497,1298,1007,1007,1032,2566,1142,394,336,664,306,4598,284,7859,46,7889,32,8016,324,236,1133,23,1296,39,27,30,400,313,207,17,1498,1200,8,8,515,208,19,1273,880,138,101,126,12,235,19,19,240,111],[682,698,807,506,412,412,416,777,432,133,117,195,102,966,94,1550,26,1552,20,1566,95,73,260,14,289,19,15,16,104,85,61,11,325,266,7,7,127,60,11,276,198,42,33,39,8,62,10,10,63,34],[3.49341724975865e-120,1.67823421235582e-104,2.18308574281622e-80,2.03153636203605e-78,2.16539339765293e-70,2.16539339765293e-70,7.55925494113875e-69,8.11545366519395e-63,5.44122511323004e-62,8.45596025352572e-15,3.02911820739919e-14,5.20591490791628e-14,2.60120228948827e-11,4.80717206533446e-11,2.37684179178145e-10,3.18096489951193e-09,3.70334206206444e-09,6.71010417736292e-09,2.23753662482306e-08,4.36325237409659e-08,1.70493636893304e-07,4.830138137427e-07,2.6496064818211e-06,4.80082717610204e-06,8.57859203202179e-06,9.09042267661985e-06,1.04457669068979e-05,1.06626043165579e-05,1.95557875198908e-05,2.00638856116351e-05,2.16501006655692e-05,2.30169720521419e-05,2.58444261599642e-05,3.08098675362986e-05,3.69557186521355e-05,3.69557186521355e-05,3.83718364858986e-05,5.00968427978002e-05,9.90476514846966e-05,0.000111473574255863,0.000134651553167959,0.000179143121967622,0.000200017440650394,0.000203270191184441,0.000242118052007588,0.000569854299398944,0.000584020859104921,0.000584020859104921,0.000596646365491369,0.000620806400430897],[2.12399768785326e-117,1.02036640111234e-101,1.32731613163226e-77,1.23517410811792e-75,1.31655918577298e-67,1.31655918577298e-67,4.59602700421236e-66,4.93419582843792e-60,3.30826486884386e-59,5.14122383414364e-12,1.8417038700987e-11,3.1651962640131e-11,1.58153099200887e-08,2.92276061572335e-08,1.44511980940312e-07,1.93402665890325e-06,2.25163197373518e-06,4.07974333983666e-06,1.36042226789242e-05,2.65285744345073e-05,0.000103660131231129,0.000293672398755561,0.00161096074094723,0.00291890292307004,0.00521578395546925,0.00552697698738487,0.00635102627939394,0.00648286342446718,0.0118899188120936,0.0121988424518741,0.0131632612046661,0.0139943190077023,0.0157134111052582,0.0187323994620695,0.0224690769404984,0.0224690769404984,0.0233300765834263,0.0304588804210625,0.0602209721026955,0.0677759331475647,0.0818681443261192,0.108919018156314,0.121610603915439,0.12358827624014,0.147207775620613,0.346471414034558,0.355084682335792,0.355084682335792,0.362760990218752,0.377450291461985],["DNA-binding transcription factor activity","transcription regulator activity","molecular function regulator","sequence-specific DNA binding","transcription regulatory region sequence-specific DNA binding","regulatory region nucleic acid binding","sequence-specific double-stranded DNA binding","DNA binding","double-stranded DNA binding","identical protein binding","DNA-binding transcription factor activity, RNA polymerase II-specific","protein dimerization activity","cis-regulatory region sequence-specific DNA binding","nucleic acid binding","RNA polymerase II cis-regulatory region sequence-specific DNA binding","heterocyclic compound binding","receptor serine/threonine kinase binding","organic cyclic compound binding","hormone binding","protein binding","RNA polymerase II transcription regulatory region sequence-specific DNA binding","molecular transducer activity","protein kinase activity","isoprenoid binding","phosphotransferase activity, alcohol group as acceptor","alcohol binding","protein kinase inhibitor activity","kinase inhibitor activity","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen","secondary active transmembrane transporter activity","signaling receptor activity","peptide receptor activity","kinase activity","transmembrane transporter activity","cyclin-dependent protein serine/threonine kinase inhibitor activity","protein serine/threonine kinase inhibitor activity","active transmembrane transporter activity","protein homodimerization activity","phosphorelay response regulator activity","transporter activity","protein serine/threonine kinase activity","kinase binding","symporter activity","signaling receptor binding","trehalose-phosphatase activity","calmodulin binding","abscisic acid binding","auxin transmembrane transporter activity","anion transmembrane transporter activity","kinase regulator activity"],["MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>CLID<\/th>\n      <th>CLSZ<\/th>\n      <th>GOID<\/th>\n      <th>NodeSize<\/th>\n      <th>SampleMatch<\/th>\n      <th>Phyper<\/th>\n      <th>Padj<\/th>\n      <th>Term<\/th>\n      <th>Ont<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"scrollX":true,"autoWidth":true,"columnDefs":[{"className":"dt-right","targets":[2,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

### Motif analysis

#### Parse DNA sequences of peak regions from genome

Enrichment analysis of known DNA binding motifs or *de novo* discovery
of novel motifs requires the DNA sequences of the identified peak
regions. To parse the corresponding sequences from the reference genome,
the `getSeq` function from the `Biostrings` package can be used. The
following example parses the sequences for each peak set and saves the
results to separate FASTA files, one for each peak set. In addition, the
sequences in the FASTA files are ranked (sorted) by increasing p-values
as expected by some motif discovery tools, such as `BCRANK`.

``` r
appendStep(sal) <- LineWise(code = {
    library(Biostrings)
    library(seqLogo)
    library(BCRANK)
    rangefiles <- getColumn(sal, step = "call_peaks_macs_withref",
        "outfiles")
    for (i in seq(along = rangefiles)) {
        df <- read.delim(rangefiles[i], comment = "#")
        peaks <- as(df, "GRanges")
        names(peaks) <- paste0(as.character(seqnames(peaks)),
            "_", start(peaks), "-", end(peaks))
        peaks <- peaks[order(values(peaks)$X.log10.pvalue., decreasing = TRUE)]
        pseq <- getSeq(FaFile("./data/tair10.fasta"), peaks)
        names(pseq) <- names(peaks)
        writeXStringSet(pseq, paste0(rangefiles[i], ".fasta"))
    }
}, step_name = "parse_peak_sequences", dependency = "call_peaks_macs_withref")
```

#### Motif discovery with `BCRANK`

The Bioconductor package `BCRANK` is one of the many tools available for
*de novo* discovery of DNA binding motifs in peak regions of ChIP-Seq
experiments. The given example applies this method on the first peak
sample set and plots the sequence logo of the highest ranking motif.

``` r
appendStep(sal) <- LineWise(code = {
    library(Biostrings)
    library(seqLogo)
    library(BCRANK)
    rangefiles <- getColumn(sal, step = "call_peaks_macs_withref",
        "outfiles")
    set.seed(0)
    BCRANKout <- bcrank(paste0(rangefiles[1], ".fasta"), restarts = 25,
        use.P1 = TRUE, use.P2 = TRUE)
    toptable(BCRANKout)
    topMotif <- toptable(BCRANKout, 1)
    weightMatrix <- pwm(topMotif, normalize = FALSE)
    weightMatrixNormalized <- pwm(topMotif, normalize = TRUE)
    pdf("results/seqlogo.pdf")
    seqLogo(weightMatrixNormalized)
    dev.off()
}, step_name = "bcrank_enrich", dependency = "call_peaks_macs_withref")
```

![](../results/seqlogo.png)

<div data-align="center">

Figure 2: One of the motifs identified by `BCRANK`

</div>

</br>

### Version Information

``` r
appendStep(sal) <- LineWise(code = {
    sessionInfo()
}, step_name = "sessionInfo", dependency = "bcrank_enrich")
```

## Running workflow

### Interactive job submissions in a single machine

For running the workflow, `runWF` function will execute all the steps store in
the workflow container. The execution will be on a single machine without
submitting to a queuing system of a computer cluster.

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
                  Njobs=8, 
                  walltime=120, ## minutes
                  ntasks=1,
                  ncpus=4, 
                  memory=1024, ## Mb
                  partition = "short"
                  )
sal <- addResources(sal, c("bowtie2_alignment"), resources = resources)
sal <- runWF(sal)
```

### Visualize workflow

*`systemPipeR`* workflows instances can be visualized with the `plotWF` function.

``` r
plotWF(sal)
```

<div id="htmlwidget-3" style="width:672px;height:480px;" class="plotwf html-widget"></div>
<script type="application/json" data-for="htmlwidget-3">{"x":{"dot":"digraph {\n    node[fontsize=20];\n    subgraph {\n        load_SPR -> preprocessing -> bowtie2_index -> bowtie2_alignment -> merge_bams -> writeTargetsRef -> call_peaks_macs_withref -> bcrank_enrich -> sessionInfo\n   }\n    preprocessing -> fastq_report\n    bowtie2_alignment -> align_stats\n    bowtie2_alignment -> bam_IGV\n    call_peaks_macs_withref -> annotation_ChIPseeker\n    call_peaks_macs_withref -> count_peak_ranges\n    call_peaks_macs_withref -> parse_peak_sequences\n    annotation_ChIPseeker -> ChIPseeker_plots\n    annotation_ChIPseeker -> go_enrich\n    count_peak_ranges -> diff_bind_analysis\n    load_SPR[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">load_SPR<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step load_SPR: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-30 18:12:57; End time: 2022-04-30 18:12:57; Duration: 00:00:00\"]\n    preprocessing[fillcolor=\"#d3d6eb\" style=\"filled, rounded\" label=<<b><font color=\"black\">preprocessing<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">7<\/font><\/b>; <font color=\"black\">0s<\/font>> , shape=\"box\"  tooltip=\"step preprocessing: 0 samples passed; 0 samples have warnings; 0 samples have errors; 7 samples in total; Start time: 2022-04-30 18:12:57; End time: 2022-04-30 18:12:57; Duration: 00:00:00\"]\n    fastq_report[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">fastq_report<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step fastq_report: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-30 18:12:57; End time: 2022-04-30 18:12:57; Duration: 00:00:00\"]\n    bowtie2_index[fillcolor=\"#d3d6eb\" style=\"filled, rounded\" label=<<b><font color=\"black\">bowtie2_index<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">6<\/font><\/b>; <font color=\"black\">0s<\/font>> , shape=\"box\"  tooltip=\"step bowtie2_index: 0 samples passed; 0 samples have warnings; 0 samples have errors; 6 samples in total; Start time: 2022-04-30 18:12:57; End time: 2022-04-30 18:12:57; Duration: 00:00:00\"]\n    bowtie2_alignment[fillcolor=\"#d3d6eb\" style=\"filled, rounded\" label=<<b><font color=\"black\">bowtie2_alignment<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">28<\/font><\/b>; <font color=\"black\">0s<\/font>> , shape=\"box\"  tooltip=\"step bowtie2_alignment: 0 samples passed; 0 samples have warnings; 0 samples have errors; 28 samples in total; Start time: 2022-04-30 18:12:57; End time: 2022-04-30 18:12:57; Duration: 00:00:00\"]\n    align_stats[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">align_stats<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step align_stats: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-30 18:12:57; End time: 2022-04-30 18:12:57; Duration: 00:00:00\"]\n    bam_IGV[style=\"solid, \"label=<<b><font color=\"black\">bam_IGV<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step bam_IGV: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-30 18:12:57; End time: 2022-04-30 18:12:57; Duration: 00:00:00\"]\n    merge_bams[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">merge_bams<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step merge_bams: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-30 18:12:57; End time: 2022-04-30 18:12:57; Duration: 00:00:00\"]\n    writeTargetsRef[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">writeTargetsRef<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step writeTargetsRef: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-30 18:12:57; End time: 2022-04-30 18:12:57; Duration: 00:00:00\"]\n    call_peaks_macs_withref[fillcolor=\"#d3d6eb\" style=\"filled, rounded\" label=<<b><font color=\"black\">call_peaks_macs_withref<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">5<\/font><\/b>; <font color=\"black\">0s<\/font>> , shape=\"box\"  tooltip=\"step call_peaks_macs_withref: 0 samples passed; 0 samples have warnings; 0 samples have errors; 5 samples in total; Start time: 2022-04-30 18:12:57; End time: 2022-04-30 18:12:57; Duration: 00:00:00\"]\n    annotation_ChIPseeker[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">annotation_ChIPseeker<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step annotation_ChIPseeker: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-30 18:12:57; End time: 2022-04-30 18:12:57; Duration: 00:00:00\"]\n    ChIPseeker_plots[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">ChIPseeker_plots<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step ChIPseeker_plots: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-30 18:12:57; End time: 2022-04-30 18:12:57; Duration: 00:00:00\"]\n    count_peak_ranges[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">count_peak_ranges<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step count_peak_ranges: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-30 18:12:57; End time: 2022-04-30 18:12:57; Duration: 00:00:00\"]\n    diff_bind_analysis[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">diff_bind_analysis<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step diff_bind_analysis: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-30 18:12:57; End time: 2022-04-30 18:12:57; Duration: 00:00:00\"]\n    go_enrich[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">go_enrich<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step go_enrich: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-30 18:12:57; End time: 2022-04-30 18:12:57; Duration: 00:00:00\"]\n    parse_peak_sequences[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">parse_peak_sequences<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step parse_peak_sequences: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-30 18:12:57; End time: 2022-04-30 18:12:57; Duration: 00:00:00\"]\n    bcrank_enrich[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">bcrank_enrich<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step bcrank_enrich: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-30 18:12:57; End time: 2022-04-30 18:12:57; Duration: 00:00:00\"]\n    sessionInfo[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">sessionInfo<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step sessionInfo: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-30 18:12:57; End time: 2022-04-30 18:12:57; Duration: 00:00:00\"]\n        subgraph cluster_legend {\n        rankdir=TB;\n        color=\"#eeeeee\";\n        style=filled;\n        ranksep =1;\n        label=\"Legends\";\n        fontsize = 30;\n        node [style=filled, fontsize=10];\n        legend_img-> step_state[color=\"#eeeeee\"];\n\n        legend_img[shape=none, image=\"plotwf_legend-src.png\", label = \" \", height=1, width=3, style=\"\"];\n\n        step_state[style=\"filled\", shape=\"box\" color=white, label =<\n            <table>\n            <tr><td><b>Step Colors<\/b><\/td><\/tr>\n            <tr><td><font color=\"black\">Pending steps<\/font>; <font color=\"#5cb85c\">Successful steps<\/font>; <font color=\"#d9534f\">Failed steps<\/font><\/td><\/tr>\n            <tr><td><b>Targets Files / Code Chunk <\/b><\/td><\/tr><tr><td><font color=\"#5cb85c\">0 (pass) <\/font> | <font color=\"#f0ad4e\">0 (warning) <\/font> | <font color=\"#d9534f\">0 (error) <\/font> | <font color=\"blue\">0 (total)<\/font>; Duration<\/td><\/tr><\/table>\n            >];\n\n    }\n\n}\n","plotid":"sprwf-75681234","responsive":true,"width":null,"height":null,"plot_method":"renderSVGElement","rmd":true,"msg":"","plot_ctr":true,"pan_zoom":false,"legend_uri":"data:image/svg+xml;base64,PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iVVRGLTgiPz4KPCEtLSBEbyBub3QgZWRpdCB0aGlzIGZpbGUgd2l0aCBlZGl0b3JzIG90aGVyIHRoYW4gZGlhZ3JhbXMubmV0IC0tPgo8IURPQ1RZUEUgc3ZnIFBVQkxJQyAiLS8vVzNDLy9EVEQgU1ZHIDEuMS8vRU4iICJodHRwOi8vd3d3LnczLm9yZy9HcmFwaGljcy9TVkcvMS4xL0RURC9zdmcxMS5kdGQiPgo8c3ZnIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyIgeG1sbnM6eGxpbms9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGxpbmsiIHZlcnNpb249IjEuMSIgd2lkdGg9IjQ5NnB4IiBoZWlnaHQ9IjI3OHB4IiB2aWV3Qm94PSItMC41IC0wLjUgNDk2IDI3OCIgY29udGVudD0iJmx0O214ZmlsZSBob3N0PSZxdW90O2FwcC5kaWFncmFtcy5uZXQmcXVvdDsgbW9kaWZpZWQ9JnF1b3Q7MjAyMS0xMS0yNFQyMDozOTo0NC45MjNaJnF1b3Q7IGFnZW50PSZxdW90OzUuMCAoWDExOyBMaW51eCB4ODZfNjQpIEFwcGxlV2ViS2l0LzUzNy4zNiAoS0hUTUwsIGxpa2UgR2Vja28pIENocm9tZS85Ni4wLjQ2NjQuNDUgU2FmYXJpLzUzNy4zNiZxdW90OyB2ZXJzaW9uPSZxdW90OzE1LjguMyZxdW90OyBldGFnPSZxdW90O1RfZEV4bkw3U0xYMmJ1WDhQYVgzJnF1b3Q7IHR5cGU9JnF1b3Q7Z29vZ2xlJnF1b3Q7Jmd0OyZsdDtkaWFncmFtIGlkPSZxdW90O1Z3OVZWaUdzeC1zTF9OaktlN0pQJnF1b3Q7Jmd0OzdWcmJjdUk0RVAwYUhwT3lKZDk0REFSbVhtWjNLMnpWUGl0WUdGVmt5V1BMZ2N6WHIyVEwrQ0lEM3NHUVRSVWtsZGd0cVMyZDAycDF0NW5BZWJ6L2xxSmsrNE9IbUU2QUZlNG44SGtDZ0EwQW1LaGZLL3pRRXN2eFMwbVVrbERMYXNHSy9NSlZSeTNOU1lpelZrZkJPUlVrYVF2WG5ERzhGaTBaU2xPK2EzZmJjTnArYW9JaWJBaFdhMFJONlQ4a0ZOdFNPcldzV3Y0ZGsyaGJQZG1yV21KVWRkYUNiSXRDdm11STRHSUM1eW5ub3J5SzkzTk1GWG9WTHVXNDVaSFd3OFJTek1TUUFacUpkMFJ6dlRZOUwvRlJMVFpLZVo1TTRFeitZeUZXNHl4NWQ1aTR1dEZLY0Nyd3ZvOEE5Rm9wczh3SjJvZGxTNFBCUE1ZaS9aQmRLa1dCSHFKdEJVQjl2NnVCQjQ1YnlyWk4wS0duQ2Rka1J3ZmROUjd5UWtQU0R3ODhEMDhibHQyV0NMeEswRnExN3FUNVM5bFd4Rkwvc3kwdk40VFNPYWM4TGNiQzBNVkI2Q2cwUmNyZmNOWENPSlBEWnhGRldkWUhkL2FHeFhwN0R2c214cUFmWTQycEE0ZEI2b3lBcUhOZFJKZkxSZkM4T0lib1ZWRU1wZ05SOUM5SDBiMHVpb3ZsRWl6bnQwWHhMR3p3Y3RpOHNXSGpUT2pUS1RCUjNManFSOG9SSlJHVE1vbzM0aGlvU2xWamJQbVJjaTRmVG9ReXNXQkVoSDNyU2c3VE54RE9PSldOWFpqbDlFVWJ5d3FrdFZ3R2xpak0xQ0tKUEhPZmRFTk13bEFObjZVNEk3LzBrYUl3U1RoaG9waXpPNXU0ejBwWExuaFc4bUlmQmJ4QlZvT0RKcDN5ZG9saVFoWDZmNU5ZUmh6QStnUHY1TjhYSGlNMmxJN1RCeHp3L1VlM2ZjUUIwNVBBSHNLcXcvUVN2Z0tEcnhCbFczd243RGhoMEhVTndnS0RNT2RLaEUwTnduNGdKdVBOV00zK3E1Rm0yN2ZiWm03dzZIZG9jODFRc3M4eGprRmJsYmswZUp2ek9Na0Z2cE4yWXFzNXhsYnI0Y3k3Rm1kMlQ3amcwZUlJVCtRU202eDVQM09Wc2hXNFBKVDRQc2tPUWJLdjIrUlZwUDdMcDhoRUZWZTY1RFJLZFdXcllROVY5LzhVeElWNGczSXFicmU5dktETGxPOFpUTG4yWTA5a0IwWUlpRzB6a1gyNTc2d1RmQVYyaHk5bzlXU0IxOXBaWmw0dHZhRmNYUGhBQ2J1N3hGTXVFWFpkSXJSc2c3amdXc1NaNlh1cTZqc3NvaVpyZFNwbFgrSzhQaUhBQTBOOEdlajFaZTRJSUEvSTdqRUxuMVFoc3piSEJwYWRTbDFoejFXZEVoNHd3cUZSNUR5TFVIUDlQU1pXeVZKTWtTRHZiZlY5a09nbi9LWDJaUE13bVQ2Q0RnUGRxa25HODNTTjlVRFFLSEIyZEVGZ0czRmZWNWRBYVlTRm9hdmc2YkQ0WWRRTnFEQmNSSjNPenNyTy8zY2VaYWJVQWQ1N0RLeHAvWEYvajFQSE1sTG1NNXBIWk5pc2NNZ0VMRVNDNjAxK1A3ZU9CSWhkcHdvZDA2bGVMWlEzNnh4L0pvSndodWlkdEZNbm9UV0F0RDRQTWdwcFpxMWpsUHhMd1RzMCticm83UTRNUGZ4NnF1QjdpK3BpKy9VWkJHYVU3L2k5a2N3SURBS3o2akVLZzR6TGxsdlJ1Q3crbjB1akdaTGVtRWl6RlBLU00wWllORkhxU3hKZTB5R0UycFprdEtDeVErb0taNW4weUYvSUgrczVkdDRFMmRiTnpNSU9nbllZQkwzS0FzNzRaemlHVVpnMWw4b1VRdkora1MyYzlSS3ZhUDBXRmJ2NllWMHlvL1NKRkxHc2dtK21mTVRBNXgyTStZZzdPWWlMaGVreCtHZE8wcTlXZVA5OG0vWGJOZ3VCR1ZIMFZlRkhzVml6N3JRU09CblBoODJMcjJ6Y3JXRzROVXhoMjRNRlpqR3I3MVhhYjFpRHZLMi9WMVZtZnZYWDArRGlYdz09Jmx0Oy9kaWFncmFtJmd0OyZsdDsvbXhmaWxlJmd0OyI+PGRlZnMvPjxnPjxyZWN0IHg9IjQiIHk9IjkwIiB3aWR0aD0iNDkwIiBoZWlnaHQ9IjkyIiBmaWxsPSIjZDVlOGQ0IiBzdHJva2U9Im5vbmUiIHBvaW50ZXItZXZlbnRzPSJub25lIi8+PHJlY3QgeD0iNCIgeT0iMTgyIiB3aWR0aD0iNDkwIiBoZWlnaHQ9Ijk0IiBmaWxsPSIjZmZlOGRlIiBzdHJva2U9Im5vbmUiIHBvaW50ZXItZXZlbnRzPSJub25lIi8+PHJlY3QgeD0iNCIgeT0iNCIgd2lkdGg9IjQ5MCIgaGVpZ2h0PSI4NiIgZmlsbD0iI2VmZjJmYyIgc3Ryb2tlPSJub25lIiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxyZWN0IHg9IjQiIHk9IjQiIHdpZHRoPSIxNDAiIGhlaWdodD0iMjcyIiBmaWxsLW9wYWNpdHk9IjAuOCIgZmlsbD0iI2Y1ZjVmNSIgc3Ryb2tlPSJub25lIiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiAxMXB4OyBtYXJnaW4tbGVmdDogMTE1cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogOHB4OyBmb250LWZhbWlseTogJnF1b3Q7VGltZXMgTmV3IFJvbWFuJnF1b3Q7OyBjb2xvcjogcmdiKDAsIDAsIDApOyBsaW5lLWhlaWdodDogMS4yOyBwb2ludGVyLWV2ZW50czogbm9uZTsgd2hpdGUtc3BhY2U6IG5vd3JhcDsiPnNvbGlkPC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjExNSIgeT0iMTMiIGZpbGw9InJnYmEoMCwgMCwgMCwgMSkiIGZvbnQtZmFtaWx5PSJUaW1lcyBOZXcgUm9tYW4iIGZvbnQtc2l6ZT0iOHB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5zb2xpZDwvdGV4dD48L3N3aXRjaD48L2c+PGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTAuNSAtMC41KXNjYWxlKDIpIj48c3dpdGNoPjxmb3JlaWduT2JqZWN0IHBvaW50ZXItZXZlbnRzPSJub25lIiB3aWR0aD0iMTAwJSIgaGVpZ2h0PSIxMDAlIiByZXF1aXJlZEZlYXR1cmVzPSJodHRwOi8vd3d3LnczLm9yZy9UUi9TVkcxMS9mZWF0dXJlI0V4dGVuc2liaWxpdHkiIHN0eWxlPSJvdmVyZmxvdzogdmlzaWJsZTsgdGV4dC1hbGlnbjogbGVmdDsiPjxkaXYgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGh0bWwiIHN0eWxlPSJkaXNwbGF5OiBmbGV4OyBhbGlnbi1pdGVtczogdW5zYWZlIGNlbnRlcjsganVzdGlmeS1jb250ZW50OiB1bnNhZmUgY2VudGVyOyB3aWR0aDogMXB4OyBoZWlnaHQ6IDFweDsgcGFkZGluZy10b3A6IDEwcHg7IG1hcmdpbi1sZWZ0OiAxOThweDsiPjxkaXYgZGF0YS1kcmF3aW8tY29sb3JzPSJjb2xvcjogcmdiYSgwLCAwLCAwLCAxKTsgIiBzdHlsZT0iYm94LXNpemluZzogYm9yZGVyLWJveDsgZm9udC1zaXplOiAwcHg7IHRleHQtYWxpZ246IGNlbnRlcjsiPjxkaXYgc3R5bGU9ImRpc3BsYXk6IGlubGluZS1ibG9jazsgZm9udC1zaXplOiA4cHg7IGZvbnQtZmFtaWx5OiAmcXVvdDtUaW1lcyBOZXcgUm9tYW4mcXVvdDs7IGNvbG9yOiByZ2IoMCwgMCwgMCk7IGxpbmUtaGVpZ2h0OiAxLjI7IHBvaW50ZXItZXZlbnRzOiBub25lOyB3aGl0ZS1zcGFjZTogbm93cmFwOyI+ZGFzaGVkPC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjE5OCIgeT0iMTIiIGZpbGw9InJnYmEoMCwgMCwgMCwgMSkiIGZvbnQtZmFtaWx5PSJUaW1lcyBOZXcgUm9tYW4iIGZvbnQtc2l6ZT0iOHB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5kYXNoZWQ8L3RleHQ+PC9zd2l0Y2g+PC9nPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiAzMnB4OyBtYXJnaW4tbGVmdDogMTE2cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTFweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3dyYXA7Ij5NYW5hZ2VtZW50PC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjExNiIgeT0iMzUiIGZpbGw9InJnYmEoMCwgMCwgMCwgMSkiIGZvbnQtZmFtaWx5PSJUaW1lcyBOZXcgUm9tYW4iIGZvbnQtc2l6ZT0iMTFweCIgdGV4dC1hbmNob3I9Im1pZGRsZSI+TWFuYWdlbWVudDwvdGV4dD48L3N3aXRjaD48L2c+PGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTAuNSAtMC41KXNjYWxlKDIpIj48c3dpdGNoPjxmb3JlaWduT2JqZWN0IHBvaW50ZXItZXZlbnRzPSJub25lIiB3aWR0aD0iMTAwJSIgaGVpZ2h0PSIxMDAlIiByZXF1aXJlZEZlYXR1cmVzPSJodHRwOi8vd3d3LnczLm9yZy9UUi9TVkcxMS9mZWF0dXJlI0V4dGVuc2liaWxpdHkiIHN0eWxlPSJvdmVyZmxvdzogdmlzaWJsZTsgdGV4dC1hbGlnbjogbGVmdDsiPjxkaXYgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGh0bWwiIHN0eWxlPSJkaXNwbGF5OiBmbGV4OyBhbGlnbi1pdGVtczogdW5zYWZlIGNlbnRlcjsganVzdGlmeS1jb250ZW50OiB1bnNhZmUgY2VudGVyOyB3aWR0aDogMXB4OyBoZWlnaHQ6IDFweDsgcGFkZGluZy10b3A6IDMycHg7IG1hcmdpbi1sZWZ0OiAxOThweDsiPjxkaXYgZGF0YS1kcmF3aW8tY29sb3JzPSJjb2xvcjogcmdiYSgwLCAwLCAwLCAxKTsgIiBzdHlsZT0iYm94LXNpemluZzogYm9yZGVyLWJveDsgZm9udC1zaXplOiAwcHg7IHRleHQtYWxpZ246IGNlbnRlcjsiPjxkaXYgc3R5bGU9ImRpc3BsYXk6IGlubGluZS1ibG9jazsgZm9udC1zaXplOiAxMXB4OyBmb250LWZhbWlseTogJnF1b3Q7VGltZXMgTmV3IFJvbWFuJnF1b3Q7OyBjb2xvcjogcmdiKDAsIDAsIDApOyBsaW5lLWhlaWdodDogMS4yOyBwb2ludGVyLWV2ZW50czogbm9uZTsgd2hpdGUtc3BhY2U6IG5vd3JhcDsiPkNvbXB1dGU8L2Rpdj48L2Rpdj48L2Rpdj48L2ZvcmVpZ25PYmplY3Q+PHRleHQgeD0iMTk4IiB5PSIzNSIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMXB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5Db21wdXRlPC90ZXh0Pjwvc3dpdGNoPjwvZz48ZWxsaXBzZSBjeD0iMjMyLjUiIGN5PSIxMjMiIHJ4PSI1MS41IiByeT0iMjciIGZpbGw9InJnYmEoMjU1LCAyNTUsIDI1NSwgMSkiIHN0cm9rZT0icmdiYSgwLCAwLCAwLCAxKSIgc3Ryb2tlLXdpZHRoPSIyIiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDUwcHg7IGhlaWdodDogMXB4OyBwYWRkaW5nLXRvcDogNjJweDsgbWFyZ2luLWxlZnQ6IDkycHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTJweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3JtYWw7IG92ZXJmbG93LXdyYXA6IG5vcm1hbDsiPjxzcGFuIHN0eWxlPSJmb250LXNpemU6IDhweCI+ZWxsaXBzZTwvc3Bhbj48L2Rpdj48L2Rpdj48L2Rpdj48L2ZvcmVpZ25PYmplY3Q+PHRleHQgeD0iMTE2IiB5PSI2NSIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMnB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5lbGxpcHNlPC90ZXh0Pjwvc3dpdGNoPjwvZz48ZyB0cmFuc2Zvcm09InRyYW5zbGF0ZSgtMC41IC0wLjUpc2NhbGUoMikiPjxzd2l0Y2g+PGZvcmVpZ25PYmplY3QgcG9pbnRlci1ldmVudHM9Im5vbmUiIHdpZHRoPSIxMDAlIiBoZWlnaHQ9IjEwMCUiIHJlcXVpcmVkRmVhdHVyZXM9Imh0dHA6Ly93d3cudzMub3JnL1RSL1NWRzExL2ZlYXR1cmUjRXh0ZW5zaWJpbGl0eSIgc3R5bGU9Im92ZXJmbG93OiB2aXNpYmxlOyB0ZXh0LWFsaWduOiBsZWZ0OyI+PGRpdiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMTk5OS94aHRtbCIgc3R5bGU9ImRpc3BsYXk6IGZsZXg7IGFsaWduLWl0ZW1zOiB1bnNhZmUgY2VudGVyOyBqdXN0aWZ5LWNvbnRlbnQ6IHVuc2FmZSBjZW50ZXI7IHdpZHRoOiAxcHg7IGhlaWdodDogMXB4OyBwYWRkaW5nLXRvcDogODVweDsgbWFyZ2luLWxlZnQ6IDExNHB4OyI+PGRpdiBkYXRhLWRyYXdpby1jb2xvcnM9ImNvbG9yOiByZ2JhKDAsIDAsIDAsIDEpOyAiIHN0eWxlPSJib3gtc2l6aW5nOiBib3JkZXItYm94OyBmb250LXNpemU6IDBweDsgdGV4dC1hbGlnbjogY2VudGVyOyI+PGRpdiBzdHlsZT0iZGlzcGxheTogaW5saW5lLWJsb2NrOyBmb250LXNpemU6IDExcHg7IGZvbnQtZmFtaWx5OiAmcXVvdDtUaW1lcyBOZXcgUm9tYW4mcXVvdDs7IGNvbG9yOiByZ2IoMCwgMCwgMCk7IGxpbmUtaGVpZ2h0OiAxLjI7IHBvaW50ZXItZXZlbnRzOiBub25lOyB3aGl0ZS1zcGFjZTogbm93cmFwOyI+UjwvZGl2PjwvZGl2PjwvZGl2PjwvZm9yZWlnbk9iamVjdD48dGV4dCB4PSIxMTQiIHk9Ijg4IiBmaWxsPSJyZ2JhKDAsIDAsIDAsIDEpIiBmb250LWZhbWlseT0iVGltZXMgTmV3IFJvbWFuIiBmb250LXNpemU9IjExcHgiIHRleHQtYW5jaG9yPSJtaWRkbGUiPlI8L3RleHQ+PC9zd2l0Y2g+PC9nPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiA4M3B4OyBtYXJnaW4tbGVmdDogMTk4cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTFweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3dyYXA7Ij5Db21tYW5kLWxpbmU8L2Rpdj48L2Rpdj48L2Rpdj48L2ZvcmVpZ25PYmplY3Q+PHRleHQgeD0iMTk4IiB5PSI4NiIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMXB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5Db21tYW5kLWxpbmU8L3RleHQ+PC9zd2l0Y2g+PC9nPjxyZWN0IHg9IjM0OSIgeT0iOTYiIHdpZHRoPSIxMDUiIGhlaWdodD0iNTAiIHJ4PSI3LjUiIHJ5PSI3LjUiIGZpbGw9InJnYmEoMjU1LCAyNTUsIDI1NSwgMSkiIHN0cm9rZT0icmdiYSgwLCAwLCAwLCAxKSIgc3Ryb2tlLXdpZHRoPSIyIiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDUxcHg7IGhlaWdodDogMXB4OyBwYWRkaW5nLXRvcDogNjFweDsgbWFyZ2luLWxlZnQ6IDE3NnB4OyI+PGRpdiBkYXRhLWRyYXdpby1jb2xvcnM9ImNvbG9yOiByZ2JhKDAsIDAsIDAsIDEpOyAiIHN0eWxlPSJib3gtc2l6aW5nOiBib3JkZXItYm94OyBmb250LXNpemU6IDBweDsgdGV4dC1hbGlnbjogY2VudGVyOyI+PGRpdiBzdHlsZT0iZGlzcGxheTogaW5saW5lLWJsb2NrOyBmb250LXNpemU6IDhweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3JtYWw7IG92ZXJmbG93LXdyYXA6IG5vcm1hbDsiPnJlY3RhbmdsZTwvZGl2PjwvZGl2PjwvZGl2PjwvZm9yZWlnbk9iamVjdD48dGV4dCB4PSIyMDEiIHk9IjYzIiBmaWxsPSJyZ2JhKDAsIDAsIDAsIDEpIiBmb250LWZhbWlseT0iVGltZXMgTmV3IFJvbWFuIiBmb250LXNpemU9IjhweCIgdGV4dC1hbmNob3I9Im1pZGRsZSI+cmVjdGFuZ2xlPC90ZXh0Pjwvc3dpdGNoPjwvZz48cGF0aCBkPSJNIDE4Mi41IDM4IEwgMjg3LjUgMzgiIGZpbGw9Im5vbmUiIHN0cm9rZT0icmdiYSgwLCAwLCAwLCAxKSIgc3Ryb2tlLXdpZHRoPSI2IiBzdHJva2UtbWl0ZXJsaW1pdD0iMTAiIHBvaW50ZXItZXZlbnRzPSJub25lIi8+PHBhdGggZD0iTSAzNTQgMzcuNjIgTCA0NTkgMzcuNjIiIGZpbGw9Im5vbmUiIHN0cm9rZT0icmdiYSgwLCAwLCAwLCAxKSIgc3Ryb2tlLXdpZHRoPSI2IiBzdHJva2UtbWl0ZXJsaW1pdD0iMTAiIHN0cm9rZS1kYXNoYXJyYXk9IjE4IDE4IiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiAxMjhweDsgbWFyZ2luLWxlZnQ6IDExNXB4OyI+PGRpdiBkYXRhLWRyYXdpby1jb2xvcnM9ImNvbG9yOiByZ2JhKDAsIDAsIDAsIDEpOyAiIHN0eWxlPSJib3gtc2l6aW5nOiBib3JkZXItYm94OyBmb250LXNpemU6IDBweDsgdGV4dC1hbGlnbjogY2VudGVyOyI+PGRpdiBzdHlsZT0iZGlzcGxheTogaW5saW5lLWJsb2NrOyBmb250LXNpemU6IDExcHg7IGZvbnQtZmFtaWx5OiAmcXVvdDtUaW1lcyBOZXcgUm9tYW4mcXVvdDs7IGNvbG9yOiByZ2IoMCwgMCwgMCk7IGxpbmUtaGVpZ2h0OiAxLjI7IHBvaW50ZXItZXZlbnRzOiBub25lOyB3aGl0ZS1zcGFjZTogbm93cmFwOyI+TWFuZGF0b3J5PC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjExNSIgeT0iMTMxIiBmaWxsPSJyZ2JhKDAsIDAsIDAsIDEpIiBmb250LWZhbWlseT0iVGltZXMgTmV3IFJvbWFuIiBmb250LXNpemU9IjExcHgiIHRleHQtYW5jaG9yPSJtaWRkbGUiPk1hbmRhdG9yeTwvdGV4dD48L3N3aXRjaD48L2c+PGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTAuNSAtMC41KXNjYWxlKDIpIj48c3dpdGNoPjxmb3JlaWduT2JqZWN0IHBvaW50ZXItZXZlbnRzPSJub25lIiB3aWR0aD0iMTAwJSIgaGVpZ2h0PSIxMDAlIiByZXF1aXJlZEZlYXR1cmVzPSJodHRwOi8vd3d3LnczLm9yZy9UUi9TVkcxMS9mZWF0dXJlI0V4dGVuc2liaWxpdHkiIHN0eWxlPSJvdmVyZmxvdzogdmlzaWJsZTsgdGV4dC1hbGlnbjogbGVmdDsiPjxkaXYgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGh0bWwiIHN0eWxlPSJkaXNwbGF5OiBmbGV4OyBhbGlnbi1pdGVtczogdW5zYWZlIGNlbnRlcjsganVzdGlmeS1jb250ZW50OiB1bnNhZmUgY2VudGVyOyB3aWR0aDogMXB4OyBoZWlnaHQ6IDFweDsgcGFkZGluZy10b3A6IDEyOHB4OyBtYXJnaW4tbGVmdDogMTk4cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTFweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3dyYXA7Ij5PcHRpb25hbDwvZGl2PjwvZGl2PjwvZGl2PjwvZm9yZWlnbk9iamVjdD48dGV4dCB4PSIxOTgiIHk9IjEzMSIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMXB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5PcHRpb25hbDwvdGV4dD48L3N3aXRjaD48L2c+PHJlY3QgeD0iMTg0IiB5PSIxOTAiIHdpZHRoPSI5NSIgaGVpZ2h0PSI0MCIgZmlsbD0iI2QzZDZlYiIgc3Ryb2tlPSJub25lIiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDQ2cHg7IGhlaWdodDogMXB4OyBwYWRkaW5nLXRvcDogMTA1cHg7IG1hcmdpbi1sZWZ0OiA5M3B4OyI+PGRpdiBkYXRhLWRyYXdpby1jb2xvcnM9ImNvbG9yOiByZ2JhKDAsIDAsIDAsIDEpOyAiIHN0eWxlPSJib3gtc2l6aW5nOiBib3JkZXItYm94OyBmb250LXNpemU6IDBweDsgdGV4dC1hbGlnbjogY2VudGVyOyI+PGRpdiBzdHlsZT0iZGlzcGxheTogaW5saW5lLWJsb2NrOyBmb250LXNpemU6IDEycHg7IGZvbnQtZmFtaWx5OiAmcXVvdDtUaW1lcyBOZXcgUm9tYW4mcXVvdDs7IGNvbG9yOiByZ2IoMCwgMCwgMCk7IGxpbmUtaGVpZ2h0OiAxLjI7IHBvaW50ZXItZXZlbnRzOiBub25lOyB3aGl0ZS1zcGFjZTogbm9ybWFsOyBvdmVyZmxvdy13cmFwOiBub3JtYWw7Ij48c3BhbiBzdHlsZT0iZm9udC1zaXplOiA4cHgiPmZpbGw8L3NwYW4+PC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjExNiIgeT0iMTA5IiBmaWxsPSJyZ2JhKDAsIDAsIDAsIDEpIiBmb250LWZhbWlseT0iVGltZXMgTmV3IFJvbWFuIiBmb250LXNpemU9IjEycHgiIHRleHQtYW5jaG9yPSJtaWRkbGUiPmZpbGw8L3RleHQ+PC9zd2l0Y2g+PC9nPjxyZWN0IHg9IjM0OSIgeT0iMTkwIiB3aWR0aD0iOTUiIGhlaWdodD0iNDAiIGZpbGw9IiNmZmZmZmYiIHN0cm9rZT0ibm9uZSIgcG9pbnRlci1ldmVudHM9Im5vbmUiLz48ZyB0cmFuc2Zvcm09InRyYW5zbGF0ZSgtMC41IC0wLjUpc2NhbGUoMikiPjxzd2l0Y2g+PGZvcmVpZ25PYmplY3QgcG9pbnRlci1ldmVudHM9Im5vbmUiIHdpZHRoPSIxMDAlIiBoZWlnaHQ9IjEwMCUiIHJlcXVpcmVkRmVhdHVyZXM9Imh0dHA6Ly93d3cudzMub3JnL1RSL1NWRzExL2ZlYXR1cmUjRXh0ZW5zaWJpbGl0eSIgc3R5bGU9Im92ZXJmbG93OiB2aXNpYmxlOyB0ZXh0LWFsaWduOiBsZWZ0OyI+PGRpdiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMTk5OS94aHRtbCIgc3R5bGU9ImRpc3BsYXk6IGZsZXg7IGFsaWduLWl0ZW1zOiB1bnNhZmUgY2VudGVyOyBqdXN0aWZ5LWNvbnRlbnQ6IHVuc2FmZSBjZW50ZXI7IHdpZHRoOiA0NnB4OyBoZWlnaHQ6IDFweDsgcGFkZGluZy10b3A6IDEwNXB4OyBtYXJnaW4tbGVmdDogMTc2cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTJweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3JtYWw7IG92ZXJmbG93LXdyYXA6IG5vcm1hbDsiPjxzcGFuIHN0eWxlPSJmb250LXNpemU6IDhweCI+bm8gZmlsbDwvc3Bhbj48L2Rpdj48L2Rpdj48L2Rpdj48L2ZvcmVpZ25PYmplY3Q+PHRleHQgeD0iMTk4IiB5PSIxMDkiIGZpbGw9InJnYmEoMCwgMCwgMCwgMSkiIGZvbnQtZmFtaWx5PSJUaW1lcyBOZXcgUm9tYW4iIGZvbnQtc2l6ZT0iMTJweCIgdGV4dC1hbmNob3I9Im1pZGRsZSI+bm8gZmlsbDwvdGV4dD48L3N3aXRjaD48L2c+PGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTAuNSAtMC41KXNjYWxlKDIpIj48c3dpdGNoPjxmb3JlaWduT2JqZWN0IHBvaW50ZXItZXZlbnRzPSJub25lIiB3aWR0aD0iMTAwJSIgaGVpZ2h0PSIxMDAlIiByZXF1aXJlZEZlYXR1cmVzPSJodHRwOi8vd3d3LnczLm9yZy9UUi9TVkcxMS9mZWF0dXJlI0V4dGVuc2liaWxpdHkiIHN0eWxlPSJvdmVyZmxvdzogdmlzaWJsZTsgdGV4dC1hbGlnbjogbGVmdDsiPjxkaXYgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGh0bWwiIHN0eWxlPSJkaXNwbGF5OiBmbGV4OyBhbGlnbi1pdGVtczogdW5zYWZlIGNlbnRlcjsganVzdGlmeS1jb250ZW50OiB1bnNhZmUgY2VudGVyOyB3aWR0aDogMXB4OyBoZWlnaHQ6IDFweDsgcGFkZGluZy10b3A6IDI0cHg7IG1hcmdpbi1sZWZ0OiAzNXB4OyI+PGRpdiBkYXRhLWRyYXdpby1jb2xvcnM9ImNvbG9yOiByZ2JhKDAsIDAsIDAsIDEpOyAiIHN0eWxlPSJib3gtc2l6aW5nOiBib3JkZXItYm94OyBmb250LXNpemU6IDBweDsgdGV4dC1hbGlnbjogY2VudGVyOyI+PGRpdiBzdHlsZT0iZGlzcGxheTogaW5saW5lLWJsb2NrOyBmb250LXNpemU6IDEwcHg7IGZvbnQtZmFtaWx5OiAmcXVvdDtUaW1lcyBOZXcgUm9tYW4mcXVvdDs7IGNvbG9yOiByZ2IoMCwgMCwgMCk7IGxpbmUtaGVpZ2h0OiAxLjI7IHBvaW50ZXItZXZlbnRzOiBub25lOyBmb250LXdlaWdodDogYm9sZDsgd2hpdGUtc3BhY2U6IG5vd3JhcDsiPlJ1bm5pbmcgPGJyIHN0eWxlPSJmb250LXNpemU6IDEwcHgiIC8+U2Vzc2lvbjwvZGl2PjwvZGl2PjwvZGl2PjwvZm9yZWlnbk9iamVjdD48dGV4dCB4PSIzNSIgeT0iMjciIGZpbGw9InJnYmEoMCwgMCwgMCwgMSkiIGZvbnQtZmFtaWx5PSJUaW1lcyBOZXcgUm9tYW4iIGZvbnQtc2l6ZT0iMTBweCIgdGV4dC1hbmNob3I9Im1pZGRsZSIgZm9udC13ZWlnaHQ9ImJvbGQiPlJ1bm5pbmcuLi48L3RleHQ+PC9zd2l0Y2g+PC9nPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiAxMTNweDsgbWFyZ2luLWxlZnQ6IDM1cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTBweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IGZvbnQtd2VpZ2h0OiBib2xkOyB3aGl0ZS1zcGFjZTogbm93cmFwOyI+PGRpdiBzdHlsZT0iZm9udC1zaXplOiAxMHB4Ij48c3BhbiBzdHlsZT0iYmFja2dyb3VuZC1jb2xvcjogdHJhbnNwYXJlbnQgOyBmb250LXNpemU6IDEwcHgiPlJ1bm5pbmc8L3NwYW4+PC9kaXY+UmVxdWlyZW1lbnQ8L2Rpdj48L2Rpdj48L2Rpdj48L2ZvcmVpZ25PYmplY3Q+PHRleHQgeD0iMzUiIHk9IjExNiIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMHB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIiBmb250LXdlaWdodD0iYm9sZCI+UnVubmluZ1JlcXVpcmUuLi48L3RleHQ+PC9zd2l0Y2g+PC9nPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiA2OHB4OyBtYXJnaW4tbGVmdDogMzVweDsiPjxkaXYgZGF0YS1kcmF3aW8tY29sb3JzPSJjb2xvcjogcmdiYSgwLCAwLCAwLCAxKTsgIiBzdHlsZT0iYm94LXNpemluZzogYm9yZGVyLWJveDsgZm9udC1zaXplOiAwcHg7IHRleHQtYWxpZ246IGNlbnRlcjsiPjxkaXYgc3R5bGU9ImRpc3BsYXk6IGlubGluZS1ibG9jazsgZm9udC1zaXplOiAxMHB4OyBmb250LWZhbWlseTogJnF1b3Q7VGltZXMgTmV3IFJvbWFuJnF1b3Q7OyBjb2xvcjogcmdiKDAsIDAsIDApOyBsaW5lLWhlaWdodDogMS4yOyBwb2ludGVyLWV2ZW50czogbm9uZTsgZm9udC13ZWlnaHQ6IGJvbGQ7IHdoaXRlLXNwYWNlOiBub3dyYXA7Ij5TdGVwIDxiciBzdHlsZT0iZm9udC1zaXplOiAxMHB4IiAvPkNsYXNzPC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjM1IiB5PSI3MSIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMHB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIiBmb250LXdlaWdodD0iYm9sZCI+U3RlcC4uLjwvdGV4dD48L3N3aXRjaD48L2c+PC9nPjxzd2l0Y2g+PGcgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5Ii8+PGEgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoMCwtNSkiIHhsaW5rOmhyZWY9Imh0dHBzOi8vd3d3LmRpYWdyYW1zLm5ldC9kb2MvZmFxL3N2Zy1leHBvcnQtdGV4dC1wcm9ibGVtcyIgdGFyZ2V0PSJfYmxhbmsiPjx0ZXh0IHRleHQtYW5jaG9yPSJtaWRkbGUiIGZvbnQtc2l6ZT0iMTBweCIgeD0iNTAlIiB5PSIxMDAlIj5WaWV3ZXIgZG9lcyBub3Qgc3VwcG9ydCBmdWxsIFNWRyAxLjE8L3RleHQ+PC9hPjwvc3dpdGNoPjwvc3ZnPg=="},"evals":[],"jsHooks":[]}</script>

### Checking workflow status

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

Alternatively, scientific reports can be rendered with the `render` function from `rmarkdown`.

``` r
rmarkdown::render("systemPipeChIPseq.Rmd", clean = TRUE, output_format = "BiocStyle::html_document")
```

## Funding

This project was supported by funds from the National Institutes of
Health (NIH) and the National Science Foundation (NSF).

## Session Info

``` r
sessionInfo()
```

    ## R version 4.3.0 (2023-04-21)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Debian GNU/Linux 11 (bullseye)
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils    
    ## [6] datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] DT_0.27                     systemPipeRdata_2.4.0      
    ##  [3] systemPipeR_2.6.0           ShortRead_1.58.0           
    ##  [5] GenomicAlignments_1.36.0    SummarizedExperiment_1.30.1
    ##  [7] Biobase_2.60.0              MatrixGenerics_1.12.0      
    ##  [9] matrixStats_0.63.0          BiocParallel_1.34.1        
    ## [11] Rsamtools_2.16.0            Biostrings_2.68.0          
    ## [13] XVector_0.40.0              GenomicRanges_1.52.0       
    ## [15] GenomeInfoDb_1.36.0         IRanges_2.34.0             
    ## [17] S4Vectors_0.38.1            BiocGenerics_0.46.0        
    ## [19] BiocStyle_2.28.0           
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.3            xfun_0.39              
    ##  [3] bslib_0.4.2             hwriter_1.3.2.1        
    ##  [5] ggplot2_3.4.2           remotes_2.4.2          
    ##  [7] htmlwidgets_1.6.2       latticeExtra_0.6-30    
    ##  [9] lattice_0.21-8          generics_0.1.3         
    ## [11] vctrs_0.6.2             tools_4.3.0            
    ## [13] bitops_1.0-7            parallel_4.3.0         
    ## [15] tibble_3.2.1            fansi_1.0.4            
    ## [17] pkgconfig_2.0.3         Matrix_1.5-4           
    ## [19] RColorBrewer_1.1-3      lifecycle_1.0.3        
    ## [21] GenomeInfoDbData_1.2.10 stringr_1.5.0          
    ## [23] compiler_4.3.0          deldir_1.0-6           
    ## [25] munsell_0.5.0           codetools_0.2-19       
    ## [27] htmltools_0.5.5         sass_0.4.6             
    ## [29] RCurl_1.98-1.12         yaml_2.3.7             
    ## [31] pillar_1.9.0            crayon_1.5.2           
    ## [33] jquerylib_0.1.4         DelayedArray_0.26.2    
    ## [35] cachem_1.0.8            tidyselect_1.2.0       
    ## [37] digest_0.6.31           stringi_1.7.12         
    ## [39] dplyr_1.1.2             bookdown_0.33          
    ## [41] fastmap_1.1.1           grid_4.3.0             
    ## [43] colorspace_2.1-0        cli_3.6.1              
    ## [45] magrittr_2.0.3          S4Arrays_1.0.1         
    ## [47] utf8_1.2.3              scales_1.2.1           
    ## [49] rmarkdown_2.21          jpeg_0.1-10            
    ## [51] interp_1.1-4            blogdown_1.16          
    ## [53] png_0.1-8               evaluate_0.21          
    ## [55] knitr_1.42              rlang_1.1.1            
    ## [57] Rcpp_1.0.10             glue_1.6.2             
    ## [59] formatR_1.14            BiocManager_1.30.20    
    ## [61] jsonlite_1.8.4          R6_2.5.1               
    ## [63] zlibbioc_1.46.0

## References

<div id="refs" class="references hanging-indent">

<div id="ref-H_Backman2016-bt">

H Backman, Tyler W, and Thomas Girke. 2016. “systemPipeR: NGS workflow and report generation environment.” *BMC Bioinformatics* 17 (1): 388. <https://doi.org/10.1186/s12859-016-1241-0>.

</div>

<div id="ref-Kaufmann2010-me">

Kaufmann, Kerstin, Frank Wellmer, Jose M Muiño, Thilia Ferrier, Samuel E Wuest, Vijaya Kumar, Antonio Serrano-Mislata, et al. 2010. “Orchestration of floral initiation by APETALA1.” *Science* 328 (5974): 85–89. <https://doi.org/10.1126/science.1185244>.

</div>

<div id="ref-Langmead2012-bs">

Langmead, Ben, and Steven L Salzberg. 2012. “Fast Gapped-Read Alignment with Bowtie 2.” *Nat. Methods* 9 (4): 357–59. <https://doi.org/10.1038/nmeth.1923>.

</div>

<div id="ref-Love2014-sh">

Love, Michael, Wolfgang Huber, and Simon Anders. 2014. “Moderated Estimation of Fold Change and Dispersion for RNA-seq Data with DESeq2.” *Genome Biol.* 15 (12): 550. <https://doi.org/10.1186/s13059-014-0550-8>.

</div>

<div id="ref-Robinson2010-uk">

Robinson, M D, D J McCarthy, and G K Smyth. 2010. “EdgeR: A Bioconductor Package for Differential Expression Analysis of Digital Gene Expression Data.” *Bioinformatics* 26 (1): 139–40. <https://doi.org/10.1093/bioinformatics/btp616>.

</div>

<div id="ref-Yu2015-xu">

Yu, Guangchuang, Li-Gen Wang, and Qing-Yu He. 2015. “ChIPseeker: An R/Bioconductor Package for ChIP Peak Annotation, Comparison and Visualization.” *Bioinformatics* 31 (14): 2382–3. <https://doi.org/10.1093/bioinformatics/btv145>.

</div>

<div id="ref-Zhang2008-pc">

Zhang, Y, T Liu, C A Meyer, J Eeckhoute, D S Johnson, B E Bernstein, C Nussbaum, et al. 2008. “Model-Based Analysis of ChIP-Seq (MACS).” *Genome Biol.* 9 (9). <https://doi.org/10.1186/gb-2008-9-9-r137>.

</div>

</div>
