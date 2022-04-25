---
title: "systemPipeR: Workflow Design and Reporting Environment" 
author: "Author: Daniela Cassol, Le Zhang and Thomas Girke"
date: "Last update: 25 April, 2022" 
output:
  BiocStyle::html_document:
    toc_float: true
    code_folding: show
package: systemPipeR
vignette: |
  %\VignetteEncoding{UTF-8}
  %\VignetteIndexEntry{systemPipeR: Workflow design and reporting generation environment}
  %\VignetteEngine{knitr::rmarkdown}
fontsize: 14pt
bibliography: bibtex.bib
editor_options: 
  chunk_output_type: console
weight: 7
type: docs
---

<script src="/rmarkdown-libs/htmlwidgets/htmlwidgets.js"></script>

<link href="/rmarkdown-libs/datatables-css/datatables-crosstalk.css" rel="stylesheet" />

<script src="/rmarkdown-libs/datatables-binding/datatables.js"></script>

<script src="/rmarkdown-libs/jquery/jquery-3.6.0.min.js"></script>

<link href="/rmarkdown-libs/dt-core/css/jquery.dataTables.min.css" rel="stylesheet" />
<link href="/rmarkdown-libs/dt-core/css/jquery.dataTables.extra.css" rel="stylesheet" />

<script src="/rmarkdown-libs/dt-core/js/jquery.dataTables.min.js"></script>

<link href="/rmarkdown-libs/dt-ext-fixedcolumns/css/fixedColumns.dataTables.min.css" rel="stylesheet" />

<script src="/rmarkdown-libs/dt-ext-fixedcolumns/js/dataTables.fixedColumns.min.js"></script>

<link href="/rmarkdown-libs/dt-ext-scroller/css/scroller.dataTables.min.css" rel="stylesheet" />

<script src="/rmarkdown-libs/dt-ext-scroller/js/dataTables.scroller.min.js"></script>

<link href="/rmarkdown-libs/crosstalk/css/crosstalk.min.css" rel="stylesheet" />

<script src="/rmarkdown-libs/crosstalk/js/crosstalk.min.js"></script>

<script src="/rmarkdown-libs/htmlwidgets/htmlwidgets.js"></script>

<link href="/rmarkdown-libs/datatables-css/datatables-crosstalk.css" rel="stylesheet" />

<script src="/rmarkdown-libs/datatables-binding/datatables.js"></script>

<script src="/rmarkdown-libs/jquery/jquery-3.6.0.min.js"></script>

<link href="/rmarkdown-libs/dt-core/css/jquery.dataTables.min.css" rel="stylesheet" />
<link href="/rmarkdown-libs/dt-core/css/jquery.dataTables.extra.css" rel="stylesheet" />

<script src="/rmarkdown-libs/dt-core/js/jquery.dataTables.min.js"></script>

<link href="/rmarkdown-libs/dt-ext-fixedcolumns/css/fixedColumns.dataTables.min.css" rel="stylesheet" />

<script src="/rmarkdown-libs/dt-ext-fixedcolumns/js/dataTables.fixedColumns.min.js"></script>

<link href="/rmarkdown-libs/dt-ext-scroller/css/scroller.dataTables.min.css" rel="stylesheet" />

<script src="/rmarkdown-libs/dt-ext-scroller/js/dataTables.scroller.min.js"></script>

<link href="/rmarkdown-libs/crosstalk/css/crosstalk.min.css" rel="stylesheet" />

<script src="/rmarkdown-libs/crosstalk/js/crosstalk.min.js"></script>

<script src="/rmarkdown-libs/kePrint/kePrint.js"></script>

<link href="/rmarkdown-libs/lightable/lightable.css" rel="stylesheet" />

<script src="/rmarkdown-libs/htmlwidgets/htmlwidgets.js"></script>

<link href="/rmarkdown-libs/vizjs/plotwf.css" rel="stylesheet" />

<script src="/rmarkdown-libs/vizjs/viz.js"></script>

<script src="/rmarkdown-libs/vizjs/full.render.js"></script>

<script src="/rmarkdown-libs/dom_to_image/dom_to_image.js"></script>

<link id="plotwf_legend-1-attachment" rel="attachment" href="systemPipeR_files/plotwf_legend/plotwf_legend.svg"/>

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
- Compile from command-line
Rscript -e "rmarkdown::render('systemPipeR.Rmd', c('BiocStyle::html_document'), clean=F); knitr::knit('systemPipeR.Rmd', tangle=TRUE)"
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
\[ [.Rmd](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/systempiper/systemPipeR.Rmd) \]    
\[ [.R](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/systempiper/systemPipeR.R) \]

</div>

### Introduction

[*`systemPipeR`*](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html)
provides flexible utilities for designing, building and running automated end-to-end
analysis workflows for a wide range of research applications, including
next-generation sequencing (NGS) experiments, such as RNA-Seq, ChIP-Seq,
VAR-Seq and Ribo-Seq (H Backman and Girke 2016). Important features include a uniform
workflow interface across different data analysis applications, automated
report generation, and support for running both R and command-line software,
such as NGS aligners or peak/variant callers, on local computers or compute
clusters (Figure 1). The latter supports interactive job submissions and batch
submissions to queuing systems of clusters.

*`systemPipeR`* has been designed to improve the reproducibility of large-scale data analysis
projects while substantially reducing the time it takes to analyze complex omics
data sets. It provides a uniform workflow interface and management system that allows
the users to run selected workflow steps, as well as customize and design entirely new workflows.
Additionally, the package take advantage of central community S4 classes of the Bioconductor
ecosystem, and enhances them with command-line software support.

The main motivation and advantages of using *`systemPipeR`* for complex data analysis tasks are:

1.  Design of complex workflows involving multiple R/Bioconductor packages
2.  Common workflow interface for different applications
3.  User-friendly access to widely used Bioconductor utilities
4.  Support of command-line software from within R
5.  Reduced complexity of using compute clusters from R
6.  Accelerated runtime of workflows via parallelization on computer systems with multiple CPU cores and/or multiple nodes
7.  Improved reproducibility by automating the generation of analysis reports

<center>

<img src="../utilities.png">

</center>

**Figure 1:** Relevant features in *`systemPipeR`*.
Workflow design concepts are illustrated under (A). Examples of `systemPipeR's` visualization functionalities are given under (B). </br>

A central concept for designing workflows within the *`systemPipeR`*
environment is the use of workflow management containers. Workflow management
containers facilitate the design and execution of complex data
analysis steps. For its command-line interface *`systemPipeR`* adopts the
widely used [Common Workflow Language](https://www.commonwl.org/) (CWL)
(Amstutz et al. 2016). The interface to CWL is established by *`systemPipeR's`*
workflow control class called *`SYSargsList`* (see Figure 2). This
design offers many advantages such as: (i) options to run workflows either
entirely from within R, from various command-line wrappers (e.g., *cwl-runner*)
or from other languages (*, e.g.,* Bash or Python). Apart from providing
support for both command-line and R/Bioconductor software, the package provides
utilities for containerization, parallel evaluations on computer clusters and
automated generation of interactive analysis reports.

<center>

<img src="../general.png">

</center>

**Figure 2:** Overview of `systemPipeR` workflows management instances. A) A
typical analysis workflow requires multiple software tools (red), as well the
description of the input (green) and output files, including analysis reports
(purple). B) `systemPipeR` provides multiple utilities to design and
build a workflow, allowing multi-instance, integration of R code and command-line
software, a simple and efficient annotation system, that allows automatic control
of the input and output data, and multiple resources to manage the entire workflow.
c) Options are provided to execute single or multiple workflow steps, while
enabling scalability, checkpoints, and generation of technical and scientific reports.

An important feature of *`systemPipeR's`* CWL interface is that it provides two
options to run command-line tools and workflows based on CWL. First, one can
run CWL in its native way via an R-based wrapper utility for *cwl-runner* or
*cwl-tools* (CWL-based approach). Second, one can run workflows using CWL’s
command-line and workflow instructions from within R (R-based approach). In the
latter case the same CWL workflow definition files (*e.g.* `*.cwl` and `*.yml`)
are used but rendered and executed entirely with R functions defined by
*`systemPipeR`*, and thus use CWL mainly as a command-line and workflow
definition format rather than execution software to run workflows. In this regard
*`systemPipeR`* also provides several convenience functions that are useful for
designing and debugging workflows, such as a command-line rendering function to
retrieve the exact command-line strings for each data set and processing step
prior to running a command-line.

### Workflow Management with *`SYSargsList`*

The `SYSargsList` S4 class is a list-like container that stores the paths to
all input and output files along with the corresponding parameters used in each
analysis step (see Figure 3). `SYSargsList` instances are constructed from an
optional targets files, and two CWL parameter files including `*.cwl` and
`*.yml` (for details, see below). When running preconfigured NGS workflows, the
only input the user needs to provide is the initial targets file containing the
paths to the input files (e.g., FASTQ) and experiment design information, such
as sample labels and biological replicates. Subsequent targets instances are
created automatically, based on the connectivity establish between each
workflow step. *`SYSargsList`* containers store all information required for
one or multiple steps. This establishes central control for running, monitoring
and debugging complex workflows from start to finish.

<center>

<img src="../sysargslist.png">

</center>

**Figure 3:** Workflow steps with input/output file operations are controlled by
the *`SYSargsList`* container. Each of its components (*`SYSargs2`*) are constructed
from an optional *targets* and two *param* files. Alternatively, *`LineWise`* instances
containing pure R code can be used.

## Getting Started

### Installation

The R software for running
[*`systemPipeR`*](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html)
can be downloaded from [*CRAN*](http://cran.at.r-project.org/). The
`systemPipeR` environment can be installed from the R console using the
[*`BiocManager::install`*](https://cran.r-project.org/web/packages/BiocManager/index.html)
command. The associated data package
[*`systemPipeRdata`*](http://www.bioconductor.org/packages/devel/data/experiment/html/systemPipeRdata.html)
can be installed the same way. The latter is a helper package for generating
`systemPipeR` workflow environments with a single command containing all
parameter files and sample data required to quickly test and run workflows.

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("systemPipeR")
BiocManager::install("systemPipeRdata")
```

Please note that if you desire to use a third-party command line tool, the
particular tool and dependencies need to be installed and executable.
See [details](#third-party-software-tools).

### Loading package and documentation

``` r
library("systemPipeR")  # Loads the package
library(help = "systemPipeR")  # Lists package info
vignette("systemPipeR")  # Opens vignette
```

### Load sample data and workflow templates

The mini sample FASTQ files used by this overview vignette as well as the
associated workflow reporting vignettes can be loaded via the
*`systemPipeRdata`* package as shown below. The chosen data set
[`SRP010938`](http://www.ncbi.nlm.nih.gov/sra/?term=SRP010938) obtains 18
paired-end (PE) read sets from *Arabidposis thaliana* (Howard et al. 2013). To
minimize processing time during testing, each FASTQ file has been subsetted to
90,000-100,000 randomly sampled PE reads that map to the first 100,000
nucleotides of each chromosome of the *A. thalina* genome. The corresponding
reference genome sequence (FASTA) and its GFF annotation files (provided in the
same download) have been truncated accordingly. This way the entire test sample
data set requires less than 200MB disk storage space. A PE read set has been
chosen for this test data set for flexibility, because it can be used for
testing both types of analysis routines requiring either SE (single-end) reads
or PE reads.

The following generates a fully populated *`systemPipeR`* workflow environment
(here for RNA-Seq) in the current working directory of an R session. At this time
the package includes workflow templates for RNA-Seq, ChIP-Seq, VAR-Seq, and Ribo-Seq.
Templates for additional NGS applications will be provided in the future.

``` r
systemPipeRdata::genWorkenvir(workflow = "rnaseq")
setwd("rnaseq")
```

### Project structure

The working environment of the sample data loaded in the previous step contains
the following pre-configured directory structure (Figure 4). Directory names are indicated
in <span style="color:grey">***green***</span>. Users can change this
structure as needed, but need to adjust the code in their workflows
accordingly.

  - <span style="color:green">***workflow/***</span> (*e.g.* *rnaseq/*)
      - This is the root directory of the R session running the workflow.
      - Run script ( *\*.Rmd*) and sample annotation (*targets.txt*) files are located here.
      - Note, this directory can have any name (*e.g.* <span style="color:green">***rnaseq***</span>, <span style="color:green">***varseq***</span>). Changing its name does not require any modifications in the run script(s).
      - **Important subdirectories**:
          - <span style="color:green">***param/***</span>
              - <span style="color:green">***param/cwl/***</span>: This subdirectory stores all the CWL parameter files. To organize workflows, each can have its own subdirectory, where all `CWL param` and `input.yml` files need to be in the same subdirectory.
          - <span style="color:green">***data/*** </span>
              - FASTQ files
              - FASTA file of reference (*e.g.* reference genome)
              - Annotation files
              - etc.
          - <span style="color:green">***results/***</span>
              - Analysis results are usually written to this directory, including: alignment, variant and peak files (BAM, VCF, BED); tabular result files; and image/plot files.
              - Note, the user has the option to organize results files for a given sample and analysis step in a separate subdirectory.

<center>

<img src="../directory.png">

</center>

**Figure 4:** *systemPipeR’s* preconfigured directory structure.

The following parameter files are included in each workflow template:

1.  *`targets.txt`*: initial one provided by user; downstream *`targets_*.txt`* files are generated automatically
2.  *`*.param/cwl`*: defines parameter for input/output file operations, *e.g.*:
      - *`hisat2/hisat2-mapping-se.cwl`*
      - *`hisat2/hisat2-mapping-se.yml`*
3.  Configuration files for computer cluster environments (skip on single machines):
      - *`.batchtools.conf.R`*: defines the type of scheduler for *`batchtools`* pointing to template file of cluster, and located in user’s home directory
      - *`*.tmpl`*: specifies parameters of scheduler used by a system, *e.g.* Torque, SGE, Slurm, etc.

### Structure of *`targets`* file

The *`targets`* file defines all input files (*e.g.* FASTQ, BAM, BCF) and
sample comparisons of an analysis workflow. The following shows the format of a
sample *`targets`* file included in the package. It also can be viewed and
downloaded from *`systemPipeR`*’s GitHub repository
[here](https://github.com/tgirke/systemPipeR/blob/master/inst/extdata/targets.txt).
In a target file with a single type of input files, here FASTQ files of
single-end (SE) reads, the first three columns are mandatory including their
column names, while it is four mandatory columns for FASTQ files of PE reads.
All subsequent columns are optional and any number of additional columns can be
added as needed. The columns in targets files are expected to be tab separated (TSV format).
The `SampleName` column contains usually short labels for
referencing samples (here FASTQ files) across many workflow steps (*e.g.*
plots and column titles). Importantly, the labels used in the `SampleName`
column need to be unique, while technical or biological replicates are
indicated by duplicated values under the `Factor` column. For readability
and transparency, it is useful to use here a short, consistent and informative
syntax for naming samples and replicates. To avoid problems with other
packages or external software, it is recommended to use the basic naming rules
for R objects and their components as outlined [here](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/rbasics/rbasics/#data-objects).
This is important since the values used under the `SampleName` and `Factor`
columns are intended to be used as labels for naming columns or plotting features
in downstream analysis steps.

Users should note here, the usage of targets files is optional when using
*systemPipeR’s* new CWL interface. They can be replaced by a standard YAML
input file used by CWL. Since for organizing experimental variables targets
files are extremely useful and user-friendly. Thus, we encourage users to keep using
them.

#### Structure of *`targets`* file for single-end (SE) samples

``` r
library(systemPipeR)
targetspath <- system.file("extdata", "targets.txt", package = "systemPipeR")
showDF(read.delim(targetspath, comment.char = "#"))
```

<div id="htmlwidget-1" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-1">{"x":{"filter":"none","vertical":false,"extensions":["FixedColumns","Scroller"],"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"],["./data/SRR446027_1.fastq.gz","./data/SRR446028_1.fastq.gz","./data/SRR446029_1.fastq.gz","./data/SRR446030_1.fastq.gz","./data/SRR446031_1.fastq.gz","./data/SRR446032_1.fastq.gz","./data/SRR446033_1.fastq.gz","./data/SRR446034_1.fastq.gz","./data/SRR446035_1.fastq.gz","./data/SRR446036_1.fastq.gz","./data/SRR446037_1.fastq.gz","./data/SRR446038_1.fastq.gz","./data/SRR446039_1.fastq.gz","./data/SRR446040_1.fastq.gz","./data/SRR446041_1.fastq.gz","./data/SRR446042_1.fastq.gz","./data/SRR446043_1.fastq.gz","./data/SRR446044_1.fastq.gz"],["M1A","M1B","A1A","A1B","V1A","V1B","M6A","M6B","A6A","A6B","V6A","V6B","M12A","M12B","A12A","A12B","V12A","V12B"],["M1","M1","A1","A1","V1","V1","M6","M6","A6","A6","V6","V6","M12","M12","A12","A12","V12","V12"],["Mock.1h.A","Mock.1h.B","Avr.1h.A","Avr.1h.B","Vir.1h.A","Vir.1h.B","Mock.6h.A","Mock.6h.B","Avr.6h.A","Avr.6h.B","Vir.6h.A","Vir.6h.B","Mock.12h.A","Mock.12h.B","Avr.12h.A","Avr.12h.B","Vir.12h.A","Vir.12h.B"],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],["23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>FileName<\/th>\n      <th>SampleName<\/th>\n      <th>Factor<\/th>\n      <th>SampleLong<\/th>\n      <th>Experiment<\/th>\n      <th>Date<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"scrollX":true,"fixedColumns":true,"deferRender":true,"scrollY":200,"scroller":true,"columnDefs":[{"className":"dt-right","targets":5},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

To work with custom data, users need to generate a *`targets`* file containing
the paths to their own FASTQ files and then provide under *`targetspath`* the
path to the corresponding *`targets`* file.

#### Structure of *`targets`* file for paired-end (PE) samples

For paired-end (PE) samples, the structure of the targets file is similar, where
users need to provide two FASTQ path columns: *`FileName1`* and *`FileName2`*
with the paths to the PE FASTQ files.

``` r
targetspath <- system.file("extdata", "targetsPE.txt", package = "systemPipeR")
showDF(read.delim(targetspath, comment.char = "#"))
```

<div id="htmlwidget-2" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-2">{"x":{"filter":"none","vertical":false,"extensions":["FixedColumns","Scroller"],"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"],["./data/SRR446027_1.fastq.gz","./data/SRR446028_1.fastq.gz","./data/SRR446029_1.fastq.gz","./data/SRR446030_1.fastq.gz","./data/SRR446031_1.fastq.gz","./data/SRR446032_1.fastq.gz","./data/SRR446033_1.fastq.gz","./data/SRR446034_1.fastq.gz","./data/SRR446035_1.fastq.gz","./data/SRR446036_1.fastq.gz","./data/SRR446037_1.fastq.gz","./data/SRR446038_1.fastq.gz","./data/SRR446039_1.fastq.gz","./data/SRR446040_1.fastq.gz","./data/SRR446041_1.fastq.gz","./data/SRR446042_1.fastq.gz","./data/SRR446043_1.fastq.gz","./data/SRR446044_1.fastq.gz"],["./data/SRR446027_2.fastq.gz","./data/SRR446028_2.fastq.gz","./data/SRR446029_2.fastq.gz","./data/SRR446030_2.fastq.gz","./data/SRR446031_2.fastq.gz","./data/SRR446032_2.fastq.gz","./data/SRR446033_2.fastq.gz","./data/SRR446034_2.fastq.gz","./data/SRR446035_2.fastq.gz","./data/SRR446036_2.fastq.gz","./data/SRR446037_2.fastq.gz","./data/SRR446038_2.fastq.gz","./data/SRR446039_2.fastq.gz","./data/SRR446040_2.fastq.gz","./data/SRR446041_2.fastq.gz","./data/SRR446042_2.fastq.gz","./data/SRR446043_2.fastq.gz","./data/SRR446044_2.fastq.gz"],["M1A","M1B","A1A","A1B","V1A","V1B","M6A","M6B","A6A","A6B","V6A","V6B","M12A","M12B","A12A","A12B","V12A","V12B"],["M1","M1","A1","A1","V1","V1","M6","M6","A6","A6","V6","V6","M12","M12","A12","A12","V12","V12"],["Mock.1h.A","Mock.1h.B","Avr.1h.A","Avr.1h.B","Vir.1h.A","Vir.1h.B","Mock.6h.A","Mock.6h.B","Avr.6h.A","Avr.6h.B","Vir.6h.A","Vir.6h.B","Mock.12h.A","Mock.12h.B","Avr.12h.A","Avr.12h.B","Vir.12h.A","Vir.12h.B"],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],["23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>FileName1<\/th>\n      <th>FileName2<\/th>\n      <th>SampleName<\/th>\n      <th>Factor<\/th>\n      <th>SampleLong<\/th>\n      <th>Experiment<\/th>\n      <th>Date<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"scrollX":true,"fixedColumns":true,"deferRender":true,"scrollY":200,"scroller":true,"columnDefs":[{"className":"dt-right","targets":6},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

#### Sample comparisons

Sample comparisons are defined in the header lines of the *`targets`* file
starting with ‘`# <CMP>`’.

``` r
readLines(targetspath)[1:4]
```

    ## [1] "# Project ID: Arabidopsis - Pseudomonas alternative splicing study (SRA: SRP010938; PMID: 24098335)"                                                                              
    ## [2] "# The following line(s) allow to specify the contrasts needed for comparative analyses, such as DEG identification. All possible comparisons can be specified with 'CMPset: ALL'."
    ## [3] "# <CMP> CMPset1: M1-A1, M1-V1, A1-V1, M6-A6, M6-V6, A6-V6, M12-A12, M12-V12, A12-V12"                                                                                             
    ## [4] "# <CMP> CMPset2: ALL"

The function *`readComp`* imports the comparison information and stores it in a
*`list`*. Alternatively, *`readComp`* can obtain the comparison information from
the corresponding *`SYSargsList`* object (see below). Note, these header lines are
optional. They are mainly useful for controlling comparative analyses according
to certain biological expectations, such as identifying differentially expressed
genes in RNA-Seq experiments based on simple pair-wise comparisons.

``` r
readComp(file = targetspath, format = "vector", delim = "-")
```

    ## $CMPset1
    ## [1] "M1-A1"   "M1-V1"   "A1-V1"   "M6-A6"   "M6-V6"   "A6-V6"   "M12-A12"
    ## [8] "M12-V12" "A12-V12"
    ## 
    ## $CMPset2
    ##  [1] "M1-A1"   "M1-V1"   "M1-M6"   "M1-A6"   "M1-V6"   "M1-M12"  "M1-A12" 
    ##  [8] "M1-V12"  "A1-V1"   "A1-M6"   "A1-A6"   "A1-V6"   "A1-M12"  "A1-A12" 
    ## [15] "A1-V12"  "V1-M6"   "V1-A6"   "V1-V6"   "V1-M12"  "V1-A12"  "V1-V12" 
    ## [22] "M6-A6"   "M6-V6"   "M6-M12"  "M6-A12"  "M6-V12"  "A6-V6"   "A6-M12" 
    ## [29] "A6-A12"  "A6-V12"  "V6-M12"  "V6-A12"  "V6-V12"  "M12-A12" "M12-V12"
    ## [36] "A12-V12"

## Project initialization

To create a workflow within *`systemPipeR`*, we can start by defining an empty
container and checking the directory structure:

``` r
sal <- SPRproject()
```

    ## Creating directory '/home/tgirke/tmp/GEN242/content/en/tutorials/systempiper/rnaseq/.SPRproject'
    ## Creating file '/home/tgirke/tmp/GEN242/content/en/tutorials/systempiper/rnaseq/.SPRproject/SYSargsList.yml'

Internally, `SPRproject` function will create a hidden folder called `.SPRproject`,
by default, to store all the log files.

In this stage, the object `sal` is a empty container, except for the project
information. The project information can be accessed by the `projectInfo` method:

``` r
sal
```

    ## Instance of 'SYSargsList': 
    ##  No workflow steps added

``` r
projectInfo(sal)
```

    ## $project
    ## [1] "/home/tgirke/tmp/GEN242/content/en/tutorials/systempiper/rnaseq"
    ## 
    ## $data
    ## [1] "data"
    ## 
    ## $param
    ## [1] "param"
    ## 
    ## $results
    ## [1] "results"
    ## 
    ## $logsDir
    ## [1] ".SPRproject"
    ## 
    ## $sysargslist
    ## [1] ".SPRproject/SYSargsList.yml"

Also, the `length` function will return how many steps this workflow contains,
and in this case, it is empty, as follow:

``` r
length(sal)
```

    ## [1] 0

## Workflow Design

*`systemPipeR`* workflows can be designed and built from start to finish with a
single command, or via importing an R Markdown file, or stepwise in interactive mode
from the R console.

In the next section, we will demonstrate how to build a workflow in interactive
mode, and in the following section how to build it from an R Markdown file.

New workflows are constructed, or existing ones modified, by connecting each step
via the `appendStep` method. Each `SYSargsList` instance contains instructions needed
for processing a set of input files with a specific command-line and the paths to
the corresponding outfiles generated.

The constructor function `Linewise` is used to construct the R code-based instructions.

### Build workflow interactively

This section demonstrates how to create a step in the workflow for running
the short read aligner HISAT2 (Kim, Langmead, and Salzberg 2015).

The constructor function renders the proper command-line strings for each sample
and software tool, appending a new step in the `SYSargsList` object defined in the
previous step. For that, the `SYSargsList` function requires data from three input files:

    - CWL command-line specification file (`wf_file` argument)
    - Input variables (`input_file` argument)
    - Targets file (`targets` argument)

In CWL, files with the extension *`.cwl`* define the parameters of a chosen
command-line step or an entire workflow, while files with the extension *`.yml`* define
the input variables of command-line steps. Note, input variables provided
by a *targets* or *targetsPE* file can be passed on to a *`SYSargsList`* instance via the *inputvars*
argument of the *SYSargsList* function.

``` r
appendStep(sal) <- SYSargsList(step_name = "hisat2_mapping", dir = TRUE, targets = "targetsPE.txt",
    wf_file = "workflow-hisat2/workflow_hisat2-pe.cwl", input_file = "workflow-hisat2/workflow_hisat2-pe.yml",
    dir_path = "param/cwl", inputvars = c(FileName1 = "_FASTQ_PATH1_", FileName2 = "_FASTQ_PATH2_",
        SampleName = "_SampleName_"))
```

An overview of the workflow can be printed as follows.

``` r
sal
```

    ## Instance of 'SYSargsList': 
    ##     WF Steps:
    ##        1. hisat2_mapping --> Status: Pending 
    ##            Total Files: 72 | Existing: 0 | Missing: 72 
    ##          1.1. hisat2
    ##              cmdlist: 18 | Pending: 18
    ##          1.2. samtools-view
    ##              cmdlist: 18 | Pending: 18
    ##          1.3. samtools-sort
    ##              cmdlist: 18 | Pending: 18
    ##          1.4. samtools-index
    ##              cmdlist: 18 | Pending: 18
    ## 

Note that the workflow status is *Pending*, which means the workflow object has
been constructed in R. However, we did not execute the workflow yet.

Several accessor methods are available to explore the *`SYSargsList`* object.

Of particular interest is the *`cmdlist()`* method. It constructs the system
commands for running command-line software as specified by a given *`.cwl`*
file combined with the paths to the input samples (*e.g.* FASTQ files), here provided
by a *`targets`* file. The example below shows the *`cmdlist()`* output for
running HISAT2 on the first PE read sample. Evaluating the output of
*`cmdlist()`* can be very helpful for designing and debugging *`.cwl`* files
of new command-line software or changing the parameter settings of existing
ones.

The rendered command-line instance for each input sample can be returned as follows.

``` r
cmdlist(sal, step = "hisat2_mapping", targets = 1)
```

    ## $hisat2_mapping
    ## $hisat2_mapping$M1A
    ## $hisat2_mapping$M1A$hisat2
    ## [1] "hisat2 -S ./results/M1A.sam  -x ./data/tair10.fasta  -k 1  --min-intronlen 30  --max-intronlen 3000  -1 ./data/SRR446027_1.fastq.gz -2 ./data/SRR446027_2.fastq.gz --threads 4"
    ## 
    ## $hisat2_mapping$M1A$`samtools-view`
    ## [1] "samtools view -bS -o ./results/M1A.bam  ./results/M1A.sam "
    ## 
    ## $hisat2_mapping$M1A$`samtools-sort`
    ## [1] "samtools sort -o ./results/M1A.sorted.bam  ./results/M1A.bam  -@ 4"
    ## 
    ## $hisat2_mapping$M1A$`samtools-index`
    ## [1] "samtools index -b results/M1A.sorted.bam  results/M1A.sorted.bam.bai  ./results/M1A.sorted.bam "

The outfiles components of *`SYSargsList`* define the expected output files for
each step in the workflow; some of which are the input for the next workflow step,
here next *`SYSargsList`* instance (see Figure 3).

``` r
outfiles(sal)
```

    ## $hisat2_mapping
    ## DataFrame with 18 rows and 4 columns
    ##              hisat2_sam       samtools_bam      samtools_sort_bam
    ##             <character>        <character>            <character>
    ## M1A   ./results/M1A.sam  ./results/M1A.bam ./results/M1A.sorted..
    ## M1B   ./results/M1B.sam  ./results/M1B.bam ./results/M1B.sorted..
    ## A1A   ./results/A1A.sam  ./results/A1A.bam ./results/A1A.sorted..
    ## A1B   ./results/A1B.sam  ./results/A1B.bam ./results/A1B.sorted..
    ## V1A   ./results/V1A.sam  ./results/V1A.bam ./results/V1A.sorted..
    ## ...                 ...                ...                    ...
    ## M12B ./results/M12B.sam ./results/M12B.bam ./results/M12B.sorte..
    ## A12A ./results/A12A.sam ./results/A12A.bam ./results/A12A.sorte..
    ## A12B ./results/A12B.sam ./results/A12B.bam ./results/A12B.sorte..
    ## V12A ./results/V12A.sam ./results/V12A.bam ./results/V12A.sorte..
    ## V12B ./results/V12B.sam ./results/V12B.bam ./results/V12B.sorte..
    ##              samtools_index
    ##                 <character>
    ## M1A  ./results/M1A.sorted..
    ## M1B  ./results/M1B.sorted..
    ## A1A  ./results/A1A.sorted..
    ## A1B  ./results/A1B.sorted..
    ## V1A  ./results/V1A.sorted..
    ## ...                     ...
    ## M12B ./results/M12B.sorte..
    ## A12A ./results/A12A.sorte..
    ## A12B ./results/A12B.sorte..
    ## V12A ./results/V12A.sorte..
    ## V12B ./results/V12B.sorte..

In an ‘R-centric’ rather than a ‘CWL-centric’ workflow design the connectivity
among workflow steps is established by creating the downstream targets files
automatically (see Figure 3). Each step uses the outfiles from the previous step as
input, thereby establishing connectivity among each step in the workflow. By chaining
several *`SYSargsList`* steps together one can construct complex workflows involving
many sample-level input/output file operations with any combination of command-line or R-based
software. Also, *`systemPipeR`* provides features to automatically build these
connections. This provides additional security to make sure that all samples will be
processed in a reproducible manner.

Alternatively, a CWL-centric workflow design can be used that defines
all or most workflow steps with CWL workflow and parameter files. Due to time and
space restrictions, the CWL-centric approach is not covered by this tutorial.

The following code generates a `data.frame` containing important read alignment statistics
such as the total number of reads in the FASTQ files, the number of total alignments,
as well as the number of primary alignments in the corresponding BAM files.

The constructor function `LineWise` requires the `step_name` and the R code
assigned to the `code` argument.
The R code should be enclosed by braces (`{}`) followed by new lines.

``` r
appendStep(sal) <- LineWise(code = {
    fqpaths <- getColumn(sal, step = "hisat2_mapping", "targetsWF", column = "FileName1")
    bampaths <- getColumn(sal, step = "hisat2_mapping", "outfiles", column = "samtools_sort_bam")
    read_statsDF <- alignStats(args = bampaths, fqpaths = fqpaths, pairEnd = TRUE)
    write.table(read_statsDF, "results/alignStats.xls", row.names = FALSE, quote = FALSE,
        sep = "\t")
}, step_name = "align_stats", dependency = "hisat2_mapping")
```

To inspect the R code of this step, the `codeLine` method can be used.

``` r
codeLine(sal, "align_stats")
```

    ## align_stats
    ##     fqpaths <- getColumn(sal, step = "hisat2_mapping", "targetsWF", column = "FileName1")
    ##     bampaths <- getColumn(sal, step = "hisat2_mapping", "outfiles", column = "samtools_sort_bam")
    ##     read_statsDF <- alignStats(args = bampaths, fqpaths = fqpaths, pairEnd = TRUE)
    ##     write.table(read_statsDF, "results/alignStats.xls", row.names = FALSE, quote = FALSE, sep = "\t")

The `getColumn` method allows to extract certain information of a workflow instance, here
output files that can be accessed from within R as demonstrated below.

``` r
getColumn(sal, step = "hisat2_mapping", "outfiles", column = "samtools_sort_bam")
```

    ##                         M1A                         M1B 
    ##  "./results/M1A.sorted.bam"  "./results/M1B.sorted.bam" 
    ##                         A1A                         A1B 
    ##  "./results/A1A.sorted.bam"  "./results/A1B.sorted.bam" 
    ##                         V1A                         V1B 
    ##  "./results/V1A.sorted.bam"  "./results/V1B.sorted.bam" 
    ##                         M6A                         M6B 
    ##  "./results/M6A.sorted.bam"  "./results/M6B.sorted.bam" 
    ##                         A6A                         A6B 
    ##  "./results/A6A.sorted.bam"  "./results/A6B.sorted.bam" 
    ##                         V6A                         V6B 
    ##  "./results/V6A.sorted.bam"  "./results/V6B.sorted.bam" 
    ##                        M12A                        M12B 
    ## "./results/M12A.sorted.bam" "./results/M12B.sorted.bam" 
    ##                        A12A                        A12B 
    ## "./results/A12A.sorted.bam" "./results/A12B.sorted.bam" 
    ##                        V12A                        V12B 
    ## "./results/V12A.sorted.bam" "./results/V12B.sorted.bam"

### Build workflow from an R Markdown

The workflow can be created by importing the steps from an R Markdown file.
As demonstrated above, it is required to initialize the project with the `SPRproject` function.

The `importWF` function will scan and import all R chunks from the R Markdown file
and build and store the workflow in a `SYSargsList` object (here named `sal`).

``` r
sal <- SPRproject()
```

    ## Creating directory '/home/tgirke/tmp/GEN242/content/en/tutorials/systempiper/rnaseq/.SPRproject'
    ## Creating file '/home/tgirke/tmp/GEN242/content/en/tutorials/systempiper/rnaseq/.SPRproject/SYSargsList.yml'

``` r
sal <- importWF(sal, file_path = "systemPipeRNAseq.Rmd")
```

    ## Reading Rmd file

    ## 
    ##  ---- Actions ----

    ## Ignore none-R chunks at line: 25

    ## Checking chunk eval values

    ## Checking chunk SPR option

    ## Ignore non-SPR chunks: 34, 43, 61, 96, 121, 188, 205, 298, 658, 686, 704, 712, 723

    ## Parse chunk code

    ## Now importing step 'load_SPR' 
    ## Now importing step 'preprocessing' 
    ## Now importing step 'trimming' 
    ## Now importing step 'fastq_report' 
    ## Now importing step 'hisat2_index' 
    ## Now importing step 'hisat2_mapping' 
    ## Now importing step 'align_stats' 
    ## Now importing step 'bam_IGV' 
    ## Now importing step 'create_db' 
    ## Now importing step 'read_counting' 
    ## Now importing step 'sample_tree' 
    ## Now importing step 'run_edger' 
    ## Now importing step 'custom_annot' 
    ## Now importing step 'filter_degs' 
    ## Now importing step 'venn_diagram' 
    ## Now importing step 'get_go_annot' 
    ## Now importing step 'go_enrich' 
    ## Now importing step 'go_plot' 
    ## Now importing step 'heatmap' 
    ## Now importing step 'sessionInfo'

Several accessor functions are available to inspect the workflow.

``` r
stepsWF(sal)
dependency(sal)
codeLine(sal)
targetsWF(sal)
```

#### Third-party software tools

Currently, *systemPipeR* provides a collection of preconfigured *`param`* file templates  
for third-party software tools. A sub-selection of supported software is provided in the  
following table.

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:500px; overflow-x: scroll; width:100%; ">

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;">

Tool Name

</th>

<th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;">

Description

</th>

<th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;">

Step

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:center;">

<a href="http://bio-bwa.sourceforge.net/bwa.shtml">bwa</a>

</td>

<td style="text-align:center;">

BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. 

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #8FBC8F !important;">Alignment</span>

</td>

</tr>

<tr>

<td style="text-align:center;">

<a href="http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml">Bowtie2</a>

</td>

<td style="text-align:center;">

Bowtie 2 is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences.

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #8FBC8F !important;">Alignment</span>

</td>

</tr>

<tr>

<td style="text-align:center;">

<a href="http://hannonlab.cshl.edu/fastx_toolkit/commandline.html">FASTX-Toolkit</a>

</td>

<td style="text-align:center;">

FASTX-Toolkit is a collection of command line tools for Short-Reads FASTA/FASTQ files preprocessing.

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #EC7770 !important;">Read Preprocessing</span>

</td>

</tr>

<tr>

<td style="text-align:center;">

<a href="http://hibberdlab.com/transrate/">TransRate</a>

</td>

<td style="text-align:center;">

Transrate is software for de-novo transcriptome assembly quality analysis.

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #D98576 !important;">Quality</span>

</td>

</tr>

<tr>

<td style="text-align:center;">

<a href="http://research-pub.gene.com/gmap/">Gsnap</a>

</td>

<td style="text-align:center;">

GSNAP is a genomic short-read nucleotide alignment program.

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #8FBC8F !important;">Alignment</span>

</td>

</tr>

<tr>

<td style="text-align:center;">

<a href="http://www.htslib.org/doc/samtools-1.2.html">Samtools</a>

</td>

<td style="text-align:center;">

Samtools is a suite of programs for interacting with high-throughput sequencing data.

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #D08C79 !important;">Post-processing</span>

</td>

</tr>

<tr>

<td style="text-align:center;">

<a href="http://www.usadellab.org/cms/?page=trimmomatic">Trimmomatic</a>

</td>

<td style="text-align:center;">

Trimmomatic is a flexible read trimming tool for Illumina NGS data.

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #EC7770 !important;">Read Preprocessing</span>

</td>

</tr>

<tr>

<td style="text-align:center;">

<a href="https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf">Rsubread</a>

</td>

<td style="text-align:center;">

Rsubread is a Bioconductor software package that provides high-performance alignment and read counting functions for RNA-seq reads.

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #8FBC8F !important;">Alignment</span>

</td>

</tr>

<tr>

<td style="text-align:center;">

<a href="https://broadinstitute.github.io/picard/">Picard</a>

</td>

<td style="text-align:center;">

Picard is a set of command line tools for manipulating high-throughput sequencing (HTS) data and formats such as SAM/BAM/CRAM and VCF.

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #B4A082 !important;">Manipulating HTS data</span>

</td>

</tr>

<tr>

<td style="text-align:center;">

<a href="https://busco.ezlab.org/">Busco</a>

</td>

<td style="text-align:center;">

BUSCO assesses genome assembly and annotation completeness with Benchmarking Universal Single-Copy Orthologs.

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #D98576 !important;">Quality</span>

</td>

</tr>

<tr>

<td style="text-align:center;">

<a href="https://ccb.jhu.edu/software/hisat2/manual.shtml">Hisat2</a>

</td>

<td style="text-align:center;">

HISAT2 is a fast and sensitive alignment program for mapping NGS reads (both DNA and RNA) to reference genomes.

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #8FBC8F !important;">Alignment</span>

</td>

</tr>

<tr>

<td style="text-align:center;">

<a href="https://ccb.jhu.edu/software/tophat/manual.shtml">Tophat2</a>

</td>

<td style="text-align:center;">

TopHat is a fast splice junction mapper for RNA-Seq reads.

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #8FBC8F !important;">Alignment</span>

</td>

</tr>

<tr>

<td style="text-align:center;">

<a href="https://gatk.broadinstitute.org/hc/en-us">GATK</a>

</td>

<td style="text-align:center;">

Variant Discovery in High-Throughput Sequencing Data.

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FF6A6A !important;">Variant Discovery</span>

</td>

</tr>

<tr>

<td style="text-align:center;">

<a href="https://github.com/alexdobin/STAR">STAR</a>

</td>

<td style="text-align:center;">

STAR is an ultrafast universal RNA-seq aligner.

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #8FBC8F !important;">Alignment</span>

</td>

</tr>

<tr>

<td style="text-align:center;">

<a href="https://github.com/FelixKrueger/TrimGalore">Trim\_galore</a>

</td>

<td style="text-align:center;">

Trim Galore is a wrapper around Cutadapt and FastQC to consistently apply adapter and quality trimming to FastQ files.

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #EC7770 !important;">Read Preprocessing</span>

</td>

</tr>

<tr>

<td style="text-align:center;">

<a href="https://github.com/TransDecoder/TransDecoder/wiki">TransDecoder</a>

</td>

<td style="text-align:center;">

TransDecoder identifies candidate coding regions within transcript sequences.

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #ABA785 !important;">Find Coding Regions</span>

</td>

</tr>

<tr>

<td style="text-align:center;">

<a href="https://github.com/trinityrnaseq/trinityrnaseq/wiki">Trinity</a>

</td>

<td style="text-align:center;">

Trinity assembles transcript sequences from Illumina RNA-Seq data.

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #A1AE88 !important;">denovo Transcriptome Assembly</span>

</td>

</tr>

<tr>

<td style="text-align:center;">

<a href="https://github.com/Trinotate/Trinotate.github.io/wiki">Trinotate</a>

</td>

<td style="text-align:center;">

Trinotate is a comprehensive annotation suite designed for automatic functional annotation of transcriptomes.

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #F5706D !important;">Transcriptome Functional Annotation</span>

</td>

</tr>

<tr>

<td style="text-align:center;">

<a href="https://macs3-project.github.io/MACS/">MACS2</a>

</td>

<td style="text-align:center;">

MACS2 identifies transcription factor binding sites in ChIP-seq data.

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #C7937C !important;">Peak calling</span>

</td>

</tr>

<tr>

<td style="text-align:center;">

<a href="https://pachterlab.github.io/kallisto/manual">Kallisto</a>

</td>

<td style="text-align:center;">

kallisto is a program for quantifying abundances of transcripts from RNA-Seq data.

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #E37E73 !important;">Read counting</span>

</td>

</tr>

<tr>

<td style="text-align:center;">

<a href="https://samtools.github.io/bcftools/howtos/index.html">BCFtools</a>

</td>

<td style="text-align:center;">

BCFtools is a program for variant calling and manipulating files in the Variant Call Format (VCF) and its binary counterpart BCF.

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FF6A6A !important;">Variant Discovery</span>

</td>

</tr>

<tr>

<td style="text-align:center;">

<a href="https://www.bioinformatics.babraham.ac.uk/projects/bismark/">Bismark</a>

</td>

<td style="text-align:center;">

Bismark is a program to map bisulfite treated sequencing reads to a genome of interest and perform methylation calls in a single step.

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #98B58B !important;">Bisulfite mapping</span>

</td>

</tr>

<tr>

<td style="text-align:center;">

<a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/">Fastqc</a>

</td>

<td style="text-align:center;">

FastQC is a quality control tool for high throughput sequence data.

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #D98576 !important;">Quality</span>

</td>

</tr>

<tr>

<td style="text-align:center;">

<a href="https://www.ncbi.nlm.nih.gov/books/NBK279690/">Blast</a>

</td>

<td style="text-align:center;">

BLAST finds regions of similarity between biological sequences.

</td>

<td style="text-align:center;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #BD997F !important;">Blast</span>

</td>

</tr>

</tbody>

</table>

</div>

Importantly, command-line software needs to be installed on a user’s system  
and available in a user’s PATH. To check whether this is the case, one can use the  
`tryCL` function.

``` r
tryCL(command = "grep")
```

## How to run a Workflow

### Setup and requirements

To execute the code in this tutorial, one needs the following software installed.

  - R (version \>=4.1.2)
  - systemPipeR package (version \>=2.0.8)
  - Hisat2 (version \>= 2.1.0)

To build custom workflows with any additional command-line software, one needs
the respective software installed and available in the user PATH. To test this,
one can use the `tryCL` function.

``` r
tryCL(command = "hisat2")  ## 'All set up, proceed!'
```

### Running the workflow

For running a workflow, the `runWF` function can be used. I executes all command-line
calls stored in the workflow container.

``` r
sal <- runWF(sal)
```

This essential function allows the user to choose one or multiple steps to be
executed using the `steps` argument. However, it is necessary to follow the
workflow dependency graph. If a selected step depends on a previous step(s) that
was not executed, then execution will fail.

``` r
sal <- runWF(sal, steps = c(1, 3))
```

The `runWF` function allows to force the execution of chosen workflow steps, even if
the status of a given step is `'Success'` and the expected `outfiles` exist already.
This is useful for updating certain steps when needed. Another option is to
force `runWF` to ignore all warnings and errors. This can be achieved by assigning
`FALSE` to the arguments `warning.stop` and `error.stop`, respectively.

``` r
sal <- runWF(sal, force = TRUE, warning.stop = FALSE, error.stop = FALSE)
```

#### Parallelization on clusters

This section of the tutorial provides an introduction to the usage of *`systemPipeR`*
on computer clusters.

The computation of time-consuming steps can be greatly accelerated by processing many files
in parallel using several compute nodes of a cluster, where a scheduling
system is used for load balancing.

The `resources` list object provides the number of independent parallel cluster
processes defined under the `Njobs` element in the list. The following example
will run 18 processes in parallel using for each 4 CPU cores.
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
add these resources later with the `addResources` function.

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
sal <- addResources(sal, c("hisat2_mapping"), resources = resources)
sal <- runWF(sal)
```

Note: The above example is submitting the job to partition called `short`. Users need
to adjust this parameter to their environment.

### Visualize workflow

*`systemPipeR`* workflows instances can be visualized with the `plotWF` function.
The resulting plot includes the following information.

  - Workflow topology graph (dependency graphs between different steps)
      - Workflow step status, *e.g.* `Success`, `Error`, `Pending`, `Warnings`
      - Sample status and statistics
      - Workflow timing: running duration time

If no argument is provided, the basic plot will automatically detect width,
height, layout, plot method, branches, *etc*.

``` r
plotWF(sal)
```

<div id="htmlwidget-3" style="width:672px;height:480px;" class="plotwf html-widget"></div>
<script type="application/json" data-for="htmlwidget-3">{"x":{"dot":"digraph {\n    node[fontsize=20];\n    subgraph {\n        load_SPR -> hisat2_index -> hisat2_mapping -> create_db -> read_counting -> run_edger -> custom_annot -> filter_degs -> get_go_annot -> go_enrich -> heatmap -> sessionInfo\n   }\n    load_SPR -> preprocessing\n    preprocessing -> hisat2_mapping\n    hisat2_mapping -> align_stats\n    hisat2_mapping -> bam_IGV\n    read_counting -> sample_tree\n    filter_degs -> venn_diagram\n    go_enrich -> go_plot\n    load_SPR -> trimming\n    preprocessing -> fastq_report\n    load_SPR[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">load_SPR<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step load_SPR: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-24 16:55:24; End time: 2022-04-24 16:55:24; Duration: 00:00:00\"]\n    preprocessing[fillcolor=\"#d3d6eb\" style=\"filled, rounded\" label=<<b><font color=\"black\">preprocessing<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">36<\/font><\/b>; <font color=\"black\">0s<\/font>> , shape=\"box\"  tooltip=\"step preprocessing: 0 samples passed; 0 samples have warnings; 0 samples have errors; 36 samples in total; Start time: 2022-04-24 16:55:24; End time: 2022-04-24 16:55:24; Duration: 00:00:00\"]\n    trimming[style=\"solid, rounded\" label=<<b><font color=\"black\">trimming<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">72<\/font><\/b>; <font color=\"black\">0s<\/font>> , shape=\"box\"  tooltip=\"step trimming: 0 samples passed; 0 samples have warnings; 0 samples have errors; 72 samples in total; Start time: 2022-04-24 16:55:24; End time: 2022-04-24 16:55:24; Duration: 00:00:00\"]\n    fastq_report[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">fastq_report<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step fastq_report: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-24 16:55:24; End time: 2022-04-24 16:55:24; Duration: 00:00:00\"]\n    hisat2_index[fillcolor=\"#d3d6eb\" style=\"filled, rounded\" label=<<b><font color=\"black\">hisat2_index<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">8<\/font><\/b>; <font color=\"black\">0s<\/font>> , shape=\"box\"  tooltip=\"step hisat2_index: 0 samples passed; 0 samples have warnings; 0 samples have errors; 8 samples in total; Start time: 2022-04-24 16:55:24; End time: 2022-04-24 16:55:24; Duration: 00:00:00\"]\n    hisat2_mapping[fillcolor=\"#d3d6eb\" style=\"filled, rounded\" label=<<b><font color=\"black\">hisat2_mapping<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">72<\/font><\/b>; <font color=\"black\">0s<\/font>> , shape=\"box\"  tooltip=\"step hisat2_mapping: 0 samples passed; 0 samples have warnings; 0 samples have errors; 72 samples in total; Start time: 2022-04-24 16:55:24; End time: 2022-04-24 16:55:24; Duration: 00:00:00\"]\n    align_stats[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">align_stats<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step align_stats: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-24 16:55:24; End time: 2022-04-24 16:55:24; Duration: 00:00:00\"]\n    bam_IGV[style=\"solid, \"label=<<b><font color=\"black\">bam_IGV<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step bam_IGV: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-24 16:55:24; End time: 2022-04-24 16:55:24; Duration: 00:00:00\"]\n    create_db[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">create_db<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step create_db: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-24 16:55:24; End time: 2022-04-24 16:55:24; Duration: 00:00:00\"]\n    read_counting[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">read_counting<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step read_counting: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-24 16:55:24; End time: 2022-04-24 16:55:24; Duration: 00:00:00\"]\n    sample_tree[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">sample_tree<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step sample_tree: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-24 16:55:24; End time: 2022-04-24 16:55:24; Duration: 00:00:00\"]\n    run_edger[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">run_edger<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step run_edger: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-24 16:55:24; End time: 2022-04-24 16:55:24; Duration: 00:00:00\"]\n    custom_annot[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">custom_annot<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step custom_annot: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-24 16:55:24; End time: 2022-04-24 16:55:24; Duration: 00:00:00\"]\n    filter_degs[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">filter_degs<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step filter_degs: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-24 16:55:24; End time: 2022-04-24 16:55:24; Duration: 00:00:00\"]\n    venn_diagram[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">venn_diagram<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step venn_diagram: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-24 16:55:24; End time: 2022-04-24 16:55:24; Duration: 00:00:00\"]\n    get_go_annot[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">get_go_annot<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step get_go_annot: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-24 16:55:24; End time: 2022-04-24 16:55:24; Duration: 00:00:00\"]\n    go_enrich[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">go_enrich<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step go_enrich: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-24 16:55:24; End time: 2022-04-24 16:55:24; Duration: 00:00:00\"]\n    go_plot[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">go_plot<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step go_plot: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-24 16:55:24; End time: 2022-04-24 16:55:24; Duration: 00:00:00\"]\n    heatmap[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">heatmap<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step heatmap: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-24 16:55:24; End time: 2022-04-24 16:55:24; Duration: 00:00:00\"]\n    sessionInfo[fillcolor=\"#d3d6eb\" style=\"filled, \"label=<<b><font color=\"black\">sessionInfo<\/font><br><\/br><font color=\"#5cb85c\">0<\/font>/<font color=\"#f0ad4e\">0<\/font>/<font color=\"#d9534f\">0<\/font>/<font color=\"blue\">1<\/font><\/b>; <font color=\"black\">0s<\/font>>  tooltip=\"step sessionInfo: 0 samples passed; 0 samples have warnings; 0 samples have errors; 1 samples in total; Start time: 2022-04-24 16:55:24; End time: 2022-04-24 16:55:24; Duration: 00:00:00\"]\n        subgraph cluster_legend {\n        rankdir=TB;\n        color=\"#eeeeee\";\n        style=filled;\n        ranksep =1;\n        label=\"Legends\";\n        fontsize = 30;\n        node [style=filled, fontsize=10];\n        legend_img-> step_state[color=\"#eeeeee\"];\n\n        legend_img[shape=none, image=\"plotwf_legend-src.png\", label = \" \", height=1, width=3, style=\"\"];\n\n        step_state[style=\"filled\", shape=\"box\" color=white, label =<\n            <table>\n            <tr><td><b>Step Colors<\/b><\/td><\/tr>\n            <tr><td><font color=\"black\">Pending steps<\/font>; <font color=\"#5cb85c\">Successful steps<\/font>; <font color=\"#d9534f\">Failed steps<\/font><\/td><\/tr>\n            <tr><td><b>Targets Files / Code Chunk <\/b><\/td><\/tr><tr><td><font color=\"#5cb85c\">0 (pass) <\/font> | <font color=\"#f0ad4e\">0 (warning) <\/font> | <font color=\"#d9534f\">0 (error) <\/font> | <font color=\"blue\">0 (total)<\/font>; Duration<\/td><\/tr><\/table>\n            >];\n\n    }\n\n}\n","plotid":"sprwf-65123748","responsive":true,"width":null,"height":null,"plot_method":"renderSVGElement","rmd":true,"msg":"","plot_ctr":true,"pan_zoom":false,"legend_uri":"data:image/svg+xml;base64,PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iVVRGLTgiPz4KPCEtLSBEbyBub3QgZWRpdCB0aGlzIGZpbGUgd2l0aCBlZGl0b3JzIG90aGVyIHRoYW4gZGlhZ3JhbXMubmV0IC0tPgo8IURPQ1RZUEUgc3ZnIFBVQkxJQyAiLS8vVzNDLy9EVEQgU1ZHIDEuMS8vRU4iICJodHRwOi8vd3d3LnczLm9yZy9HcmFwaGljcy9TVkcvMS4xL0RURC9zdmcxMS5kdGQiPgo8c3ZnIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyIgeG1sbnM6eGxpbms9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGxpbmsiIHZlcnNpb249IjEuMSIgd2lkdGg9IjQ5NnB4IiBoZWlnaHQ9IjI3OHB4IiB2aWV3Qm94PSItMC41IC0wLjUgNDk2IDI3OCIgY29udGVudD0iJmx0O214ZmlsZSBob3N0PSZxdW90O2FwcC5kaWFncmFtcy5uZXQmcXVvdDsgbW9kaWZpZWQ9JnF1b3Q7MjAyMS0xMS0yNFQyMDozOTo0NC45MjNaJnF1b3Q7IGFnZW50PSZxdW90OzUuMCAoWDExOyBMaW51eCB4ODZfNjQpIEFwcGxlV2ViS2l0LzUzNy4zNiAoS0hUTUwsIGxpa2UgR2Vja28pIENocm9tZS85Ni4wLjQ2NjQuNDUgU2FmYXJpLzUzNy4zNiZxdW90OyB2ZXJzaW9uPSZxdW90OzE1LjguMyZxdW90OyBldGFnPSZxdW90O1RfZEV4bkw3U0xYMmJ1WDhQYVgzJnF1b3Q7IHR5cGU9JnF1b3Q7Z29vZ2xlJnF1b3Q7Jmd0OyZsdDtkaWFncmFtIGlkPSZxdW90O1Z3OVZWaUdzeC1zTF9OaktlN0pQJnF1b3Q7Jmd0OzdWcmJjdUk0RVAwYUhwT3lKZDk0REFSbVhtWjNLMnpWUGl0WUdGVmt5V1BMZ2N6WHIyVEwrQ0lEM3NHUVRSVWtsZGd0cVMyZDAycDF0NW5BZWJ6L2xxSmsrNE9IbUU2QUZlNG44SGtDZ0EwQW1LaGZLL3pRRXN2eFMwbVVrbERMYXNHSy9NSlZSeTNOU1lpelZrZkJPUlVrYVF2WG5ERzhGaTBaU2xPK2EzZmJjTnArYW9JaWJBaFdhMFJONlQ4a0ZOdFNPcldzV3Y0ZGsyaGJQZG1yV21KVWRkYUNiSXRDdm11STRHSUM1eW5ub3J5SzkzTk1GWG9WTHVXNDVaSFd3OFJTek1TUUFacUpkMFJ6dlRZOUwvRlJMVFpLZVo1TTRFeitZeUZXNHl4NWQ1aTR1dEZLY0Nyd3ZvOEE5Rm9wczh3SjJvZGxTNFBCUE1ZaS9aQmRLa1dCSHFKdEJVQjl2NnVCQjQ1YnlyWk4wS0duQ2Rka1J3ZmROUjd5UWtQU0R3ODhEMDhibHQyV0NMeEswRnExN3FUNVM5bFd4Rkwvc3kwdk40VFNPYWM4TGNiQzBNVkI2Q2cwUmNyZmNOWENPSlBEWnhGRldkWUhkL2FHeFhwN0R2c214cUFmWTQycEE0ZEI2b3lBcUhOZFJKZkxSZkM4T0lib1ZWRU1wZ05SOUM5SDBiMHVpb3ZsRWl6bnQwWHhMR3p3Y3RpOHNXSGpUT2pUS1RCUjNManFSOG9SSlJHVE1vbzM0aGlvU2xWamJQbVJjaTRmVG9ReXNXQkVoSDNyU2c3VE54RE9PSldOWFpqbDlFVWJ5d3FrdFZ3R2xpak0xQ0tKUEhPZmRFTk13bEFObjZVNEk3LzBrYUl3U1RoaG9waXpPNXU0ejBwWExuaFc4bUlmQmJ4QlZvT0RKcDN5ZG9saVFoWDZmNU5ZUmh6QStnUHY1TjhYSGlNMmxJN1RCeHp3L1VlM2ZjUUIwNVBBSHNLcXcvUVN2Z0tEcnhCbFczd243RGhoMEhVTndnS0RNT2RLaEUwTnduNGdKdVBOV00zK3E1Rm0yN2ZiWm03dzZIZG9jODFRc3M4eGprRmJsYmswZUp2ek9Na0Z2cE4yWXFzNXhsYnI0Y3k3Rm1kMlQ3amcwZUlJVCtRU202eDVQM09Wc2hXNFBKVDRQc2tPUWJLdjIrUlZwUDdMcDhoRUZWZTY1RFJLZFdXcllROVY5LzhVeElWNGczSXFicmU5dktETGxPOFpUTG4yWTA5a0IwWUlpRzB6a1gyNTc2d1RmQVYyaHk5bzlXU0IxOXBaWmw0dHZhRmNYUGhBQ2J1N3hGTXVFWFpkSXJSc2c3amdXc1NaNlh1cTZqc3NvaVpyZFNwbFgrSzhQaUhBQTBOOEdlajFaZTRJSUEvSTdqRUxuMVFoc3piSEJwYWRTbDFoejFXZEVoNHd3cUZSNUR5TFVIUDlQU1pXeVZKTWtTRHZiZlY5a09nbi9LWDJaUE13bVQ2Q0RnUGRxa25HODNTTjlVRFFLSEIyZEVGZ0czRmZWNWRBYVlTRm9hdmc2YkQ0WWRRTnFEQmNSSjNPenNyTy8zY2VaYWJVQWQ1N0RLeHAvWEYvajFQSE1sTG1NNXBIWk5pc2NNZ0VMRVNDNjAxK1A3ZU9CSWhkcHdvZDA2bGVMWlEzNnh4L0pvSndodWlkdEZNbm9UV0F0RDRQTWdwcFpxMWpsUHhMd1RzMCticm83UTRNUGZ4NnF1QjdpK3BpKy9VWkJHYVU3L2k5a2N3SURBS3o2akVLZzR6TGxsdlJ1Q3crbjB1akdaTGVtRWl6RlBLU00wWllORkhxU3hKZTB5R0UycFprdEtDeVErb0taNW4weUYvSUgrczVkdDRFMmRiTnpNSU9nbllZQkwzS0FzNzRaemlHVVpnMWw4b1VRdkora1MyYzlSS3ZhUDBXRmJ2NllWMHlvL1NKRkxHc2dtK21mTVRBNXgyTStZZzdPWWlMaGVreCtHZE8wcTlXZVA5OG0vWGJOZ3VCR1ZIMFZlRkhzVml6N3JRU09CblBoODJMcjJ6Y3JXRzROVXhoMjRNRlpqR3I3MVhhYjFpRHZLMi9WMVZtZnZYWDArRGlYdz09Jmx0Oy9kaWFncmFtJmd0OyZsdDsvbXhmaWxlJmd0OyI+PGRlZnMvPjxnPjxyZWN0IHg9IjQiIHk9IjkwIiB3aWR0aD0iNDkwIiBoZWlnaHQ9IjkyIiBmaWxsPSIjZDVlOGQ0IiBzdHJva2U9Im5vbmUiIHBvaW50ZXItZXZlbnRzPSJub25lIi8+PHJlY3QgeD0iNCIgeT0iMTgyIiB3aWR0aD0iNDkwIiBoZWlnaHQ9Ijk0IiBmaWxsPSIjZmZlOGRlIiBzdHJva2U9Im5vbmUiIHBvaW50ZXItZXZlbnRzPSJub25lIi8+PHJlY3QgeD0iNCIgeT0iNCIgd2lkdGg9IjQ5MCIgaGVpZ2h0PSI4NiIgZmlsbD0iI2VmZjJmYyIgc3Ryb2tlPSJub25lIiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxyZWN0IHg9IjQiIHk9IjQiIHdpZHRoPSIxNDAiIGhlaWdodD0iMjcyIiBmaWxsLW9wYWNpdHk9IjAuOCIgZmlsbD0iI2Y1ZjVmNSIgc3Ryb2tlPSJub25lIiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiAxMXB4OyBtYXJnaW4tbGVmdDogMTE1cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogOHB4OyBmb250LWZhbWlseTogJnF1b3Q7VGltZXMgTmV3IFJvbWFuJnF1b3Q7OyBjb2xvcjogcmdiKDAsIDAsIDApOyBsaW5lLWhlaWdodDogMS4yOyBwb2ludGVyLWV2ZW50czogbm9uZTsgd2hpdGUtc3BhY2U6IG5vd3JhcDsiPnNvbGlkPC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjExNSIgeT0iMTMiIGZpbGw9InJnYmEoMCwgMCwgMCwgMSkiIGZvbnQtZmFtaWx5PSJUaW1lcyBOZXcgUm9tYW4iIGZvbnQtc2l6ZT0iOHB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5zb2xpZDwvdGV4dD48L3N3aXRjaD48L2c+PGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTAuNSAtMC41KXNjYWxlKDIpIj48c3dpdGNoPjxmb3JlaWduT2JqZWN0IHBvaW50ZXItZXZlbnRzPSJub25lIiB3aWR0aD0iMTAwJSIgaGVpZ2h0PSIxMDAlIiByZXF1aXJlZEZlYXR1cmVzPSJodHRwOi8vd3d3LnczLm9yZy9UUi9TVkcxMS9mZWF0dXJlI0V4dGVuc2liaWxpdHkiIHN0eWxlPSJvdmVyZmxvdzogdmlzaWJsZTsgdGV4dC1hbGlnbjogbGVmdDsiPjxkaXYgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGh0bWwiIHN0eWxlPSJkaXNwbGF5OiBmbGV4OyBhbGlnbi1pdGVtczogdW5zYWZlIGNlbnRlcjsganVzdGlmeS1jb250ZW50OiB1bnNhZmUgY2VudGVyOyB3aWR0aDogMXB4OyBoZWlnaHQ6IDFweDsgcGFkZGluZy10b3A6IDEwcHg7IG1hcmdpbi1sZWZ0OiAxOThweDsiPjxkaXYgZGF0YS1kcmF3aW8tY29sb3JzPSJjb2xvcjogcmdiYSgwLCAwLCAwLCAxKTsgIiBzdHlsZT0iYm94LXNpemluZzogYm9yZGVyLWJveDsgZm9udC1zaXplOiAwcHg7IHRleHQtYWxpZ246IGNlbnRlcjsiPjxkaXYgc3R5bGU9ImRpc3BsYXk6IGlubGluZS1ibG9jazsgZm9udC1zaXplOiA4cHg7IGZvbnQtZmFtaWx5OiAmcXVvdDtUaW1lcyBOZXcgUm9tYW4mcXVvdDs7IGNvbG9yOiByZ2IoMCwgMCwgMCk7IGxpbmUtaGVpZ2h0OiAxLjI7IHBvaW50ZXItZXZlbnRzOiBub25lOyB3aGl0ZS1zcGFjZTogbm93cmFwOyI+ZGFzaGVkPC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjE5OCIgeT0iMTIiIGZpbGw9InJnYmEoMCwgMCwgMCwgMSkiIGZvbnQtZmFtaWx5PSJUaW1lcyBOZXcgUm9tYW4iIGZvbnQtc2l6ZT0iOHB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5kYXNoZWQ8L3RleHQ+PC9zd2l0Y2g+PC9nPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiAzMnB4OyBtYXJnaW4tbGVmdDogMTE2cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTFweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3dyYXA7Ij5NYW5hZ2VtZW50PC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjExNiIgeT0iMzUiIGZpbGw9InJnYmEoMCwgMCwgMCwgMSkiIGZvbnQtZmFtaWx5PSJUaW1lcyBOZXcgUm9tYW4iIGZvbnQtc2l6ZT0iMTFweCIgdGV4dC1hbmNob3I9Im1pZGRsZSI+TWFuYWdlbWVudDwvdGV4dD48L3N3aXRjaD48L2c+PGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTAuNSAtMC41KXNjYWxlKDIpIj48c3dpdGNoPjxmb3JlaWduT2JqZWN0IHBvaW50ZXItZXZlbnRzPSJub25lIiB3aWR0aD0iMTAwJSIgaGVpZ2h0PSIxMDAlIiByZXF1aXJlZEZlYXR1cmVzPSJodHRwOi8vd3d3LnczLm9yZy9UUi9TVkcxMS9mZWF0dXJlI0V4dGVuc2liaWxpdHkiIHN0eWxlPSJvdmVyZmxvdzogdmlzaWJsZTsgdGV4dC1hbGlnbjogbGVmdDsiPjxkaXYgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGh0bWwiIHN0eWxlPSJkaXNwbGF5OiBmbGV4OyBhbGlnbi1pdGVtczogdW5zYWZlIGNlbnRlcjsganVzdGlmeS1jb250ZW50OiB1bnNhZmUgY2VudGVyOyB3aWR0aDogMXB4OyBoZWlnaHQ6IDFweDsgcGFkZGluZy10b3A6IDMycHg7IG1hcmdpbi1sZWZ0OiAxOThweDsiPjxkaXYgZGF0YS1kcmF3aW8tY29sb3JzPSJjb2xvcjogcmdiYSgwLCAwLCAwLCAxKTsgIiBzdHlsZT0iYm94LXNpemluZzogYm9yZGVyLWJveDsgZm9udC1zaXplOiAwcHg7IHRleHQtYWxpZ246IGNlbnRlcjsiPjxkaXYgc3R5bGU9ImRpc3BsYXk6IGlubGluZS1ibG9jazsgZm9udC1zaXplOiAxMXB4OyBmb250LWZhbWlseTogJnF1b3Q7VGltZXMgTmV3IFJvbWFuJnF1b3Q7OyBjb2xvcjogcmdiKDAsIDAsIDApOyBsaW5lLWhlaWdodDogMS4yOyBwb2ludGVyLWV2ZW50czogbm9uZTsgd2hpdGUtc3BhY2U6IG5vd3JhcDsiPkNvbXB1dGU8L2Rpdj48L2Rpdj48L2Rpdj48L2ZvcmVpZ25PYmplY3Q+PHRleHQgeD0iMTk4IiB5PSIzNSIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMXB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5Db21wdXRlPC90ZXh0Pjwvc3dpdGNoPjwvZz48ZWxsaXBzZSBjeD0iMjMyLjUiIGN5PSIxMjMiIHJ4PSI1MS41IiByeT0iMjciIGZpbGw9InJnYmEoMjU1LCAyNTUsIDI1NSwgMSkiIHN0cm9rZT0icmdiYSgwLCAwLCAwLCAxKSIgc3Ryb2tlLXdpZHRoPSIyIiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDUwcHg7IGhlaWdodDogMXB4OyBwYWRkaW5nLXRvcDogNjJweDsgbWFyZ2luLWxlZnQ6IDkycHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTJweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3JtYWw7IG92ZXJmbG93LXdyYXA6IG5vcm1hbDsiPjxzcGFuIHN0eWxlPSJmb250LXNpemU6IDhweCI+ZWxsaXBzZTwvc3Bhbj48L2Rpdj48L2Rpdj48L2Rpdj48L2ZvcmVpZ25PYmplY3Q+PHRleHQgeD0iMTE2IiB5PSI2NSIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMnB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5lbGxpcHNlPC90ZXh0Pjwvc3dpdGNoPjwvZz48ZyB0cmFuc2Zvcm09InRyYW5zbGF0ZSgtMC41IC0wLjUpc2NhbGUoMikiPjxzd2l0Y2g+PGZvcmVpZ25PYmplY3QgcG9pbnRlci1ldmVudHM9Im5vbmUiIHdpZHRoPSIxMDAlIiBoZWlnaHQ9IjEwMCUiIHJlcXVpcmVkRmVhdHVyZXM9Imh0dHA6Ly93d3cudzMub3JnL1RSL1NWRzExL2ZlYXR1cmUjRXh0ZW5zaWJpbGl0eSIgc3R5bGU9Im92ZXJmbG93OiB2aXNpYmxlOyB0ZXh0LWFsaWduOiBsZWZ0OyI+PGRpdiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMTk5OS94aHRtbCIgc3R5bGU9ImRpc3BsYXk6IGZsZXg7IGFsaWduLWl0ZW1zOiB1bnNhZmUgY2VudGVyOyBqdXN0aWZ5LWNvbnRlbnQ6IHVuc2FmZSBjZW50ZXI7IHdpZHRoOiAxcHg7IGhlaWdodDogMXB4OyBwYWRkaW5nLXRvcDogODVweDsgbWFyZ2luLWxlZnQ6IDExNHB4OyI+PGRpdiBkYXRhLWRyYXdpby1jb2xvcnM9ImNvbG9yOiByZ2JhKDAsIDAsIDAsIDEpOyAiIHN0eWxlPSJib3gtc2l6aW5nOiBib3JkZXItYm94OyBmb250LXNpemU6IDBweDsgdGV4dC1hbGlnbjogY2VudGVyOyI+PGRpdiBzdHlsZT0iZGlzcGxheTogaW5saW5lLWJsb2NrOyBmb250LXNpemU6IDExcHg7IGZvbnQtZmFtaWx5OiAmcXVvdDtUaW1lcyBOZXcgUm9tYW4mcXVvdDs7IGNvbG9yOiByZ2IoMCwgMCwgMCk7IGxpbmUtaGVpZ2h0OiAxLjI7IHBvaW50ZXItZXZlbnRzOiBub25lOyB3aGl0ZS1zcGFjZTogbm93cmFwOyI+UjwvZGl2PjwvZGl2PjwvZGl2PjwvZm9yZWlnbk9iamVjdD48dGV4dCB4PSIxMTQiIHk9Ijg4IiBmaWxsPSJyZ2JhKDAsIDAsIDAsIDEpIiBmb250LWZhbWlseT0iVGltZXMgTmV3IFJvbWFuIiBmb250LXNpemU9IjExcHgiIHRleHQtYW5jaG9yPSJtaWRkbGUiPlI8L3RleHQ+PC9zd2l0Y2g+PC9nPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiA4M3B4OyBtYXJnaW4tbGVmdDogMTk4cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTFweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3dyYXA7Ij5Db21tYW5kLWxpbmU8L2Rpdj48L2Rpdj48L2Rpdj48L2ZvcmVpZ25PYmplY3Q+PHRleHQgeD0iMTk4IiB5PSI4NiIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMXB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5Db21tYW5kLWxpbmU8L3RleHQ+PC9zd2l0Y2g+PC9nPjxyZWN0IHg9IjM0OSIgeT0iOTYiIHdpZHRoPSIxMDUiIGhlaWdodD0iNTAiIHJ4PSI3LjUiIHJ5PSI3LjUiIGZpbGw9InJnYmEoMjU1LCAyNTUsIDI1NSwgMSkiIHN0cm9rZT0icmdiYSgwLCAwLCAwLCAxKSIgc3Ryb2tlLXdpZHRoPSIyIiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDUxcHg7IGhlaWdodDogMXB4OyBwYWRkaW5nLXRvcDogNjFweDsgbWFyZ2luLWxlZnQ6IDE3NnB4OyI+PGRpdiBkYXRhLWRyYXdpby1jb2xvcnM9ImNvbG9yOiByZ2JhKDAsIDAsIDAsIDEpOyAiIHN0eWxlPSJib3gtc2l6aW5nOiBib3JkZXItYm94OyBmb250LXNpemU6IDBweDsgdGV4dC1hbGlnbjogY2VudGVyOyI+PGRpdiBzdHlsZT0iZGlzcGxheTogaW5saW5lLWJsb2NrOyBmb250LXNpemU6IDhweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3JtYWw7IG92ZXJmbG93LXdyYXA6IG5vcm1hbDsiPnJlY3RhbmdsZTwvZGl2PjwvZGl2PjwvZGl2PjwvZm9yZWlnbk9iamVjdD48dGV4dCB4PSIyMDEiIHk9IjYzIiBmaWxsPSJyZ2JhKDAsIDAsIDAsIDEpIiBmb250LWZhbWlseT0iVGltZXMgTmV3IFJvbWFuIiBmb250LXNpemU9IjhweCIgdGV4dC1hbmNob3I9Im1pZGRsZSI+cmVjdGFuZ2xlPC90ZXh0Pjwvc3dpdGNoPjwvZz48cGF0aCBkPSJNIDE4Mi41IDM4IEwgMjg3LjUgMzgiIGZpbGw9Im5vbmUiIHN0cm9rZT0icmdiYSgwLCAwLCAwLCAxKSIgc3Ryb2tlLXdpZHRoPSI2IiBzdHJva2UtbWl0ZXJsaW1pdD0iMTAiIHBvaW50ZXItZXZlbnRzPSJub25lIi8+PHBhdGggZD0iTSAzNTQgMzcuNjIgTCA0NTkgMzcuNjIiIGZpbGw9Im5vbmUiIHN0cm9rZT0icmdiYSgwLCAwLCAwLCAxKSIgc3Ryb2tlLXdpZHRoPSI2IiBzdHJva2UtbWl0ZXJsaW1pdD0iMTAiIHN0cm9rZS1kYXNoYXJyYXk9IjE4IDE4IiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiAxMjhweDsgbWFyZ2luLWxlZnQ6IDExNXB4OyI+PGRpdiBkYXRhLWRyYXdpby1jb2xvcnM9ImNvbG9yOiByZ2JhKDAsIDAsIDAsIDEpOyAiIHN0eWxlPSJib3gtc2l6aW5nOiBib3JkZXItYm94OyBmb250LXNpemU6IDBweDsgdGV4dC1hbGlnbjogY2VudGVyOyI+PGRpdiBzdHlsZT0iZGlzcGxheTogaW5saW5lLWJsb2NrOyBmb250LXNpemU6IDExcHg7IGZvbnQtZmFtaWx5OiAmcXVvdDtUaW1lcyBOZXcgUm9tYW4mcXVvdDs7IGNvbG9yOiByZ2IoMCwgMCwgMCk7IGxpbmUtaGVpZ2h0OiAxLjI7IHBvaW50ZXItZXZlbnRzOiBub25lOyB3aGl0ZS1zcGFjZTogbm93cmFwOyI+TWFuZGF0b3J5PC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjExNSIgeT0iMTMxIiBmaWxsPSJyZ2JhKDAsIDAsIDAsIDEpIiBmb250LWZhbWlseT0iVGltZXMgTmV3IFJvbWFuIiBmb250LXNpemU9IjExcHgiIHRleHQtYW5jaG9yPSJtaWRkbGUiPk1hbmRhdG9yeTwvdGV4dD48L3N3aXRjaD48L2c+PGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTAuNSAtMC41KXNjYWxlKDIpIj48c3dpdGNoPjxmb3JlaWduT2JqZWN0IHBvaW50ZXItZXZlbnRzPSJub25lIiB3aWR0aD0iMTAwJSIgaGVpZ2h0PSIxMDAlIiByZXF1aXJlZEZlYXR1cmVzPSJodHRwOi8vd3d3LnczLm9yZy9UUi9TVkcxMS9mZWF0dXJlI0V4dGVuc2liaWxpdHkiIHN0eWxlPSJvdmVyZmxvdzogdmlzaWJsZTsgdGV4dC1hbGlnbjogbGVmdDsiPjxkaXYgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGh0bWwiIHN0eWxlPSJkaXNwbGF5OiBmbGV4OyBhbGlnbi1pdGVtczogdW5zYWZlIGNlbnRlcjsganVzdGlmeS1jb250ZW50OiB1bnNhZmUgY2VudGVyOyB3aWR0aDogMXB4OyBoZWlnaHQ6IDFweDsgcGFkZGluZy10b3A6IDEyOHB4OyBtYXJnaW4tbGVmdDogMTk4cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTFweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3dyYXA7Ij5PcHRpb25hbDwvZGl2PjwvZGl2PjwvZGl2PjwvZm9yZWlnbk9iamVjdD48dGV4dCB4PSIxOTgiIHk9IjEzMSIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMXB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5PcHRpb25hbDwvdGV4dD48L3N3aXRjaD48L2c+PHJlY3QgeD0iMTg0IiB5PSIxOTAiIHdpZHRoPSI5NSIgaGVpZ2h0PSI0MCIgZmlsbD0iI2QzZDZlYiIgc3Ryb2tlPSJub25lIiBwb2ludGVyLWV2ZW50cz0ibm9uZSIvPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDQ2cHg7IGhlaWdodDogMXB4OyBwYWRkaW5nLXRvcDogMTA1cHg7IG1hcmdpbi1sZWZ0OiA5M3B4OyI+PGRpdiBkYXRhLWRyYXdpby1jb2xvcnM9ImNvbG9yOiByZ2JhKDAsIDAsIDAsIDEpOyAiIHN0eWxlPSJib3gtc2l6aW5nOiBib3JkZXItYm94OyBmb250LXNpemU6IDBweDsgdGV4dC1hbGlnbjogY2VudGVyOyI+PGRpdiBzdHlsZT0iZGlzcGxheTogaW5saW5lLWJsb2NrOyBmb250LXNpemU6IDEycHg7IGZvbnQtZmFtaWx5OiAmcXVvdDtUaW1lcyBOZXcgUm9tYW4mcXVvdDs7IGNvbG9yOiByZ2IoMCwgMCwgMCk7IGxpbmUtaGVpZ2h0OiAxLjI7IHBvaW50ZXItZXZlbnRzOiBub25lOyB3aGl0ZS1zcGFjZTogbm9ybWFsOyBvdmVyZmxvdy13cmFwOiBub3JtYWw7Ij48c3BhbiBzdHlsZT0iZm9udC1zaXplOiA4cHgiPmZpbGw8L3NwYW4+PC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjExNiIgeT0iMTA5IiBmaWxsPSJyZ2JhKDAsIDAsIDAsIDEpIiBmb250LWZhbWlseT0iVGltZXMgTmV3IFJvbWFuIiBmb250LXNpemU9IjEycHgiIHRleHQtYW5jaG9yPSJtaWRkbGUiPmZpbGw8L3RleHQ+PC9zd2l0Y2g+PC9nPjxyZWN0IHg9IjM0OSIgeT0iMTkwIiB3aWR0aD0iOTUiIGhlaWdodD0iNDAiIGZpbGw9IiNmZmZmZmYiIHN0cm9rZT0ibm9uZSIgcG9pbnRlci1ldmVudHM9Im5vbmUiLz48ZyB0cmFuc2Zvcm09InRyYW5zbGF0ZSgtMC41IC0wLjUpc2NhbGUoMikiPjxzd2l0Y2g+PGZvcmVpZ25PYmplY3QgcG9pbnRlci1ldmVudHM9Im5vbmUiIHdpZHRoPSIxMDAlIiBoZWlnaHQ9IjEwMCUiIHJlcXVpcmVkRmVhdHVyZXM9Imh0dHA6Ly93d3cudzMub3JnL1RSL1NWRzExL2ZlYXR1cmUjRXh0ZW5zaWJpbGl0eSIgc3R5bGU9Im92ZXJmbG93OiB2aXNpYmxlOyB0ZXh0LWFsaWduOiBsZWZ0OyI+PGRpdiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMTk5OS94aHRtbCIgc3R5bGU9ImRpc3BsYXk6IGZsZXg7IGFsaWduLWl0ZW1zOiB1bnNhZmUgY2VudGVyOyBqdXN0aWZ5LWNvbnRlbnQ6IHVuc2FmZSBjZW50ZXI7IHdpZHRoOiA0NnB4OyBoZWlnaHQ6IDFweDsgcGFkZGluZy10b3A6IDEwNXB4OyBtYXJnaW4tbGVmdDogMTc2cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTJweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IHdoaXRlLXNwYWNlOiBub3JtYWw7IG92ZXJmbG93LXdyYXA6IG5vcm1hbDsiPjxzcGFuIHN0eWxlPSJmb250LXNpemU6IDhweCI+bm8gZmlsbDwvc3Bhbj48L2Rpdj48L2Rpdj48L2Rpdj48L2ZvcmVpZ25PYmplY3Q+PHRleHQgeD0iMTk4IiB5PSIxMDkiIGZpbGw9InJnYmEoMCwgMCwgMCwgMSkiIGZvbnQtZmFtaWx5PSJUaW1lcyBOZXcgUm9tYW4iIGZvbnQtc2l6ZT0iMTJweCIgdGV4dC1hbmNob3I9Im1pZGRsZSI+bm8gZmlsbDwvdGV4dD48L3N3aXRjaD48L2c+PGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTAuNSAtMC41KXNjYWxlKDIpIj48c3dpdGNoPjxmb3JlaWduT2JqZWN0IHBvaW50ZXItZXZlbnRzPSJub25lIiB3aWR0aD0iMTAwJSIgaGVpZ2h0PSIxMDAlIiByZXF1aXJlZEZlYXR1cmVzPSJodHRwOi8vd3d3LnczLm9yZy9UUi9TVkcxMS9mZWF0dXJlI0V4dGVuc2liaWxpdHkiIHN0eWxlPSJvdmVyZmxvdzogdmlzaWJsZTsgdGV4dC1hbGlnbjogbGVmdDsiPjxkaXYgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGh0bWwiIHN0eWxlPSJkaXNwbGF5OiBmbGV4OyBhbGlnbi1pdGVtczogdW5zYWZlIGNlbnRlcjsganVzdGlmeS1jb250ZW50OiB1bnNhZmUgY2VudGVyOyB3aWR0aDogMXB4OyBoZWlnaHQ6IDFweDsgcGFkZGluZy10b3A6IDI0cHg7IG1hcmdpbi1sZWZ0OiAzNXB4OyI+PGRpdiBkYXRhLWRyYXdpby1jb2xvcnM9ImNvbG9yOiByZ2JhKDAsIDAsIDAsIDEpOyAiIHN0eWxlPSJib3gtc2l6aW5nOiBib3JkZXItYm94OyBmb250LXNpemU6IDBweDsgdGV4dC1hbGlnbjogY2VudGVyOyI+PGRpdiBzdHlsZT0iZGlzcGxheTogaW5saW5lLWJsb2NrOyBmb250LXNpemU6IDEwcHg7IGZvbnQtZmFtaWx5OiAmcXVvdDtUaW1lcyBOZXcgUm9tYW4mcXVvdDs7IGNvbG9yOiByZ2IoMCwgMCwgMCk7IGxpbmUtaGVpZ2h0OiAxLjI7IHBvaW50ZXItZXZlbnRzOiBub25lOyBmb250LXdlaWdodDogYm9sZDsgd2hpdGUtc3BhY2U6IG5vd3JhcDsiPlJ1bm5pbmcgPGJyIHN0eWxlPSJmb250LXNpemU6IDEwcHgiIC8+U2Vzc2lvbjwvZGl2PjwvZGl2PjwvZGl2PjwvZm9yZWlnbk9iamVjdD48dGV4dCB4PSIzNSIgeT0iMjciIGZpbGw9InJnYmEoMCwgMCwgMCwgMSkiIGZvbnQtZmFtaWx5PSJUaW1lcyBOZXcgUm9tYW4iIGZvbnQtc2l6ZT0iMTBweCIgdGV4dC1hbmNob3I9Im1pZGRsZSIgZm9udC13ZWlnaHQ9ImJvbGQiPlJ1bm5pbmcuLi48L3RleHQ+PC9zd2l0Y2g+PC9nPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiAxMTNweDsgbWFyZ2luLWxlZnQ6IDM1cHg7Ij48ZGl2IGRhdGEtZHJhd2lvLWNvbG9ycz0iY29sb3I6IHJnYmEoMCwgMCwgMCwgMSk7ICIgc3R5bGU9ImJveC1zaXppbmc6IGJvcmRlci1ib3g7IGZvbnQtc2l6ZTogMHB4OyB0ZXh0LWFsaWduOiBjZW50ZXI7Ij48ZGl2IHN0eWxlPSJkaXNwbGF5OiBpbmxpbmUtYmxvY2s7IGZvbnQtc2l6ZTogMTBweDsgZm9udC1mYW1pbHk6ICZxdW90O1RpbWVzIE5ldyBSb21hbiZxdW90OzsgY29sb3I6IHJnYigwLCAwLCAwKTsgbGluZS1oZWlnaHQ6IDEuMjsgcG9pbnRlci1ldmVudHM6IG5vbmU7IGZvbnQtd2VpZ2h0OiBib2xkOyB3aGl0ZS1zcGFjZTogbm93cmFwOyI+PGRpdiBzdHlsZT0iZm9udC1zaXplOiAxMHB4Ij48c3BhbiBzdHlsZT0iYmFja2dyb3VuZC1jb2xvcjogdHJhbnNwYXJlbnQgOyBmb250LXNpemU6IDEwcHgiPlJ1bm5pbmc8L3NwYW4+PC9kaXY+UmVxdWlyZW1lbnQ8L2Rpdj48L2Rpdj48L2Rpdj48L2ZvcmVpZ25PYmplY3Q+PHRleHQgeD0iMzUiIHk9IjExNiIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMHB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIiBmb250LXdlaWdodD0iYm9sZCI+UnVubmluZ1JlcXVpcmUuLi48L3RleHQ+PC9zd2l0Y2g+PC9nPjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0wLjUgLTAuNSlzY2FsZSgyKSI+PHN3aXRjaD48Zm9yZWlnbk9iamVjdCBwb2ludGVyLWV2ZW50cz0ibm9uZSIgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5IiBzdHlsZT0ib3ZlcmZsb3c6IHZpc2libGU7IHRleHQtYWxpZ246IGxlZnQ7Ij48ZGl2IHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hodG1sIiBzdHlsZT0iZGlzcGxheTogZmxleDsgYWxpZ24taXRlbXM6IHVuc2FmZSBjZW50ZXI7IGp1c3RpZnktY29udGVudDogdW5zYWZlIGNlbnRlcjsgd2lkdGg6IDFweDsgaGVpZ2h0OiAxcHg7IHBhZGRpbmctdG9wOiA2OHB4OyBtYXJnaW4tbGVmdDogMzVweDsiPjxkaXYgZGF0YS1kcmF3aW8tY29sb3JzPSJjb2xvcjogcmdiYSgwLCAwLCAwLCAxKTsgIiBzdHlsZT0iYm94LXNpemluZzogYm9yZGVyLWJveDsgZm9udC1zaXplOiAwcHg7IHRleHQtYWxpZ246IGNlbnRlcjsiPjxkaXYgc3R5bGU9ImRpc3BsYXk6IGlubGluZS1ibG9jazsgZm9udC1zaXplOiAxMHB4OyBmb250LWZhbWlseTogJnF1b3Q7VGltZXMgTmV3IFJvbWFuJnF1b3Q7OyBjb2xvcjogcmdiKDAsIDAsIDApOyBsaW5lLWhlaWdodDogMS4yOyBwb2ludGVyLWV2ZW50czogbm9uZTsgZm9udC13ZWlnaHQ6IGJvbGQ7IHdoaXRlLXNwYWNlOiBub3dyYXA7Ij5TdGVwIDxiciBzdHlsZT0iZm9udC1zaXplOiAxMHB4IiAvPkNsYXNzPC9kaXY+PC9kaXY+PC9kaXY+PC9mb3JlaWduT2JqZWN0Pjx0ZXh0IHg9IjM1IiB5PSI3MSIgZmlsbD0icmdiYSgwLCAwLCAwLCAxKSIgZm9udC1mYW1pbHk9IlRpbWVzIE5ldyBSb21hbiIgZm9udC1zaXplPSIxMHB4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIiBmb250LXdlaWdodD0iYm9sZCI+U3RlcC4uLjwvdGV4dD48L3N3aXRjaD48L2c+PC9nPjxzd2l0Y2g+PGcgcmVxdWlyZWRGZWF0dXJlcz0iaHR0cDovL3d3dy53My5vcmcvVFIvU1ZHMTEvZmVhdHVyZSNFeHRlbnNpYmlsaXR5Ii8+PGEgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoMCwtNSkiIHhsaW5rOmhyZWY9Imh0dHBzOi8vd3d3LmRpYWdyYW1zLm5ldC9kb2MvZmFxL3N2Zy1leHBvcnQtdGV4dC1wcm9ibGVtcyIgdGFyZ2V0PSJfYmxhbmsiPjx0ZXh0IHRleHQtYW5jaG9yPSJtaWRkbGUiIGZvbnQtc2l6ZT0iMTBweCIgeD0iNTAlIiB5PSIxMDAlIj5WaWV3ZXIgZG9lcyBub3Qgc3VwcG9ydCBmdWxsIFNWRyAxLjE8L3RleHQ+PC9hPjwvc3dpdGNoPjwvc3ZnPg=="},"evals":[],"jsHooks":[]}</script>

### Technical reports

*`systemPipeR`* compiles all the workflow execution logs in one central location,
making it easier to check any standard output (`stdout`) or standard error
(`stderr`) for any command-line tool used in a workflow.
Also, the information is appended to the workflow plot making it easy to click on
respective steps.

``` r
sal <- renderLogs(sal)
```

### Scientific reports

*`systemPipeR`* auto-generates scientific analysis reports in HTML format. These reports
compiling all the results of all workflow steps including plots and tables.

``` r
sal <- renderReport(sal)
```

## Workflow templates

Workflow templates are provided via the affiliated Bioconductor data package named `systemPipeRdata`, as well as on a dedicated GitHub repository.
Instances of these workflows can be created with a single command. The following give several examples.

### RNA-Seq WF template

Load the RNA-Seq sample workflow into your current working directory.

``` r
library(systemPipeRdata)
genWorkenvir(workflow = "rnaseq")
setwd("rnaseq")
```

#### Create the workflow

This template includes the most common analysis steps of `RNAseq` workflows. One can add, remove, modify
workflow steps by modifying the `sal` object.

``` r
sal <- SPRproject()
sal <- importWF(sal, file_path = "systemPipeRNAseq.Rmd", verbose = FALSE)
```

**Workflow steps:**

1.  Read preprocessing
      - Quality filtering (trimming)
      - FASTQ quality report
2.  Alignments: *`HISAT2`* (or any other RNA-Seq aligner)
3.  Alignment stats
4.  Read counting
5.  Sample-wise correlation analysis
6.  Analysis of differentially expressed genes (DEGs)
7.  GO term enrichment analysis
8.  Gene-wise clustering

#### Run workflow

``` r
sal <- runWF(sal)
```

Workflow visualization

``` r
plotWF(sal)
```

Report generation

``` r
sal <- renderReport(sal)
sal <- renderLogs(sal)
```

### ChIP-Seq WF template

Load the ChIP-Seq sample workflow into your current working directory.

``` r
library(systemPipeRdata)
genWorkenvir(workflow = "chipseq")
setwd("chipseq")
```

**Workflow steps:**

1.  Read preprocessing
      - Quality filtering (trimming)
      - FASTQ quality report
2.  Alignments: *`Bowtie2`* or *`rsubread`*
3.  Alignment stats
4.  Peak calling: *`MACS2`*
5.  Peak annotation with genomic context
6.  Differential binding analysis
7.  GO term enrichment analysis
8.  Motif analysis

#### Create the workflow

This template provides some common steps for a `ChIPseq` workflow. One can add, remove, modify
workflow steps by operating on the `sal` object.

``` r
sal <- SPRproject()
sal <- importWF(sal, file_path = "systemPipeChIPseq.Rmd", verbose = FALSE)
```

#### Run workflow

``` r
sal <- runWF(sal)
```

Workflow visualization

``` r
plotWF(sal)
```

Report generation

``` r
sal <- renderReport(sal)
sal <- renderLogs(sal)
```

### VAR-Seq WF template

Load the VAR-Seq sample workflow into your current working directory.

``` r
library(systemPipeRdata)
genWorkenvir(workflow = "varseq")
setwd("varseq")
```

**Workflow steps:**

1.  Read preprocessing
      - Quality filtering (trimming)
      - FASTQ quality report
2.  Alignments: *`gsnap`*, *`bwa`*
3.  Variant calling: *`VariantTools`*, *`GATK`*, *`BCFtools`*
4.  Variant filtering: *`VariantTools`* and *`VariantAnnotation`*
5.  Variant annotation: *`VariantAnnotation`*
6.  Combine results from many samples
7.  Summary statistics of samples

#### Create the workflow

This template provides some common steps for a `VARseq` workflow. One can add, remove, modify
workflow steps by operating on the `sal` object.

``` r
sal <- SPRproject()
sal <- importWF(sal, file_path = "systemPipeVARseq.Rmd", verbose = FALSE)
```

#### Run workflow

``` r
sal <- runWF(sal)
```

Workflow visualization

``` r
plotWF(sal)
```

Report generation

``` r
sal <- renderReport(sal)
sal <- renderLogs(sal)
```

### Ribo-Seq WF template

Load the Ribo-Seq sample workflow into your current working directory.

``` r
library(systemPipeRdata)
genWorkenvir(workflow = "riboseq")
setwd("riboseq")
```

**Workflow steps:**

1.  Read preprocessing
      - Adaptor trimming and quality filtering
      - FASTQ quality report
2.  Alignments: *`HISAT2`* (or any other RNA-Seq aligner)
3.  Alignment stats
4.  Compute read distribution across genomic features
5.  Adding custom features to workflow (e.g. uORFs)
6.  Genomic read coverage along transcripts
7.  Read counting
8.  Sample-wise correlation analysis
9.  Analysis of differentially expressed genes (DEGs)
10. GO term enrichment analysis
11. Gene-wise clustering
12. Differential ribosome binding (translational efficiency)

#### Create the workflow

This template provides some common steps for a `RIBOseq` workflow. One can add, remove, modify
workflow steps by operating on the `sal` object.

``` r
sal <- SPRproject()
sal <- importWF(sal, file_path = "systemPipeRIBOseq.Rmd", verbose = FALSE)
```

#### Run workflow

``` r
sal <- runWF(sal)
```

Workflow visualization

``` r
plotWF(sal, rstudio = TRUE)
```

Report generation

``` r
sal <- renderReport(sal)
sal <- renderLogs(sal)
```

## Version information

**Note:** the most recent version of this tutorial can be found <a href="http://www.bioconductor.org/packages/devel/bioc/vignettes/systemPipeR/inst/doc/systemPipeR.html">here</a>.

``` r
sessionInfo()
```

    ## R version 4.1.3 (2022-03-10)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Debian GNU/Linux 10 (buster)
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.8.0
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.8.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] magrittr_2.0.2              batchtools_0.9.15          
    ##  [3] ape_5.5                     ggplot2_3.3.5              
    ##  [5] systemPipeR_2.0.8           ShortRead_1.52.0           
    ##  [7] GenomicAlignments_1.30.0    SummarizedExperiment_1.24.0
    ##  [9] Biobase_2.54.0              MatrixGenerics_1.6.0       
    ## [11] matrixStats_0.61.0          BiocParallel_1.28.2        
    ## [13] Rsamtools_2.10.0            Biostrings_2.62.0          
    ## [15] XVector_0.34.0              GenomicRanges_1.46.1       
    ## [17] GenomeInfoDb_1.30.0         IRanges_2.28.0             
    ## [19] S4Vectors_0.32.3            BiocGenerics_0.40.0        
    ## [21] BiocStyle_2.22.0           
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] nlme_3.1-155           bitops_1.0-7           webshot_0.5.3         
    ##  [4] httr_1.4.2             RColorBrewer_1.1-2     progress_1.2.2        
    ##  [7] tools_4.1.3            backports_1.4.0        bslib_0.3.1           
    ## [10] utf8_1.2.2             R6_2.5.1               DBI_1.1.1             
    ## [13] colorspace_2.0-2       withr_2.4.3            tidyselect_1.1.1      
    ## [16] prettyunits_1.1.1      compiler_4.1.3         rvest_1.0.2           
    ## [19] cli_3.1.0              formatR_1.11           xml2_1.3.3            
    ## [22] DelayedArray_0.20.0    bookdown_0.24          sass_0.4.0            
    ## [25] scales_1.1.1           checkmate_2.0.0        rappdirs_0.3.3        
    ## [28] systemfonts_1.0.4      stringr_1.4.0          digest_0.6.29         
    ## [31] svglite_2.1.0          rmarkdown_2.13         jpeg_0.1-9            
    ## [34] pkgconfig_2.0.3        htmltools_0.5.2        fastmap_1.1.0         
    ## [37] htmlwidgets_1.5.4      rlang_1.0.2            rstudioapi_0.13       
    ## [40] jquerylib_0.1.4        generics_0.1.1         hwriter_1.3.2         
    ## [43] jsonlite_1.8.0         dplyr_1.0.7            RCurl_1.98-1.5        
    ## [46] kableExtra_1.3.4       GenomeInfoDbData_1.2.7 Matrix_1.4-0          
    ## [49] Rcpp_1.0.8.2           munsell_0.5.0          fansi_0.5.0           
    ## [52] lifecycle_1.0.1        stringi_1.7.6          yaml_2.3.5            
    ## [55] zlibbioc_1.40.0        grid_4.1.3             parallel_4.1.3        
    ## [58] crayon_1.4.2           lattice_0.20-45        hms_1.1.1             
    ## [61] knitr_1.37             pillar_1.6.4           base64url_1.4         
    ## [64] codetools_0.2-18       glue_1.6.2             evaluate_0.15         
    ## [67] blogdown_1.8.2         latticeExtra_0.6-29    data.table_1.14.2     
    ## [70] BiocManager_1.30.16    png_0.1-7              vctrs_0.3.8           
    ## [73] gtable_0.3.0           purrr_0.3.4            assertthat_0.2.1      
    ## [76] xfun_0.30              viridisLite_0.4.0      tibble_3.1.6          
    ## [79] ellipsis_0.3.2         brew_1.0-6

## Funding

This project is funded by NSF award [ABI-1661152](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1661152).

## References

<div id="refs" class="references">

<div id="ref-Amstutz2016-ka">

Amstutz, Peter, Michael R Crusoe, Nebojša Tijanić, Brad Chapman, John Chilton, Michael Heuer, Andrey Kartashov, et al. 2016. “Common Workflow Language, V1.0,” July. <https://doi.org/10.6084/m9.figshare.3115156.v2>.

</div>

<div id="ref-H_Backman2016-bt">

H Backman, Tyler W, and Thomas Girke. 2016. “systemPipeR: NGS workflow and report generation environment.” *BMC Bioinformatics* 17 (1): 388. <https://doi.org/10.1186/s12859-016-1241-0>.

</div>

<div id="ref-Howard2013-fq">

Howard, Brian E, Qiwen Hu, Ahmet Can Babaoglu, Manan Chandra, Monica Borghi, Xiaoping Tan, Luyan He, et al. 2013. “High-Throughput RNA Sequencing of Pseudomonas-Infected Arabidopsis Reveals Hidden Transcriptome Complexity and Novel Splice Variants.” *PLoS One* 8 (10): e74183. <https://doi.org/10.1371/journal.pone.0074183>.

</div>

<div id="ref-Kim2015-ve">

Kim, Daehwan, Ben Langmead, and Steven L Salzberg. 2015. “HISAT: A Fast Spliced Aligner with Low Memory Requirements.” *Nat. Methods* 12 (4): 357–60.

</div>

</div>
