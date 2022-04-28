---
title: "HW7 - RNA-Seq Analysis"
linkTitle: "HW6: RNA-Seq Analysis"
description: >
type: docs
weight: 307
---

<br></br>

<div style="text-align: right"> 
Source code downloads: &nbsp; <a href="https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/assignments/Homework/HW07/HW07.R" target="_blank">[ .R ]</a>
</div>

## A. Unstranded and strand-specific read counting 

- __Task 1__: Rerun the [RNA-Seq
  workflow](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/sprnaseq/sprnaseq/)
  with the toy data sets up to the read quantification step
  [here](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/sprnaseq/sprnaseq/#read-quantification). Note,
  the toy data set gets automatically loaded when intializing a workflow environment (directory structure) with the `genWorkenvir` 
  function (see tutorial [here](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/sprnaseq/sprnaseq/#experimental-design)). 
  In the read quantification step with `summarizeOverlaps` generate count
  tables for exons by genes (`eByg`) of the following three strand modes:

   1. Unstranded 
   2. Strand-specific for positive (sense) strand
   3. Strand-specific for negative (antisense) strand
   
   The solution for generating the unstranded read counts is given below. Note,
   the upstream steps of the RNA-Seq workflow only need to be rerun to generate
   the proper inputs for the read counting. Thus, they are not required to be
   included in the homework results (see `HW7.R` below).


```r
summarizeOverlaps(eByg, bfl, mode="Union", 
                    ignore.strand=TRUE, 
                    # preprocess.reads=invertStrand,
                    inter.feature=FALSE, 
                    singleEnd=FALSE)
```

Before attempting to solve this homework task please read the vignette
_Counting reads with `summarizeOverlaps`_
([here](http://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html))
from the `GenomicAlignments` package that defines the `summarizeOverlap`
function. In addition, the help file for `?summarizeOverlaps` provides useful information.

- __Task 2__: Provide R code that demonstrates that the two strand-specific count tables sum up to the values of the unstranded count table. 

- __Task 3__: Explain the utility (biological relevance) of the different strand counting modes used under Task 1. Include your explanation as comment text in your homework script (see `HW7.R` below). 

Note, for Tasks 1-3 only the code and/or text needs to be included in the homework submission (no data/result files). 

## B. Read counting for different feature types
- __Task 4__: Compute strand-specific count tables for the positive (sense) strand of the following feature types. The help files of `?exonsBy` and `?transcripts` provide useful information for solving these tasks. 

   1. Genes
   2. Exons
   3. Exons by genes 
   4. Introns by transcripts
   5. 5'-UTRs by transcripts

Note, for Tasks 4 only include the code and/or text in your homework submission (no data/result files). For details see below.

## C. DEG analysis

- __Task 5__: Perform the DEG analysis with `edgeR` as outlined under section 6 of the RNA-Seq workflow [here](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/sprnaseq/sprnaseq/#run-edger). 
Use in one case for the DEG analysis the unstranded count table as input (from Task 1.1) and in another the sense strand count table (from Task 1.2). 
Compare the DEG result of the two methods in two separate 4-way Venn diagrams for the same sample comparisons used in the workflow example 
[here](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/sprnaseq/sprnaseq/#venn-diagrams-of-deg-sets).

   1. 4-way Venn diagram for unstranded count table
   2. 4-way Venn diagram for sense strand count table

Note, for Tasks 5 include both the code and the resulting images in your homework submission. 

## Homework submission

Please submit the homework results in one well structured and annotated R
script to your private GitHub repository under `Homework/HW7/HW7.R`. Instead 
of an R script the homework can be submitted in form of an R Markdown (*Rmd) file.


### Due date

This homework is due on Thu, May 5th at 6:00 PM.

## Homework Solutions

See [here](https://raw.githubusercontent.com/tgirke/GEN242/main/static/custom/hw_solutions/hw7_solution.R).

