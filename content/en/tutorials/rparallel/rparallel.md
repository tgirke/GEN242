---
title: "Parallel Evaluations in R"
author: Thomas Girke
date: "Last update: 28 May, 2024" 
output:
  html_document:
    toc: true
    toc_float:
        collapsed: true
        smooth_scroll: true
    toc_depth: 3
    fig_caption: yes
    code_folding: show
    number_sections: true

fontsize: 14pt
bibliography: bibtex.bib
weight: 5
type: docs
---

<!--
- Compile from command-line
Rscript -e "rmarkdown::render('rparallel.Rmd', c('html_document'), clean=F); knitr::knit('rparallel.Rmd', tangle=TRUE)"
-->

<div style="text-align: right">

Source code downloads:    
\[ [Slides](https://girke.bioinformatics.ucr.edu/GEN242/slides/slides_12/) \]    
\[ [.Rmd](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/rparallel/rparallel.Rmd) \]    
\[ [.R](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/rparallel/rparallel.R) \]

</div>

## Overview

  - A general introduction to this topic is in the [Linux and HPCC Cluster](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/linux/linux/#queuing-system-slurm) manual of the GEN242 site.
  - R provides a large number of packages for parallel evaluations on multi-core, multi-socket and multi-node systems. The latter are usually referred to as computer clusters.
  - MPI is also supported
  - For an overview of parallelization packages available for R see [here](https://cran.r-project.org/web/views/HighPerformanceComputing.html)
  - One of the most comprehensive parallel computing environments for R is
    [`batchtools`](https://mllg.github.io/batchtools/articles/batchtools.html#migration). Older versions of this package were released under the name `BatchJobs` (Bischl et al. 2015).
  - `batchtools` supports both multi-core and multi-node computations with and without schedulers. By making use of
    cluster template files, most schedulers and queueing systems are supported (*e.g.* Torque, Sun Grid Engine, Slurm).
  - The `BiocParallel` package (see [here](https://bioconductor.org/packages/release/bioc/html/BiocParallel.html))
    provides similar functionalities as `batchtools`, but is tailored to use Bioconductor objects.

## Reminder: Traditional Job Submission for R

This topic is covered in more detail in other tutorials. The following only provides a very brief overview of this submission method.

**1.** Create Slurm submission script, here called [script\_name.sh](https://raw.githubusercontent.com/tgirke/GEN242/main/static/custom/slides/R_for_HPC/demo_files/script_name.sh) with:

``` bash
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:15:00 # 1 day and 15 minutes
#SBATCH --mail-user=useremail@address.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="some_test"
#SBATCH --partition="gen242" # Choose alternative partitions from: intel, batch, highmem, gpu, short
#SBATCH --account="gen242" # Same as above

Rscript my_script.R
```

**2.** Submit R script called [my\_script.R](https://raw.githubusercontent.com/tgirke/GEN242/main/static/custom/slides/R_for_HPC/demo_files/my_script.R) by above Slurm script with:

``` bash
sbatch script_name.sh
```

## Parallel Evaluations on Clusters with `batchtools`

  - The following introduces the usage of `batchtools` for a computer cluster
    using SLURM as scheduler (workload manager). SLURM is the scheduler used by
    the HPCC at UCR.
  - Similar instructions are provided in HPCC’s manual section covering
    `batchtools`
    [here](https://hpcc.ucr.edu/manuals_linux-cluster_parallelR.html)
  - To simplify the evaluation of the R code on the following slides, the
    corresponding text version is available for download from
    [here](https://raw.githubusercontent.com/tgirke/GEN242/main/static/custom/slides/R_for_HPC/demo_files/R_for_HPC_demo.R).

## Hands-on Demo of `batchtools`

### Set up working directory for SLURM

First login to your cluster account, open R and execute the following lines. This will
create a test directory (here `mytestdir`), redirect R into this directory and then download
the required files:

  - [`slurm.tmpl`](https://github.com/tgirke/GEN242/blob/main/content/en/tutorials/rparallel/demo_files/slurm.tmpl)
  - [`.batchtools.conf.R`](https://github.com/tgirke/GEN242/blob/main/content/en/tutorials/rparallel/demo_files/.batchtools.conf.R)

<!-- end list -->

``` r
dir.create("mytestdir")
setwd("mytestdir")
download.file("https://bit.ly/3Oh9dRO", "slurm.tmpl")
download.file("https://bit.ly/3KPBwou", ".batchtools.conf.R") 
```

### Load package and define some custom function

The following code defines a test function (here `myFct`) that will be run on the cluster for demonstration
purposes.

The test function (`myFct`) subsets the `iris` data frame by rows, and appends the host name and R version of each
node where the function was executed. The R version to be used on each node can be
specified in the `slurm.tmpl` file (under `module load`).

``` r
library('RenvModule')
module('load','slurm') # Loads slurm among other modules
library(batchtools)
myFct <- function(x) {
    Sys.sleep(10) # to see job in queue, pause for 10 sec
	result <- cbind(iris[x, 1:4,],
                    Node=system("hostname", intern=TRUE),
	                Rversion=paste(R.Version()[6:7], collapse="."))
	return(result)
    }
```

### Submit jobs from R to cluster

The following creates a `batchtools` registry, defines the number of jobs and resource requests, and then submits the jobs to the cluster
via SLURM.

``` r
reg <- makeRegistry(file.dir="myregdir", conf.file=".batchtools.conf.R")
Njobs <- 1:4 # Define number of jobs (here 4)
ids <- batchMap(fun=myFct, x=Njobs) 
done <- submitJobs(ids, reg=reg, resources=list(partition="short", walltime=120, ntasks=1, ncpus=1, memory=1024))
waitForJobs() # Wait until jobs are completed
```

### Summarize job status

After the jobs are completed one can inspect their status as follows.

``` r
getStatus() # Summarize job status
showLog(Njobs[1])
# killJobs(Njobs) # # Possible from within R or outside with scancel
```

### Access/assemble results

The results are stored as `.rds` files in the registry directory (here `myregdir`). One
can access them manually via `readRDS` or use various convenience utilities provided
by the `batchtools` package.

``` r
readRDS("myregdir/results/1.rds") # reads from rds file first result chunk
loadResult(1) 
lapply(Njobs, loadResult)
reduceResults(rbind) # Assemble result chunks in single data.frame
do.call("rbind", lapply(Njobs, loadResult))
```

### Remove registry directory from file system

By default existing registries will not be overwritten. If required one can explicitly
clean and delete them with the following functions.

``` r
clearRegistry() # Clear registry in R session
removeRegistry(wait=0, reg=reg) # Delete registry directory
# unlink("myregdir", recursive=TRUE) # Same as previous line
```

### Load registry into R

Loading a registry can be useful when accessing the results at a later state or
after moving them to a local system.

``` r
from_file <- loadRegistry("myregdir", conf.file=".batchtools.conf.R")
reduceResults(rbind)
```

## Conclusions

### Advantages of `batchtools`

  - many parallelization methods multiple cores, and across both multiple CPU sockets and nodes
  - most schedulers supported
  - takes full advantage of a cluster
  - robust job management by organizing results in registry file-based database
  - simplifies submission, monitoring and restart of jobs
  - well supported and maintained package

## Session Info

``` r
sessionInfo()
```

    ## R version 4.4.0 (2024-04-24)
    ## Platform: x86_64-pc-linux-gnu
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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] digest_0.6.35     R6_2.5.1          bookdown_0.39     fastmap_1.1.1    
    ##  [5] xfun_0.43         blogdown_1.19     cachem_1.0.8      knitr_1.46       
    ##  [9] htmltools_0.5.8.1 rmarkdown_2.26    lifecycle_1.0.4   cli_3.6.2        
    ## [13] sass_0.4.9        jquerylib_0.1.4   compiler_4.4.0    tools_4.4.0      
    ## [17] evaluate_0.23     bslib_0.7.0       yaml_2.3.8        jsonlite_1.8.8   
    ## [21] rlang_1.1.3

## References

<div id="refs" class="references hanging-indent">

<div id="ref-Bischl2015-rf">

Bischl, Bernd, Michel Lang, Olaf Mersmann, Jörg Rahnenführer, and Claus Weihs. 2015. “BatchJobs and BatchExperiments: Abstraction Mechanisms for Using R in Batch Environments.” *Journal of Statistical Software*. <http://www.jstatsoft.org/v64/i11/>.

</div>

</div>
