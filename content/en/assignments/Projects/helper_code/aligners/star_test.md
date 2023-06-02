#############
#### STAR ###
#############

### Read mapping with `STAR`


```r
library(systemPipeR)
```

```
## Loading required package: Rsamtools
```

```
## Loading required package: GenomeInfoDb
```

```
## Loading required package: BiocGenerics
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unsplit, which.max, which.min
```

```
## Loading required package: S4Vectors
```

```
## Loading required package: stats4
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:utils':
## 
##     findMatches
```

```
## The following objects are masked from 'package:base':
## 
##     expand.grid, I, unname
```

```
## Loading required package: IRanges
```

```
## Loading required package: GenomicRanges
```

```
## Loading required package: Biostrings
```

```
## Loading required package: XVector
```

```
## 
## Attaching package: 'Biostrings'
```

```
## The following object is masked from 'package:base':
## 
##     strsplit
```

```
## Loading required package: ShortRead
```

```
## Loading required package: BiocParallel
```

```
## Loading required package: GenomicAlignments
```

```
## Loading required package: SummarizedExperiment
```

```
## Loading required package: MatrixGenerics
```

```
## Loading required package: matrixStats
```

```
## 
## Attaching package: 'MatrixGenerics'
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
##     colWeightedMeans, colWeightedMedians, colWeightedSds,
##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
##     rowWeightedSds, rowWeightedVars
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## 
## Attaching package: 'Biobase'
```

```
## The following object is masked from 'package:MatrixGenerics':
## 
##     rowMedians
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     anyMissing, rowMedians
```

```r
# sal_test <- SPRproject(logs.dir= ".SPRproject_test") # use this line when .SPRproject_test doesn't exist yet
sal_test <- SPRproject(overwrite = TRUE, logs.dir= ".SPRproject_test")
```

```
## Creating directory:  /home/tgirke/tmp/GEN242/content/en/assignments/Projects/helper_code/aligners/data 
## Creating directory:  /home/tgirke/tmp/GEN242/content/en/assignments/Projects/helper_code/aligners/results 
## Creating directory '/home/tgirke/tmp/GEN242/content/en/assignments/Projects/helper_code/aligners/.SPRproject_test'
## Creating file '/home/tgirke/tmp/GEN242/content/en/assignments/Projects/helper_code/aligners/.SPRproject_test/SYSargsList.yml'
```


```r
appendStep(sal_test) <- LineWise(code = {
                library(systemPipeR)
                }, step_name = "load_SPR")
```

## Read preprocessing

### Preprocessing with `preprocessReads` function


```r
appendStep(sal_test) <- SYSargsList(
    step_name = "preprocessing",
    targets = "targetsPE.txt", dir = TRUE,
    wf_file = "preprocessReads/preprocessReads-pe.cwl",
    input_file = "preprocessReads/preprocessReads-pe.yml",
    dir_path = system.file("extdata/cwl", package = "systemPipeR"),
    inputvars = c(
        FileName1 = "_FASTQ_PATH1_",
        FileName2 = "_FASTQ_PATH2_",
        SampleName = "_SampleName_"
    ),
    dependency = c("load_SPR"))
```

## Alignments with `STAR`

### `STAR` Indexing


```r
appendStep(sal_test) <- SYSargsList(
    step_name = "star_index", 
    dir = FALSE, 
    targets=NULL, 
    wf_file = "star/star-index.cwl", 
    input_file="star/star-index.yml",
    dir_path="param/cwl", 
    dependency = "load_SPR"
)
```

### `STAR` mapping


```r
appendStep(sal_test) <- SYSargsList(
    step_name = "star_mapping",
    dir = TRUE, 
    targets ="preprocessing", 
    wf_file = "star-mapping-pe.cwl",
    input_file = "star-mapping-pe.yml",
    dir_path = "param/star_test",
    inputvars = c(preprocessReads_1 = "_FASTQ_PATH1_", preprocessReads_2 = "_FASTQ_PATH2_", 
                  SampleName = "_SampleName_"),
    rm_targets_col = c("FileName1", "FileName2"), 
    dependency = c("preprocessing", "star_index")
)
```


```r
## Return command-line calls for STAR
cmdlist(sal_test, step="star_mapping", targets=1)

## BAM outpaths required for read counting below
outpaths <- getColumn(sal_test, step = "star_mapping", "outfiles", column = "Aligned_toTranscriptome_out_bam")
file.exists(outpaths) # Will not return TRUE until STAR completed sucessfully
```


```r
## To run sal_test stepwise, make sure you have constructed your 
## sal_test object step-by-step starting from an empty sal_test
## as shown above under chunk: intialize_sal_for_testing 
sal_test <- runWF(sal_test, steps=c(1)) # increment step number one by one just for checking
sal_test
outpaths <- getColumn(sal_test, step = "star_mapping", "outfiles", column = "Aligned_toTranscriptome_out_bam")
outpaths
file.exists(outpaths) # Will not return TRUE until STAR completed sucessfully
```


```r
## The following can be used for setting up things initial testing
starPE <- loadWorkflow(targets = "targetsPE.txt", wf_file = "star-mapping-pe.cwl", 
                       input_file = "star-mapping-pe.yml", dir_path = "./param/star_test")
starPE <- renderWF(starPE, inputvars = c(FileName1 = "_FASTQ_PATH1_", FileName2 = "_FASTQ_PATH2_", 
                                         SampleName = "_SampleName_"))
cmdlist(starPE)
runCommandline(starPE, make_bam = FALSE)
```
