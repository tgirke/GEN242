---
title: Data Management for Course Projects
linkTitle: "Project Data"
description: >
type: docs
weight: 410
---

<br></br>


## Big data space on HPCC

All larger data sets of the coure projects will be organized in a big data space under
`/bigdata/gen242/<user_name>`. Within this space, each student will work in a subdirectory named after their project:

+ `/bigdata/gen242/<user_name>/<github_user_name>_project`

## Project GitHub repositories 

Students will work on their course projects within GitHub repositories, one for each student.
These project repositories are private and have been shared with each student.
To populate a course project with an initial project workflow, please follow the instructions
given in the following section.

## Generate workflow environment with real project data

1. Log in to the HPCC cluster and set your working directory to `bigdata` or (`/bigdata/gen242/<user_name>`)
2. Clone the GitHub repository for your project with `git clone ...` (URLs listed [here](https://bit.ly/3tJ3KuZ)) and then `cd` into this directory. As mentioned above, the project GitHub repos follow this naming convention: `<github_user_name>_project`.
2. Generate the workflow environment for your project on the HPCC cluster with `genWorkenvir` from `systemPipeRdata`.
3. Next, `cd` into the directory of your workflow, delete its default `data` and `results` directories, and then substitute them with empty directories outside of your project GitHub repos as follows:
   ```sh 
   mkdir ../../<workflow>_data
   mkdir ../../<workflow>_results
   ```
4. Within your workflow directory create symbolic links pointing to the new directories created in the previous step. For instance, the projects using the RNA-Seq workflow should create the symbolic links for their `data` and `results` directories like this:
   ```sh 
   ln -s /bigdata/gen242/<user_name>/rnaseq_data data
   ln -s /bigdata/gen242/<user_name>/rnaseq_results results
   ```
5. Add the workflow directory to the GitHub repository of your project with `git add -A` and then run `commit` and `push` as outlined in the GitHub instructions of this course [here](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/github/github/#github-basics-from-command-line). After this check whether the workflow directory and its content shows up on your project's online repos on GitHub. Very important: make sure that the `data` and `results` are empty at this point. If not investigate why and fix the problem in the corresponding step above.  
6. Download the FASTQ files of your project with `getSRAfastq` (see below) to the `data` directory of your project. 
7. Generate a proper `targets` file for your project where the first column(s) point(s) to the downloaded FASTQ files. In addition, provide sample names matching the experimental design (columns: `SampleNames` and `Factor`). More details about the structure of targets files are provided [here](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/systempiper/systempiper/#structure-of-targets-file). Ready to use targets files for both the RNA-Seq and ChIP-Seq project can be downloaded as tab separated (TSV) files from [here](https://github.com/tgirke/GEN242/tree/main/content/en/assignments/Projects/targets_files). Alternatively, one can download the corresponding Google Sheets with the `read_sheet` function from the `googlesheets4` package ([RNA-Seq GSheet](https://bit.ly/2QH19Ry) and [ChIP-Seq GSheet](https://bit.ly/2QFjTAV)). 
8. Inspect and adjust the `.param` files you will be using. For instance, make sure the software modules you are loading and the path to the reference genome are correct. 
9. Every time you start working on your project you `cd` into the directory of the repository and then run `git pull` to get the latest change. When you are done, you commit and push your changes back to GitHub with `git commit -am "some edits"; git push -u origin main`.

## Download of project data

After logging in to one of the computer nodes via `srun`, open R from within the GitHub respository of your project and then run the following code section, but only those that apply to your project.

### FASTQ files from SRA

#### Choose FASTQ data for your project

+ The FASTQ files for the ChIP-Seq project are from SRA study [SRP002174](http://www.ncbi.nlm.nih.gov/sra?term=SRP002174) ([Kaufman et al. 2010](http://www.ncbi.nlm.nih.gov/pubmed/20360106))
```r
sraidv <- paste("SRR0388", 45:51, sep="") 
```

+ The FASTQ files for the RNA-Seq project are from SRA study [SRP010938](http://www.ncbi.nlm.nih.gov/sra?term=SRP010938) ([Howard et al. 2013](http://www.ncbi.nlm.nih.gov/pubmed/24098335))
```r
sraidv <- paste("SRR4460", 27:44, sep="")
```

#### Load libraries and modules

```r
library(systemPipeR)                                                                                                                                                                
moduleload("sratoolkit/3.0.0")                                                                                                                                                      
system("vdb-config --prefetch-to-cwd") # sets download default to current directory                                                                                                 
# system('prefetch --help') # helps to speed up fastq-dump                                                                                                                          
# system('fastq-dump --help') # below uses this one for backwards compatibility                                                                                                     
# system('fasterq-dump --help') # faster than fastq-dump
```

#### Define download function
The following function downloads and extracts the FASTQ files for each project from SRA.
Internally, it uses the `prefetch` and `fastq-dump` utilities of the SRA Toolkit from NCBI.
The faster `fasterq-dump` alternative (see comment line below) is not used here for historical reasons.

```r
getSRAfastq <- function(sraid, threads=1) {                                                                                                                                         
    system(paste("prefetch", sraid)) # makes download faster                                                                                                                        
    system(paste("fastq-dump --split-files --gzip", sraid)) # gzip option makes it slower but saves storage space                                                                   
    # system(paste("fasterq-dump --threads 4 --split-files --progress ", sraid, "--outdir .")) # Faster alternative to fastq-dump                                                   
    unlink(x=sraid, recursive = TRUE, force = TRUE) # deletes sra download directory                                                                                                
}    
```

#### Run download

Note the following performs the download in serialized mode for the chosen data set and saves the extracted FASTQ files to 
the current working directory.
```r
mydir <- getwd(); setwd("data")
for(i in sraidv) getSRAfastq(sraid=i)
setwd(mydir)
```

Alternatively, the download can be performed in parallelized mode with `BiocParallel`. Please run this version only on one of the compute nodes.
```r
mydir <- getwd(); setwd("data")
# bplapply(sraidv, getSRAfastq, BPPARAM = MulticoreParam(workers=4))
setwd(mydir)
```

### Download reference genome and annotation

The following `downloadRefs` function downloads the _Arabidopsis thaliana_ genome sequence and GFF file from the [TAIR FTP site](ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/). 
It also assigns consistent chromosome identifiers to make them the same among both the genome sequence and the GFF file. This is
important for many analysis routines such as the read counting in the RNA-Seq workflow.  

```r
downloadRefs <- function(rerun=FALSE) {
    if(rerun==TRUE) {
        library(Biostrings)
        download.file("https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas", "./data/tair10.fasta")
        dna <- readDNAStringSet("./data/tair10.fasta")
        names(dna) <- paste(rep("Chr", 7), c(1:5, "M", "C"), sep="") # Fixes chromomse ids
        writeXStringSet(dna, "./data/tair10.fasta")
        download.file("https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff", "./data/tair10.gff")
        download.file("https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_functional_descriptions", "./data/tair10_functional_descriptions")
        }
}
```

After importing/sourcing the above function, execute it as follows:
```r
downloadRefs(rerun=TRUE) 
```

## Workflow Rmd file

To run the actual data analysis workflows, the RNA-Seq project can use the `systemPipeRNAseq.Rmd` file obtained from the `genWorkenvir(workflow='rnaseq')` call directly. However,
the ChIP-Seq group should use the Rmd linked on top right corner of the ChIP-Seq tutorial [here](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/spchipseq/spchipseq/) 
and then name the downloaded file `systemPipeChIPseq.Rmd`. To simplify this, the ChIP-Seq group members can run from the command-line within their `chipseq` workflow directory the 
following download command. 

```r
wget https://raw.githubusercontent.com/tgirke/GEN242/main/static/custom/spWFtemplates/spchipseq.Rmd -O systemPipeChIPseq.Rmd
```

This will assign the proper file name and overwrite the preloaded version of this file that has the same name. 

## Recommendations for running workflows

### Run instructions

The following provides recommendations and additional options to consider for
running and modifying workflows. This also includes parallelization settings
for the specific data used by the class projects. Note, additional details can
be found in this and other sections of the workflow introduction tutorial
[here](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/systempiper/systempiper/#loading-workflows-from-an-r-markdown). 

```r
library(systemPipeR)                                                                                                                                                                
sal <- SPRproject() # when running a WF for first time                                                                                                                                      
sal                                                                                                                                                                                 
sal <- importWF(sal, file_path = "systemPipeRNAseq.Rmd") # populates sal with WF steps defined in Rmd                                                                                                                      
sal
sal <- SPRproject(resume=TRUE) # when restarting a WF, skip above steps and resume WF with this command                                                                                                                                               
getRversion() # should be 4.1.2 or 4.2.0. Note, R version can be changed with `module load ...`                                                                                                                                                     
system("hostname") # should return number of a compute node; if not close Nvim-R session, log in to a compute node with srun and then restart Nvim-R session                                                                                                                                                                     
# sal <- runWF(sal) # runs WF serialized. Not recommended since this will take much longer than parallel mode introduced below by taking advantage of resource allocation
resources <- list(conffile=".batchtools.conf.R",                                                                                                                                    
                  template="batchtools.slurm.tmpl",                                                                                                                                 
                  Njobs=18, # chipseq should use here number of fastq files (7 or 8)                                                                                                                                                        
                  walltime=180, ## minutes                                                                                                                                          
                  ntasks=1,                                                                                                                                                         
                  ncpus=4,                                                                                                                                                          
                  memory=4096, ## Mb                                                                                                                                                
                  partition = "gen242"                                                                                                                                              
                  )                                                                                                                                                                 
## For RNA-Seq project use:
sal <- addResources(sal, step = c("preprocessing", "trimming", "hisat2_mapping"), resources = resources) # parallelizes time consuming computations assigned to `step` argument                                                                           
## For ChIP-Seq project use:
sal <- addResources(sal, step = c("preprocessing", "bowtie2_alignment"), resources = resources)
sal <- runWF(sal) # runs entire workflow; specific steps can be executed by assigning their corresponding position numbers within the workflow to the `steps` argument (see ?runWF)                                                                                                                                                               
sal <- renderReport(sal) # after workflow has completed render Rmd to HTML report (default name is SPR_Report.html) and view it via web browser which requires symbolic link in your ~/.html folder. 
rmarkdown::render("systemPipeRNAseq.Rmd", clean=TRUE, output_format="BiocStyle::html_document") # Alternative approach for rendering report from Rmd file instead of sal object
```

### Modify a workflow

If needed one can modify existing workflow steps in a pre-populated `SYSargsList` object, and potentially already executed WF, with the `replaceStep(sal) <-` replacement function. 
The following gives an example where step number 3 in a `SYSargsList` (sal) object will be updated with modified or new code. Note, this is a generalized example where the user
needs to insert the code lines and also adjust the values assigned to the arguments: `step_name` and `dependency`. 

```r
replaceStep(sal, step=3) <- LineWise(                                                                                                                                                        
    code = {                                                                                                                                                                        
        << my modified code lines >>
        },                                                                                                                                                                          
    step_name = << "my_step_name" >>,                                                                                                                                                        
    dependency = << "my_dependency" >>)
```

Subsequently, one can rerun the corresponding step (here 3) as follows: 

```r
runWF(sal, steps=3)
```

Note, any step in a workflow can only be run in isolation if its expected input exists (see `dependency`).

### Adding steps to a workflow

New steps can be added to the Rmd file of a workflow by inserting new R Markdown code chunks starting and ending with the usual `appendStep<-` syntax and then creating a new 
`SYSargsList` instance with `importWF` that contain the new step(s). To add steps to a pre-populated `SYSargsList` object, one can use the `after` argument of the `appendStep<-` 
function. The following example will add a new step after position 3 to the corresponding `sal` object. This can be useful if a longer workflow has already been completed and 
one only wants to make some refinements without re-running the entire workflow.

```r
appendStep(sal, after=3) <- << my_step_code >>
```

