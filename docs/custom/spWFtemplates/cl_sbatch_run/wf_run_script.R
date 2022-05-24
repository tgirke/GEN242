#!/usr/bin/env Rscript

## Set directory to wf root 

setwd("../")

## Run RNA-Seq workflow via sbatch

library(systemPipeR)
sal <- SPRproject()
sal <- importWF(sal, file_path = "systemPipeRNAseq.Rmd") # populates sal with WF steps defined in Rmd
resources <- list(conffile=".batchtools.conf.R",                                                                                                                                    
                  template="batchtools.slurm.tmpl",                                                                                                                                 
                  Njobs=18,
                  walltime=180,
                  ntasks=1,                                                                                                                                                         
                  ncpus=4,                                                                                                                                                          
                  memory=4096,
                  partition = "gen242"
                  )
sal <- addResources(sal, step = c("preprocessing", "trimming", "hisat2_mapping"), resources = resources)
sal <- runWF(sal) # runs entire workflow
sal <- renderReport(sal) # after workflow has completed render Rmd to HTML report 

