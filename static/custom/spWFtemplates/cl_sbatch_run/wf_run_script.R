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
                  partition = "girkelab"
                  )
sal <- addResources(sal, step = c("preprocessing", "trimming", "hisat2_mapping"), resources = resources)
sal <- runWF(sal) # runs entire workflow; specific steps can be executed by assigning their corresponding position numbers within the workflow to the `steps` argument (see ?runWF)                                                                                                                                                               
sal <- renderReport(sal) # after workflow has completed render Rmd to HTML report (default name is SPR_Report.html) and view it via web browser which requires symbolic link in your ~/.html folder. 

