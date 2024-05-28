## ----setenvir, eval=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## dir.create("mytestdir")
## setwd("mytestdir")
## download.file("https://bit.ly/3Oh9dRO", "slurm.tmpl")
## download.file("https://bit.ly/3KPBwou", ".batchtools.conf.R")


## ----custom_fct1, eval=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## library('RenvModule')
## module('load','slurm') # Loads slurm among other modules
## library(batchtools)
## myFct <- function(x) {
##     Sys.sleep(10) # to see job in queue, pause for 10 sec
## 	result <- cbind(iris[x, 1:4,],
##                     Node=system("hostname", intern=TRUE),
## 	                Rversion=paste(R.Version()[6:7], collapse="."))
## 	return(result)
##     }


## ----submit_jobs, eval=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## reg <- makeRegistry(file.dir="myregdir", conf.file=".batchtools.conf.R")
## Njobs <- 1:4 # Define number of jobs (here 4)
## ids <- batchMap(fun=myFct, x=Njobs)
## done <- submitJobs(ids, reg=reg, resources=list(partition="short", walltime=120, ntasks=1, ncpus=1, memory=1024))
## waitForJobs() # Wait until jobs are completed


## ----job_status, eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## getStatus() # Summarize job status
## showLog(Njobs[1])
## # killJobs(Njobs) # # Possible from within R or outside with scancel


## ----result_management, eval=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## readRDS("myregdir/results/1.rds") # reads from rds file first result chunk
## loadResult(1)
## lapply(Njobs, loadResult)
## reduceResults(rbind) # Assemble result chunks in single data.frame
## do.call("rbind", lapply(Njobs, loadResult))


## ----remove_registry, eval=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## clearRegistry() # Clear registry in R session
## removeRegistry(wait=0, reg=reg) # Delete registry directory
## # unlink("myregdir", recursive=TRUE) # Same as previous line


## ----load_registry, eval=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## from_file <- loadRegistry("myregdir", conf.file=".batchtools.conf.R")
## reduceResults(rbind)


## ----sessionInfo----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sessionInfo()

