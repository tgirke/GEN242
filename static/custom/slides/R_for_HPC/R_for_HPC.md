---
title: "HPC: Module System, Big Data and Parallel Processing"
author: Thomas Girke
date: April 21, 2022
output: 
  ioslides_presentation:
    keep_md: yes
    logo: ./images/ucr_logo.png
    widescreen: yes
    df_print: paged
    smaller: true
subtitle: "Tutorial for R Users" 
bibliography: bibtex.bib
---

<!---
- ioslides manual: 
   https://bookdown.org/yihui/rmarkdown/ioslides-presentation.html

- Compile from command-line
Rscript -e "rmarkdown::render('R_for_HPC.Rmd'); knitr::knit('R_for_HPC.Rmd', tangle=TRUE)"
-->

<!---
  Note: following css chunks are required for scrolling support beyond slide boundaries
-->

<style>
slides > slide {
  overflow-x: auto !important;
  overflow-y: auto !important;
}
</style>

<style type="text/css">
pre {
  max-height: 300px;
  overflow-y: auto;
}

pre[class] {
  max-height: 300px;
}
</style>

<style type="text/css">
.scroll-300 {
  max-height: 300px;
  overflow-y: auto;
  background-color: inherit;
}
</style>

## How to Navigate this Slide Show?

<br/>

- This __ioslides__ presentation contains scrollable slides. 
- Which slides are scrollable, is indicated by a tag at the bottom of the corresponding slides stating: 

<p style='text-align: center;'> __[ Scroll down to continue ]__ </p>

- The following single character keyboard shortcuts enable alternate display modes of __ioslides__:
    - `f`: enable fullscreen mode
    - `w`: toggle widescreen mode
    - `o`: enable overview mode
    - `h`: enable code highlight mode
- Pressing Esc exits all of these modes. Additional details can be found [here](https://bookdown.org/yihui/rmarkdown/ioslides-presentation.html).

# Outline

- <div class="white">__Another Nvim tip: mouse support__</div>
- Tmux review 
- Module system
- Big data storage
- Parallel processing and queuing system
- Parallel R with _batchtools_
- References

## Another Nvim tip: mouse support

The following shows how to enable in Nvim mouse support. When enabled one can position 
the cursor anywhere with the mouse as well as resize split windows, and switch 
the scope from one window split to another. 

- To enable mouse support, type in Nvim's command mode:
    - `:set mouse=a`

- To toggle back to no mouse support, type in command mode:
    - `:set mouse-=a`

- To enable mouse support by default, add `set mouse=a` to Nvim's config file located in a user's home 
under: `~/.config/nvim/init.vim`

- To find more help on this topic, type in nvim’s command mode:
    - `:help mouse`

# Outline

- Another Nvim tip: mouse support
- <div class="white">__Tmux review__</div>
- Module system
- Big data storage
- Parallel processing and queuing system
- Parallel R with _batchtools_
- References

## Tmux for managing terminal sessions

### What is Tmux?

- Tmux is a virtual terminal multiplexer providing re-attachable terminal sessions
- Advantage: work in a terminal session cannot get lost due to internet disruptions or even when switching computers
- Combined with the `Nvim-r` plugin it provides a flexible working environment for R 
- Users can send code from a script to the R console or command-line.
- On HPCC both Nvim-R and Tmux are pre-configured and easy to install 

<center><img title="Nvim-R" src="https://raw.githubusercontent.com/jalvesaq/Nvim-R/master/Nvim-R.gif" style="width:400px;"></center>


## Typical Usage Workflow for Nvim-R-Tmux


__1. Start tmux session from login node (not compute node!)__

Running Nvim from tmux provides reattachment functionality. Skip this step if this is not required.


```bash
tmux # starts a new tmux session 
tmux a # attaches to an existing or preconfigured session 
```

__2. Open nvim-connected R session__ 

Open a `*.R` or `*.Rmd` file with `nvim` and initialize a connected R session
with `\rf`. Note, the resulting split window among Nvim and R behaves like a split
viewport in `nvim` or `vim` meaning the usage of `Ctrl-w w` followed by `i` and
`Esc` is important for session navigation.


```bash
nvim myscript.R # or *.Rmd file
```

__3. Send R code from nvim to the R pane__

Single lines of code can be sent from nvim to the R console by pressing the space bar. To send 
several lines at once, one can select them in nvim's visual mode and then hit the space bar. 

## Keybindings to Control Environment

### Important keybindings for nvim

* `\rf`: opens vim-connected R session. If you do this the first time in your user account, you might be asked to create an `R` directory under `~/`. If so approve this action by pressing `y`. 
* `spacebar`: sends code from vim to R; here remapped in `init.vim` from default `\l`
* `:split` or `:vsplit`: splits viewport (similar to pane split in tmux)
* `gz`: maximizes size of viewport in normal mode (similar to Tmux's `Ctrl-a z` zoom utility) 
* `Ctrl-w w`: jumps cursor to R viewport and back; toggle between insert (`i`) and command (`Esc`) mode is required for navigation and controlling the environment.
* `Ctrl-w r`: swaps viewports
* `Ctrl-w =`: resizes splits to equal size
* `:resize <+5 or -5>`: resizes height by specified value

<br/><br/>
<p style='text-align: right;'> __[ Scroll down to continue ]__ </p>
<br/><br/><br/><br/>

* `:vertical resize <+5 or -5>`: resizes width by specified value
* `Ctrl-w H` or `Ctrl-w K`: toggles between horizontal/vertical splits
* `Ctrl-spacebar`: omni completion for R objects/functions when nvim is in insert mode. Note, this has been remapped in `init.vim` from difficult to type default `Ctrl-x Ctrl-o`. 
* `:h nvim-R`: opens nvim-R's user manual; navigation works the same as for any Vim/Nvim help document
* `:Rhelp fct_name`: opens help for a function from nvim's command mode with text completion support
* `Ctrl-s and Ctrl-x`: freezes/unfreezes vim (some systems)


### Important keybindings for tmux

__Pane-level commands__

* `Ctrl-a %`: splits pane vertically
* `Ctrl-a "`: splits pane horizontally
* `Ctrl-a o`: jumps cursor to next pane
* `Ctrl-a Ctrl-o`: swaps panes
* `Ctrl-a <space bar>`: rotates pane arrangement
* `Ctrl-a Alt <left or right>`: resizes to left or right
* `Ctrl-a Esc <up or down>`: resizes to left or right

__Window-level comands__

* `Ctrl-a n`: switches to next tmux window 
* `Ctrl-a Ctrl-a`: switches to previous tmux window
* `Ctrl-a c`: creates a new tmux window 
* `Ctrl-a 1`: switches to specific tmux window selected by number

__Session-level comands__

* `Ctrl-a d`: detaches from current session
* `Ctrl-a s`: switch between available tmux sessions
* `$ tmux new -s <name>`: starts new session with a specific name
* `$ tmux ls`: lists available tmux session(s)
* `$ tmux attach -t <id>`: attaches to specific tmux session  
* `$ tmux attach`: reattaches to session 
* `$ tmux kill-session -t <id>`: kills a specific tmux session
* `Ctrl-a : kill-session`: kills a session from tmux command mode 

# Outline

- Another Nvim tip: mouse support
- Tmux review
- <div class="white">__Module system__</div>
- Big data storage
- Parallel processing and queuing system
- Parallel R with _batchtools_
- References

## Software and module system on HPCC

* Over 2,000 software tools are currently installed on HPCC Cluster
* Custom installs in user accounts via various mechanisms, e.g. environment management systems such as [conda](https://conda.io/projects/conda/en/latest/index.html)
* Most common research databases used in bioinformatics are available
* Support of most common programming languages used in research computing
* A module system is used to facilitate the management of software tools. This includes any number of versions of each software.
* New software install requests can be sent to support@hpcc.ucr.edu.
* To use software manged under the module system, users need to learn using some basic commands. The most common commands are listed below.

Print available modules
```sh
module avail
```

<br/><br/>
<p style='text-align: right;'> __[ Scroll down to continue ]__ </p>
<br/><br/><br/><br/>

Print available modules starting with R
```sh
module avail R
```

Load default module R
```sh
module load R
```

Load specific module R version
```sh
module load R/4.1.2
```

List loaded modules
```sh
module list
```

Unload module R
```sh
module unload R
```

Unload specific module R
```sh
module unload R/4.1.3
```

# Outline

- Another Nvim tip: mouse support
- Tmux review
- Module system
- <div class="white">__Big data storage__</div>
- Parallel processing and queuing system
- Parallel R with _batchtools_
- References

## Big data storage

Each user account on HPCC Cluster comes only with 20GB of disk space. Much more disk space is 
available in a dedicated `bigdata` directory. How much space depends on the subscription 
of each user group. The path of `bigdata` and `bigdata-shared` is as follows:

* `/bigdata/labname/username`
* `/bigdata/labname/shared`

All lab members share the same bigdata pool. The course number `gen242` is used as `labname`
for user accounts adminstered under GEN242 (here /bigdata/gen242/shared).

The disk usage of `home` and `bigdata` can be monitored on the [HPCC Cluster Dashboard](https://dashboard.hpcc.ucr.edu/).

Additional details can be found on the Project Data page of GEN242 [here](https://girke.bioinformatics.ucr.edu/GEN242/assignments/projects/project_data/).

# Outline

- Another Nvim tip: mouse support
- Tmux review
- Module system
- Big data storage
- <div class="white">__Parallel processing and queuing system__</div>
- Parallel R with _batchtools_
- References

## Queuing system: `Slurm` 

HPCC Cluster uses `Slurm` as queuing and load balancing system. To control user traffic, any 
type of compute intensive jobs need to be submitted via `sbatch` or `srun` (see below) to the computer
nodes. Much more detailed information on this topic can be found on these sites: 

+ [UCR HPCC Manual](https://hpcc.ucr.edu/manuals_linux-cluster_jobs.html)
+ [Slurm Documentation](https://slurm.schedmd.com/documentation.html)
+ [Torque/Slurm Comparison](http://www.nersc.gov/users/computational-systems/cori/running-jobs/for-edison-users/torque-moab-vs-slurm-comparisons/)
+ [Switching from Torque to Slurm](https://sites.google.com/a/case.edu/hpc-upgraded-cluster/slurm-cluster-commands)
+ [Slurm Quick Start Tutorial](http://www.ceci-hpc.be/slurm_tutorial.html)

### Job submission with `sbatch`

Print information about queues/partitions available on a cluster.
```sh
sinfo
```

<br/><br/>
<p style='text-align: right;'> __[ Scroll down to continue ]__ </p>
<br/><br/><br/><br/>

Compute jobs are submitted with `sbatch` via a submission script (here `script_name.sh`).

```sh
sbatch script_name.sh
```

The following sample submission script (`script_name.sh`) executes an R script named `my_script.R`.

```sh
#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:15:00 # 1 day and 15 minutes
#SBATCH --mail-user=useremail@address.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="some_test"
#SBATCH -p batch # Choose queue/parition from: intel, batch, highmem, gpu, short

Rscript my_script.R
```

`STDOUT` and `STDERROR` of jobs will be written to files named
`slurm-<jobid>.out` or to a custom file specified under `#SBATCH --output` in
the submission script. 

### Interactive sessions with `srun`

This option logs a user in to a computer node of a specified partition (queue), while Slurm monitors and controls the resource request.

```sh
srun --pty bash -l
```

Interactive session with specific resource requests
```sh
srun --x11 --partition=short --mem=2gb --cpus-per-task 4 --ntasks 1 --time 1:00:00 --pty bash -l
```

The argument `--mem` limits the amount of RAM, `--cpus` the number of CPU
cores, `--time` the time how long a session will be active. Under
`--parition` one can choose among different queues and node architectures.
Current options under `--partition` for most users of the HPCC cluster are: `intel`, `batch`, `highmem`, `gpu`,
and `short`. The latter has a time limit of 2 hours. 


### Monitoring jobs with `squeue`

List all jobs in queue
```sh
squeue
```

List jobs of a specific user
```sh
squeue -u <user>
```

Print more detailed information about a job
```sh
scontrol show job <JOBID>
```

Custom command to summarize and visualize cluster activity
```sh
jobMonitor
```

### Deleting and altering jobs 

Delete a single job
```sh
scancel -i <JOBID>
```

Delete all jobs of a user
```sh
scancel -u <username> 
```

Delete all jobs of a certain name
```sh
scancel --name <myJobName>
```

Altering jobs with `scontrol update`. The below example changes the walltime (`<NEW_TIME>`) of a specific job (`<JOBID>`). 
```sh
scontrol update jobid=<JOBID> TimeLimit=<NEW_TIME>
```

### Resource limits

Resourse limits for users can be viewed as follows. 
```sh
sacctmgr show account $GROUP format=Account,User,Partition,GrpCPUs,GrpMem,GrpNodes --ass | grep $USER
```

Similarly, one can view the limits of the group a user belongs to. 
```sh
sacctmgr show account $GROUP format=Account,User,Partition,GrpCPUs,GrpMem,GrpNodes,GrpTRES%30 --ass | head -3
```

# Outline

- Another Nvim tip: mouse support
- Tmux review
- Module system
- Big data storage
- Parallel processing and queuing system
- <div class="white">__Parallel R with _batchtools_ __</div>
- References

## Parallel Evaluations in R 

- <b><font color="red">Note:</font></b> the content on the following slides is also available in this tutorial section [here](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/rparallel/rparallel/).
- R provides a large number of packages for parallel evaluations on multi-core, multi-socket and multi-node systems. The latter are usually referred to as computer clusters.
- MPI is also supported
- For an overview of parallelization packages available for R see [here](https://cran.r-project.org/web/views/HighPerformanceComputing.html)
- One of the most comprehensive parallel computing environments for R is
  [`batchtools`](https://mllg.github.io/batchtools/articles/batchtools.html#migration). Older versions of this package were released under the name `BatchJobs` [@Bischl2015-rf]. 
- `batchtools` supports both multi-core and multi-node computations with and without schedulers. By making use of
  cluster template files, most schedulers and queueing systems are supported (_e.g._ Torque, Sun Grid Engine, Slurm). 
- The `BiocParallel` package (see [here](https://bioconductor.org/packages/release/bioc/html/BiocParallel.html)) 
  provides similar functionalities as `batchtools`, but is tailored to use Bioconductor objects.

## Reminder: Traditional Job Submission for R 

This topic is covered in more detail in the basic Linux/HPC tutorial [here](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/linux/linux/#queuing-system-slurm). 
Briefly, the following shows how to submit a script for precessing to the computing nodes. 

__1.__ Create Slurm submission script, here called [script_name.sh](https://raw.githubusercontent.com/tgirke/GEN242/main/static/custom/slides/R_for_HPC/demo_files/script_name.sh) with:


```bash
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:15:00 # 1 day and 15 minutes
#SBATCH --mail-user=useremail@address.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="some_test"
#SBATCH -p short # Choose queue/partition from: intel, batch, highmem, gpu, short

Rscript my_script.R
```

__2.__ Submit R script called [my_script.R](https://raw.githubusercontent.com/tgirke/GEN242/main/static/custom/slides/R_for_HPC/demo_files/my_script.R) by above Slurm script with:


```bash
sbatch script_name.sh
```

## Parallel Evaluations on Clusters with `batchtools` 

- The following introduces the usage of `batchtools` for a computer cluster
  using SLURM as scheduler (workload manager). SLURM is the scheduler used by
  the HPCC at UCR.
- Similar instructions are provided this tutorial section covering
  `batchtools`
  [here](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/rparallel/rparallel/)
- To simplify the evaluation of the R code on the following slides, the
  corresponding text version is available for download from
  [here](https://raw.githubusercontent.com/tgirke/GEN242/main/static/custom/slides/R_for_HPC/demo_files/R_for_HPC_demo.R).

## Hands-on Demo of `batchtools`

### Set up working directory for SLURM
First login to your cluster account, open R and execute the following lines. This will
create a test directory (here `mytestdir`), redirect R into this directory and then download
the required files: 

+ [`slurm.tmpl`](https://github.com/ucr-hpcc/ucr-hpcc.github.io/blob/master/presentations/2020-12-18_Workshop/R_for_HPC/demo_files/slurm.tmpl)
+ [`.batchtools.conf.R`](https://github.com/ucr-hpcc/ucr-hpcc.github.io/blob/master/presentations/2020-12-18_Workshop/R_for_HPC/demo_files/.batchtools.conf.R)


```r
dir.create("mytestdir")
setwd("mytestdir")
download.file("https://bit.ly/3gZJBsy", "slurm.tmpl")
download.file("https://bit.ly/3nvSNHA", ".batchtools.conf.R")
```

### Load package and define some custom function

The following code defines a test function (here `myFct`) that will be run on the cluster for demonstration
purposes. 

<p style='text-align: right;'> __[ Scroll down to continue ]__ </p>
<br/><br/>

The test function (`myFct`) subsets the `iris` data frame by rows, and appends the host name and R version of each
node where the function was executed. The R version to be used on each node can be
specified in the `slurm.tmpl` file (under `module load`).


```r
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


```r
reg <- makeRegistry(file.dir="myregdir", conf.file=".batchtools.conf.R")
Njobs <- 1:4 # Define number of jobs (here 4)
ids <- batchMap(fun=myFct, x=Njobs) 
done <- submitJobs(ids, reg=reg, resources=list(partition="short", walltime=120, ntasks=1, ncpus=1, memory=1024))
waitForJobs() # Wait until jobs are completed
```

### Summarize job status 
After the jobs are completed one can inspect their status as follows.


```r
getStatus() # Summarize job status
showLog(Njobs[1])
# killJobs(Njobs) # # Possible from within R or outside with scancel
```

### Access/assemble results

The results are stored as `.rds` files in the registry directory (here `myregdir`). One
can access them manually via `readRDS` or use various convenience utilities provided
by the `batchtools` package.


```r
readRDS("myregdir/results/1.rds") # reads from rds file first result chunk
loadResult(1) 
lapply(Njobs, loadResult)
reduceResults(rbind) # Assemble result chunks in single data.frame
do.call("rbind", lapply(Njobs, loadResult))
```

### Remove registry directory from file system

By default existing registries will not be overwritten. If required one can explicitly
clean and delete them with the following functions. 


```r
clearRegistry() # Clear registry in R session
removeRegistry(wait=0, reg=reg) # Delete registry directory
# unlink("myregdir", recursive=TRUE) # Same as previous line
```

### Load registry into R 

Loading a registry can be useful when accessing the results at a later state or 
after moving them to a local system. 


```r
from_file <- loadRegistry("myregdir", conf.file=".batchtools.conf.R")
reduceResults(rbind)
```

## Conclusions

### Advantages of `batchtools`

- many parallelization methods: multiple cores, and across both multiple CPU sockets and nodes
- most schedulers supported
- takes full advantage of a cluster
- robust job management by organizing results in registry file-based database
- simplifies submission, monitoring and restart of jobs 
- well supported and maintained package


# Outline

- Another Nvim tip: mouse support
- Tmux review
- Module system
- Big data storage
- Parallel processing and queuing system
- Parallel R with _batchtools_ 
- <div class="white">__References__</div>


## References 

