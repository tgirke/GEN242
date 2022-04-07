## pre {

##   max-height: 300px;

##   overflow-y: auto;

## }

## 
## pre[class] {

##   max-height: 300px;

## }


## .scroll-300 {

##   max-height: 300px;

##   overflow-y: auto;

##   background-color: inherit;

## }


## nvim myfile.txt # for neovim (or 'vim myfile.txt' for vim)


## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Install_Nvim-R_Tmux


## tmux # starts a new tmux session

## tmux a # attaches to an existing or preconfigured session


## nvim myscript.R # or *.Rmd file


## srun --x11 --partition=short --mem=2gb --cpus-per-task 4 --ntasks 1 --time 1:00:00 --pty bash -l


## wget https://raw.githubusercontent.com/tgirke/GEN242/main/static/custom/slides/R_for_HPC/demo_files/R_for_HPC_demo.R


## ----nvim-r-tmux-demo_show, eval=FALSE, message=FALSE, warning=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------
## library(tidyverse)
## write_tsv(iris, "iris.txt") # Creates sample file
## read_tsv("iris.txt") %>% # Import with read_tbv from readr package
##     as_tibble() %>%
##     group_by(Species) %>%
##     summarize_all(mean) %>%
##     reshape2::melt(id.vars=c("Species"), variable.name = "Samples", value.name="Values") %>%
##     ggplot(aes(Samples, Values, fill = Species)) +
##     geom_bar(position="dodge", stat="identity")


## ----nvim-r-tmux-demo_run, echo=FALSE, eval=TRUE, message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)                                                                                                                                                            
write_tsv(iris, "iris.txt") # Creates sample file                                                                                                                             
read_tsv("iris.txt") %>% # Import with read_tbv from readr package                                                                                                            
    as_tibble() %>%                                                                                                                                                           
    group_by(Species) %>%                                                                                                                                                     
    summarize_all(mean) %>%                                                                                                                                                   
    reshape2::melt(id.vars=c("Species"), variable.name = "Samples", value.name="Values") %>%                                                                                  
    ggplot(aes(Samples, Values, fill = Species)) +                                                                                                                            
    geom_bar(position="dodge", stat="identity")


## module avail R


## module load R/4.1.2


## module list

