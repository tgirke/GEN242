## ----style, echo = FALSE, results = 'asis'--------------------------------------------------------
BiocStyle::markdown()
options(width=100, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")))


## ----setup, echo=FALSE, messages=FALSE, warnings=FALSE--------------------------------------------
suppressPackageStartupMessages({
    library(ggplot2)
    library(fgsea)
})


## ----fluidpage, eval=FALSE------------------------------------------------------------------------
## ui <- fluidPage()


## ----server, eval=FALSE---------------------------------------------------------------------------
## server <- function(input, output) {}


## ----shinyapp, eval=FALSE-------------------------------------------------------------------------
## shinyApp(ui = ui, server = server)


## ----runshinyapp1, eval=FALSE---------------------------------------------------------------------
## library(shiny)
## runApp("myappdir") # To show code in app, add argument: display.mode="showcase"


## ----deployshinyapp1, eval=FALSE------------------------------------------------------------------
## setwd("myappdir")
## library(rsconnect)
## deployApp()


## ----embedshiny, eval=FALSE-----------------------------------------------------------------------
## <iframe src="https://tgirke.shinyapps.io/diamonds/" style="border: none; width: 880px; height: 900px"></iframe>


## ----learnshiny, eval=FALSE-----------------------------------------------------------------------
## mydir <- system.file("examples", package="shiny")
## dir.create('my_shiny_test_dir')
## file.copy(mydir, "my_shiny_test_dir", recursive=TRUE)
## setwd("my_shiny_test_dir/examples")
## runApp("01_hello") # Runs first example app in directory
## dir() # Lists available Shiny examples (directories).


## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()

