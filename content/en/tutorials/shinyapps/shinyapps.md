---
title: Shiny Web Apps
author: "Thomas Girke"
date: "Last update: 25 May, 2023" 
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
weight: 16
type: docs
editor_options: 
  chunk_output_type: console
---

<!--
- Compile from command-line
Rscript -e "rmarkdown::render('shinyapps.Rmd', c('html_document'), clean=F); knitr::knit('shinyapps.Rmd', tangle=TRUE)"
-->

<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  document.querySelector("h1").className = "title";
});
</script>

<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  var links = document.links;  
  for (var i = 0, linksLength = links.length; i < linksLength; i++)
    if (links[i].hostname != window.location.hostname)
      links[i].target = '_blank';
});
</script>

<div style="text-align: right">

Source code downloads:    
[.Rmd](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/shinyapps/shinyapps.Rmd)    
[.R](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/shinyapps/shinyapps.R)

</div>

## What is Shiny?

[Shiny](https://shiny.rstudio.com/gallery/) is an R-based environment for building interactive web applications for
data analysis and exploration (“Shiny - Tutorial,” n.d.; Chang et al. 2023). Since most JavaScript code is autogenerated by
the environment, basic R knowledge is sufficient for developing Shiny apps.
They can be deployed on local computers or web servers including custom and cloud-based servers (e.g.
AWS, GCP, [shinyapp.io](http://www.shinyapps.io/) service). The basic structure of a Shiny app is an
`app.R` script containing the following components:

1.  User interface
    
    ``` r
    ui <- fluidPage()
    ```

2.  Server function
    
    ``` r
    server <- function(input, output) {}
    ```

3.  Statement to run shiny app
    
    ``` r
    shinyApp(ui = ui, server = server)
    ```

Alternatively, the `ui` and `server` functions can be organized in two script files, a `ui.R` and a `server.R` script, respectively.

## Develop and test Shiny app locally

Open R and set session to parent directory (here `myappdir`) containing shiny script `app.R`, and the
run it with the `runApp()` function. A sample `app.R` script for testing can be downloaded from [here](https://raw.githubusercontent.com/tgirke/GEN242/main/static/custom/scripts/app.R).

``` r
library(shiny)
dir.create("myappdir")
download.file("https://raw.githubusercontent.com/tgirke/GEN242/main/static/custom/scripts/app.R", "./myappdir/app.R")
runApp("myappdir") # To show code in app, add argument: display.mode="showcase" 
```

This will open the app in a web browser.

## Deploy on web server

This can be done on local or cloud systems. An easy solution is to get an account on [shinyapps.io](http://www.shinyapps.io/)
and then deploy Shiny apps there. For details, see [here](https://shiny.rstudio.com/deploy/).

``` r
setwd("myappdir")
library(rsconnect)
deployApp()
```

## Example Shiny app

The following Shiny app is hosted on `shinyapps.io` and embedded into the markdown (or html) source of this page
using the following iframe syntax:

``` r
<iframe src="https://tgirke.shinyapps.io/diamonds/" style="border: none; width: 880px; height: 900px"></iframe>
```

<iframe src="https://tgirke.shinyapps.io/diamonds/" style="border: none; width: 880px; height: 900px">

</iframe>

## Learning Shiny

The Shiny section on the Rstudio site contains excellent [tutorials](https://shiny.rstudio.com/tutorial/).
In addition, users may want to explore the example apps included in the `shiny` package. This can be
done by loading the individual examples (see [here](https://shiny.rstudio.com/tutorial/written-tutorial/lesson1/)) or saving
the code to a user writable directory like this:

``` r
mydir <- system.file("examples", package="shiny")
dir.create('my_shiny_test_dir')
file.copy(mydir, "my_shiny_test_dir", recursive=TRUE)
setwd("my_shiny_test_dir/examples")
runApp("01_hello") # Runs first example app in directory 
dir() # Lists available Shiny examples (directories). 
```

## Resources to learn Shiny

### Tutorial and books

  - Long video [tutorials](https://shiny.rstudio.com/tutorial/).
  - Shiny [official Lessons](https://shiny.rstudio.com/tutorial/written-tutorial/lesson1/).
  - Shiny [official gallery and source code](https://shiny.rstudio.com/gallery/)
  - Advanced Shiny book - [Mastering Shiny](https://mastering-shiny.org/index.html)
  - Advanced web application in R book - [Javascript for R](https://book.javascript-for-r.com/)
  - Shiny Tutorial by Le Zhang (UCR) - [Shiny Tutorial](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/shinyappslezhang/shinyapps)

### Extension packages

  - Catalog of cool extension packages - [Awesome Shiny](https://github.com/nanxstats/awesome-shiny-extensions)
  - [shinyWidgets](https://github.com/dreamRs/shinyWidgets) - UI components
  - [systemPipeShiny](https://systempipe.org/sps/) - A framework for workflow management and data visualization.
  - [spsComps](https://systempipe.org/sps/dev/spscomps/) - UI components, animations, server components
  - [shinyjs](https://deanattali.com/shinyjs/) - server end JavaScript communications

## Session Info

``` r
sessionInfo()
```

    ## R version 4.3.0 (2023-04-21)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Debian GNU/Linux 11 (bullseye)
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
    ##  [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
    ## [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] fgsea_1.26.0     ggplot2_3.4.2    BiocStyle_2.28.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Matrix_1.5-4        gtable_0.3.3        jsonlite_1.8.4      dplyr_1.1.2        
    ##  [5] compiler_4.3.0      BiocManager_1.30.20 tidyselect_1.2.0    Rcpp_1.0.10        
    ##  [9] parallel_4.3.0      jquerylib_0.1.4     scales_1.2.1        BiocParallel_1.34.1
    ## [13] yaml_2.3.7          fastmap_1.1.1       lattice_0.21-8      R6_2.5.1           
    ## [17] generics_0.1.3      knitr_1.42          tibble_3.2.1        munsell_0.5.0      
    ## [21] bslib_0.4.2         pillar_1.9.0        rlang_1.1.1         fastmatch_1.1-3    
    ## [25] utf8_1.2.3          cachem_1.0.8        xfun_0.39           sass_0.4.6         
    ## [29] cli_3.6.1           withr_2.5.0         magrittr_2.0.3      digest_0.6.31      
    ## [33] grid_4.3.0          cowplot_1.1.1       lifecycle_1.0.3     vctrs_0.6.2        
    ## [37] data.table_1.14.8   evaluate_0.21       glue_1.6.2          codetools_0.2-19   
    ## [41] fansi_1.0.4         colorspace_2.1-0    rmarkdown_2.21      tools_4.3.0        
    ## [45] pkgconfig_2.0.3     htmltools_0.5.5

## References

<div id="refs" class="references hanging-indent">

<div id="ref-shiny1">

Chang, Winston, Joe Cheng, JJ Allaire, Carson Sievert, Barret Schloerke, Yihui Xie, Jeff Allen, Jonathan McPherson, Alan Dipert, and Barbara Borges. 2023. *Shiny: Web Application Framework for R*. <https://shiny.rstudio.com/>.

</div>

<div id="ref-noauthor_undated-mi">

“Shiny - Tutorial.” n.d. <https://shiny.rstudio.com/tutorial/>. <https://shiny.rstudio.com/tutorial/>.

</div>

</div>
