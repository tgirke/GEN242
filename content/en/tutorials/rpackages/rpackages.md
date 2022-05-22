---
title: "Building R Packages" 
author: "Author: Thomas Girke"
date: "Last update: 22 May, 2022" 
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
weight: 17
type: docs
---

<!--
- Compile from command-line
Rscript -e "rmarkdown::render('rpackages.Rmd', c('html_document'), clean=F); knitr::knit('rpackages.Rmd', tangle=TRUE)"
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
\[ [Slides](https://girke.bioinformatics.ucr.edu/GEN242/slides/slides_19/) \]    
\[ [.Rmd](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/rpackages/rpackages.Rmd) \]    
\[ [.R](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/rpackages/rpackages.R) \]

</div>

## Overview

### Motivation for building R packages

1.  Organization
      - Consolidate functions with related utilties in single place  
      - Interdepencies among less complex functions make coding more efficient
      - Minimizes duplications
2.  Documentation
      - Help page infrastructure improves documentation of functions
      - Big picture of utilties provided by package vignettes (manuals)
3.  Sharability
      - Package can be easier shared with colleagues and public
      - Increases code accessibilty for other users
4.  Extendibility
      - Makes software more extentible and maintainable

### Package development environments

This page introduces two approaches for building R packages:

1.  `R Base` and related functionalities
2.  `devtools` and related packages (*e.g.* `usethis`, `roxygen2` and `sinew`)

The sample code provided below creates for each of the two methods a simple test package
that can be installed and loaded on a user’s system. The instructions for the
second appoach are more detailed since it is likely to provide the most
practical solution for newer users of R.

## 1\. R Base Approach

R packages can be built with the `package.skeleton` function. The most
comprehensive documentation on package development is provided by the [Writing
R
Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Package-Dependencies)
page on CRAN. The basic workflow example below will create a directory named
`mypackage` containing the skeleton of the package for all functions, methods
and classes defined in the R script(s) passed on to the `code_files` argument.
The basic structure of the package directory is described
[here](http://cran.fhcrc.org/doc/manuals/R-exts.html#Package-structure).
The package directory will also contain a file named `Read-and-delete-me` with
instructions for completing the package:

### 1.1 Create package skeleton

``` r
## Download R script (here pkg_build_fct.R) containing two sample functions
download.file("https://raw.githubusercontent.com/tgirke/GEN242/main/content/en/tutorials/rpackages/helper_functions/pkg_build_fct.R", "pkg_build_fct.R")
## Build package skeleton based on functions in pkg_build_fct.R
package.skeleton(name="mypackage", code_files=c("pkg_build_fct.R"))
```

The given example will create a directory named `mypackage` containing the
skeleton of the package for all functions, methods and classes defined in the R
script(s) passed on to the `code_files` argument. The basic structure of the
package directory is described [here](http://cran.fhcrc.org/doc/manuals/R-exts.html#Package-structure).
The package directory will also contain a file named ‘Read-and-delete-me’ with the following
instructions for completing the package:

  - Edit the help file skeletons in man, possibly combining help files for multiple functions.
  - Edit the exports in NAMESPACE, and add necessary imports.
  - Put any C/C++/Fortran code in src.
  - If you have compiled code, add a `useDynLib()` directive to `NAMESPACE`.
  - Run R CMD build to build the package tarball.
  - Run R CMD check to check the package tarball.
  - Read Writing R Extensions for more information.

### 1.2 Build package

Once a package skeleton is available one can build the package from the
command-line (Linux/OS X), or from within R by executing the command-line calls
with R’s `system("...")` command. For instance, the command-line call `R CMD build ...`
can be run from within R with `system("R CMD build ...")`. The following examples
refer to the R console.

``` r
system("R CMD build mypackage")
```

### 1.3 Check package

This will create a tarball of the package with its version number encoded in
the file name (here `mypackage_1.0.tar.gz`). Subsequently, the package tarball
needs to be checked for errors with `R CMD check`.

``` r
system("R CMD check mypackage_1.0.tar.gz")
```

All issues in a package’s source code and documentation should be addressed
until `R CMD check` returns no error or warning messages anymore.

### 1.4 Install package

Install package from source on Linux or OS X systems.

``` r
install.packages("mypackage_1.0.tar.gz", repos=NULL) # On OS X include the argument type="source"
```

Windows requires a zip archive for installing R packages, which can be most
conveniently created from the command-line (Linux/OS X) by installing the
package in a local directory (here `tempdir`) and then creating a zip archive
from the installed package directory:

``` r
$ mkdir tempdir
$ R CMD INSTALL -l tempdir mypackage_1.0.tar.gz
$ cd tempdir
$ zip -r mypackage mypackage
## The resulting mypackage.zip archive can be installed from R under Windows like this:
install.packages("mypackage.zip", repos=NULL)
```

This procedure only works for packages which do not rely on compiled code (C/C++).
Instructions to fully build an R package under Windows can be found [here](http://www.murdoch-sutherland.com/Rtools/) and
[here](https://www.r-bloggers.com/2010/06/how-to-create-an-r-package-in-windows/).

### 1.5 Maintain package

Several useful helper utilities exist for maintaing and extending packages. Typical package development routines
include:

  - Adding new functions, methods and classes to the script files in the ./R directory in your package
  - Adding their names to the NAMESPACE file of the package
  - Additional `.Rd` help templates can be generated with the `prompt()` function family like this:

<!-- end list -->

``` r
source("myscript.R") # imports functions, methods and classes from myscript.R
prompt(myfct) # writes help file myfct.Rd
promptClass("myclass") # writes file myclass-class.Rd
promptMethods("mymeth") # writes help file mymeth.Rd
```

The resulting `.Rd` help files can be edited in a text editor, and properly rendered and viewed from within R with help of
the following functions.

``` r
library(tools)
Rd2txt("./mypackage/man/myfct.Rd") # renders *.Rd files as they look in final help pages
checkRd("./mypackage/man/myfct.Rd") # checks *.Rd help file for problems
```

### 1.6 Submit package to public repository

The best way of sharing an R package with the community is to submit it to one
of the main R package repositories, such as CRAN or Bioconductor. The details
about the submission process are given on the corresponding repository
submission pages:

  - [Submitting to Bioconductor](http://www.bioconductor.org/developers/package-submission/)
  - [Submitting to CRAN](https://cran.r-project.org/)

## 2\. R `devtools` Approach

Several package develpment routines of the traditional method outlined above
are manual, such as updating the NAMESPACE file and documenting functions in
separate help (`*.Rd`) files. This process can be simplified and partially
automated by taking advantage of a more recent R package development
environment composed of several helper packages including `devtools`,
`usethis`, `roxygen2` and `sinew` (Wickham and Bryan, n.d.). Many books and web sites document this process
in more detail. Here is a small selection of useful online documentation about
R package development:

  - Book: [R Packages](https://r-pkgs.org/index.html) by Hadley Wickham and Jenny Bryan  
  - [My First R Package](https://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html) by Fong Chun Chan
  - [How to Creat an R Package, Easy Mode](https://www.amitkohli.com/2020/01/07/2020-01-07-how-to-create-an-r-package-my-way/) by Amit Kohli
  - [Package Development Cheat Sheet](https://rawgit.com/rstudio/cheatsheets/master/package-development.pdf)
  - Automating `roxygen2` documentation with `sinew` by Jonathan Sidi: [Blog](https://yonicd.github.io/2017-09-18-sinew/) and [CRAN](https://cran.r-project.org/web/packages/sinew/index.html)

## Workflow for building R packages

The following outlines the basic workflow for building, testing and extending R packages with
the package development environment functionalities outlined above.

### 2.1 Create package skeleton

``` r
library("devtools"); library("roxygen2"); library("usethis"); library(sinew) # If not availble install these packages with 'install.packages(...)'
create("myfirstpkg") # Creates package skeleton. The chosen name (here myfirstpkg) will be the name of the package.
setwd("myfirstpkg") # Set working directory of R session to package directory 'myfirstpkg'
use_mit_license() # Add license information to description file (here MIT). To look up alternatives, do ?use_mit_license
```

### 2.2 Add R functions

Next, R functions can be added to `*.R` file(s) under the R directory of the new
package. Several functions can be organized in one `*.R` file, each in its own
file or any combination. For demonstration purposes, the following will
download an R file (`pkg_build_fct.R` from [here](https://raw.githubusercontent.com/tgirke/GEN242/main/content/en/tutorials/rpackages/helper_functions/pkg_build_fct.R)) defining two functions (named:`myMAcomp`
and `talkToMe`) and save it to the R directory of the package.

``` r
download.file("https://raw.githubusercontent.com/tgirke/GEN242/main/content/en/tutorials/rpackages/helper_functions/pkg_build_fct.R", "R/pkg_build_fct.R")
```

### 2.3 Auto-generate roxygen comment lines

The `makeOxygen` function from the `sinew` package creates `roxygen2` comment
skeletons based on the information from each function (below for `myMAcomp` example).
The roxygen comment lines need to be added above the code of each function.
This can be done by copy and paste from the R console or by writing the output
to a temporary file (below via `writeLines`). Alternatively, the `makeOxyFile`
function can be used to create a roxygenized copy of an R source file, where the
roxygen comment lines have been added above all functions automatically. Next,
the default text in the comment lines needs to be replaced by meaningful text
describing the utility and usage of each function. This editing process of documentation
can be completed and/or revised any time.

``` r
load_all() # Loads package in a simulated way without installing it. 
writeLines(makeOxygen(myMAcomp), "myroxylines") # This creates a 'myroxylines' file in current directory. Delete this file after adding its content to the corresponding functions.
```

### 2.4 Autogenerate help files

The `document` function autogenerates for each function one `*.Rd` file in the
`man` directory of the package. The content in the `*.Rd` help files is based
on the information in the roxygen comment lines generated in the
previous step. In addition, all relevant export/import instructions are added
to the `NAMESPACE` file. Importantly, when using roxygen-based documentation in a package
then the `NAMESPACE` and `*.Rd` files should not be manually edited since this information
will be lost during the automation routines provided by `roxygen2`.

``` r
document() # Auto-generates/updates *.Rd files under man directory (here: myMAcomp.Rd and talkToMe.Rd)
tools::Rd2txt("man/myMAcomp.Rd") # Renders Rd file from source
tools::checkRd("man/myMAcomp.Rd") # Checks Rd file for problems
```

### 2.5 Add a vignette

A vignette template can be auto-generated with the `use_vignette` function from
the `usethis` package. The `*.Rmd` source file of the vignette will be located
under a new `vignette` directory. Additional vignettes can be manually added to
this directory as needed.

``` r
use_vignette("introduction", title="Introduction to this package")
```

### 2.6 Check, install and build package

Now the package can be checked for problems. All warnings and errors should be
addressed prior to submission to a public repository. After this it can be
installed on a user’s system with the `install` command. In addition, the
`build` function allows to assemble the package in a `*.tar.gz` file. The
latter is often important for sharing packages and/or submitting them to public
repositories.

``` r
setwd("..") # Redirect R session to parent directory
check("myfirstpkg") # Check package for problems, when in pkg dir one can just use check()
# remove.packages("myfirstpkg") # Optional. Removes test package if already installed
install("myfirstpkg", build_vignettes=TRUE) # Installs package  
build("myfirstpkg") # Creates *.tar.gz file for package required to for submission to CRAN/Bioc
```

### 2.7 Using the new package

After installing and loading the package its functions, help files and
vignettes can be accessed as follows.

``` r
library("myfirstpkg")
library(help="myfirstpkg")
?myMAcomp
vignette("introduction", "myfirstpkg")
```

Another very useful development function is `test` for evaluating the test code of a package.

### 2.8 Share package on GitHub

To host and share the new package `myfirstpkg` on GitHub, one can use the following steps:

1.  Create an empty target GitHub repos online (*e.g.* named `mypkg_repos`) as outlined [here](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/github/github/#exercise).
2.  Clone the new GitHub repos to local system with `git clone https://github.com/<github_username>/<repo name>` (here from command-line)
3.  Copy the root directory of the package into the downloaded repos with `cp -r myfirstpkg mypkg_repos`  
4.  Next `cd` into `mypkg_repos`, and then add all files and directories of the package to the staging area with `git add -A :/`.
5.  Commit and push the changes to GitHub with: `git commit -am "first commit"; git push`.
6.  After this the package should be life on the corresponding GitHub web page.
7.  Assuming the package is public, it can be installed directly from GitHub by anyone as shown below (from within R). Installs of
    private packages require a personal access token (PAT) that needs to be assigned to the `auth_token` argument. PATs can be created
    [here](https://github.com/settings/tokens).

<!-- end list -->

``` r
devtools::install_github("<github_user_name>/<mypkg_repos>", subdir="myfirstpkg") # If the package is in the root directory of the repos, then the 'subdir' argument can be dropped.
```

## Session Info

``` r
sessionInfo()
```

    ## R version 4.2.0 (2022-04-22)
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
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] ggplot2_3.3.6    limma_3.52.0     BiocStyle_2.24.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bslib_0.3.1         compiler_4.2.0      pillar_1.7.0        BiocManager_1.30.17
    ##  [5] jquerylib_0.1.4     tools_4.2.0         digest_0.6.29       jsonlite_1.8.0     
    ##  [9] evaluate_0.15       lifecycle_1.0.1     tibble_3.1.7        gtable_0.3.0       
    ## [13] pkgconfig_2.0.3     rlang_1.0.2         DBI_1.1.2           cli_3.3.0          
    ## [17] yaml_2.3.5          blogdown_1.9        xfun_0.30           fastmap_1.1.0      
    ## [21] withr_2.5.0         dplyr_1.0.9         stringr_1.4.0       knitr_1.39         
    ## [25] generics_0.1.2      sass_0.4.1          vctrs_0.4.1         tidyselect_1.1.2   
    ## [29] grid_4.2.0          glue_1.6.2          R6_2.5.1            fansi_1.0.3        
    ## [33] rmarkdown_2.14      bookdown_0.26       purrr_0.3.4         magrittr_2.0.3     
    ## [37] codetools_0.2-18    scales_1.2.0        htmltools_0.5.2     ellipsis_0.3.2     
    ## [41] assertthat_0.2.1    colorspace_2.0-3    utf8_1.2.2          stringi_1.7.6      
    ## [45] munsell_0.5.0       crayon_1.5.1

## References

<div id="refs" class="references hanging-indent">

<div id="ref-Wickham_undated-ei">

Wickham, Hadley, and Jennifer Bryan. n.d. “R Packages.” <https://r-pkgs.org/index.html>. <https://r-pkgs.org/index.html>.

</div>

</div>
