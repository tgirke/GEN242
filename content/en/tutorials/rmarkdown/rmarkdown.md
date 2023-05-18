---
title: R Markdown Tutorial
author: "Author: Your Name"
date: "Last update: 17 May, 2023" 
output:
  BiocStyle::html_document:
    toc: true
    toc_float:
        collapsed: true
        smooth_scroll: true
    toc_depth: 3
    fig_caption: yes
    code_folding: show
    number_sections: false
fontsize: 14pt
bibliography: bibtex.bib
type: docs
weight: 11
---

<script src="/rmarkdown-libs/htmlwidgets/htmlwidgets.js"></script>

<link href="/rmarkdown-libs/datatables-css/datatables-crosstalk.css" rel="stylesheet" />

<script src="/rmarkdown-libs/datatables-binding/datatables.js"></script>

<script src="/rmarkdown-libs/jquery/jquery-3.6.0.min.js"></script>

<link href="/rmarkdown-libs/dt-core/css/jquery.dataTables.min.css" rel="stylesheet" />
<link href="/rmarkdown-libs/dt-core/css/jquery.dataTables.extra.css" rel="stylesheet" />

<script src="/rmarkdown-libs/dt-core/js/jquery.dataTables.min.js"></script>

<link href="/rmarkdown-libs/crosstalk/css/crosstalk.min.css" rel="stylesheet" />

<script src="/rmarkdown-libs/crosstalk/js/crosstalk.min.js"></script>

<!---
- Compile from command-line
Rscript -e "rmarkdown::render('rmarkdown.Rmd', c('BiocStyle::html_document'), clean=FALSE); knitr::knit('rmarkdown.Rmd', tangle=TRUE)"
-->

<div style="text-align: right">

Source code downloads:    
\[ [.Rmd](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/rmarkdown/rmarkdown.Rmd) \]    
\[ [.html](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/rmarkdown/rmarkdown.html) \]    
\[ [.R](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/rmarkdown/rmarkdown.R) \]

</div>

## R Markdown Overview

R Markdown combines markdown (an easy to write plain text format) with embedded
R code chunks. When compiling R Markdown documents, the code components can be
evaluated so that both the code and its output can be included in the final
document. This makes analysis reports highly reproducible by allowing to automatically
regenerate them when the underlying R code or data changes. R Markdown
documents (`.Rmd` files) can be rendered to various formats including HTML and
PDF. The R code in an `.Rmd` document is processed by `knitr`, while the
resulting `.md` file is rendered by `pandoc` to the final output formats
(*e.g.* HTML or PDF). Historically, R Markdown is an extension of the older
`Sweave/Latex` environment. Rendering of mathematical expressions and reference
management is also supported by R Markdown using embedded Latex syntax and
Bibtex, respectively. A new and related publishing environemt is [Quarto](https://quarto.org/docs/tools/neovim.html) (not covered here).

## Quick Start

### Install R Markdown

To work with this tutorial, the `rmarkdown` package needs to be installed on a system.

``` r
install.packages("rmarkdown")
```

### Initialize a new R Markdown (`Rmd`) script

To minimize typing, it can be helful to start with an R Markdown template and
then modify it as needed. Note the file name of an R Markdown scirpt needs to
have the extension `.Rmd`. Template files for the following examples are available
here:

  - R Markdown sample script: [`sample.Rmd`](https://raw.githubusercontent.com/tgirke/GEN242/main/static/custom/rmarkdown/sample.Rmd)
  - Bibtex file for handling citations and reference section: [`bibtex.bib`](https://raw.githubusercontent.com/tgirke/GEN242/main/content/en/tutorials/rmarkdown/bibtex.bib)

Users want to download these files, open the `sample.Rmd` file with their preferred R IDE
(*e.g.* RStudio, vim or emacs), initilize an R session and then direct their R session to
the location of these two files.

### Metadata section

The metadata section (YAML header) in an R Markdown script defines how it will be processed and
rendered. The metadata section also includes both title, author, and date information as well as
options for customizing the output format. For instance, PDF and HTML output can be defined
with `pdf_document` and `html_document`, respectively. The `BiocStyle::` prefix will use the
formatting style of the [`BiocStyle`](http://bioconductor.org/packages/release/bioc/html/BiocStyle.html)
package from Bioconductor.

``` 
 ---
title: "My First R Markdown Document"
author: "Author: First Last"
date: "Last update: 17 May, 2023"
output:
  BiocStyle::html_document:
    toc: true
    toc_depth: 3
    fig_caption: yes

fontsize: 14pt
bibliography: bibtex.bib
 ---
```

### Render `Rmd` script

An R Markdown script can be evaluated and rendered with the following `render` command or by pressing the `knit` button in RStudio.
The `output_format` argument defines the format of the output (*e.g.* `html_document` or `pdf_document`). The setting `output_format="all"` will generate
all supported output formats. Alternatively, one can specify several output formats in the metadata section.

``` r
rmarkdown::render("sample.Rmd", clean=TRUE, output_format="BiocStyle::html_document")
```

The following shows two options how to run the rendering from the command-line. To render to PDF format, use the argument setting: `output_format="pdf_document"`.

``` sh
$ Rscript -e "rmarkdown::render('sample.Rmd', output_format='BiocStyle::html_document', clean=TRUE)"
```

Alternatively, one can use a Makefile to evaluate and render an R Markdown
script. A sample Makefile for rendering the above `sample.Rmd` can be
downloaded [`here`](https://raw.githubusercontent.com/tgirke/GEN242-2018/gh-pages/_vignettes/07_Rbasics/Makefile).
To apply it to a custom `Rmd` file, one needs open the Makefile in a text
editor and change the value assigned to `MAIN` (line 13) to the base name of
the corresponding `.Rmd` file (*e.g.* assign `systemPipeRNAseq` if the file
name is `systemPipeRNAseq.Rmd`). To execute the `Makefile`, run the following
command from the command-line.

``` sh
$ make -B
```

### R code chunks

R Code Chunks can be embedded in an R Markdown script by using three backticks
at the beginning of a new line along with arguments enclosed in curly braces
controlling the behavior of the code. The following lines contain the
plain R code. A code chunk is terminated by a new line starting with three backticks.
The following shows an example of such a code chunk. Note the backslashes are
not part of it. They have been added to print the code chunk syntax in this document.

```` 
	```\{r code_chunk_name, eval=FALSE\}
	x <- 1:10
	```
````

The following lists the most important arguments to control the behavior of R code chunks:

  - `r`: specifies language for code chunk, here R
  - `chode_chunk_name`: name of code chunk; this name needs to be unique within an Rmd
  - `eval`: if assigned `TRUE` the code will be evaluated
  - `warning`: if assigned `FALSE` warnings will not be shown
  - `message`: if assigned `FALSE` messages will not be shown
  - `cache`: if assigned `TRUE` results will be cached to reuse in future rendering instances
  - `fig.height`: allows to specify height of figures in inches
  - `fig.width`: allows to specify width of figures in inches

For more details on code chunk options see [here](https://www.rstudio.com/wp-content/uploads/2015/03/rmarkdown-reference.pdf).

### Learning Markdown

The basic syntax of Markdown and derivatives like kramdown is extremely easy to learn. Rather
than providing another introduction on this topic, here are some useful sites for learning Markdown:

  - [R Markdown Online Book](https://bookdown.org/yihui/rmarkdown/)
  - [Markdown Intro on GitHub](https://guides.github.com/features/mastering-markdown/)
  - [Markdown Cheet Sheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet)
  - [Markdown Basics from RStudio](http://rmarkdown.rstudio.com/authoring_basics.html)
  - [R Markdown Cheat Sheet](http://www.rstudio.com/wp-content/uploads/2015/02/rmarkdown-cheatsheet.pdf)
  - [kramdown Syntax](http://kramdown.gettalong.org/syntax.html)

### Tables

There are several ways to render tables. First, they can be printed within the R code chunks. Second,
much nicer formatted tables can be generated with the functions `kable`, `kableExtra`, `pander` or `xtable`. The following
example uses `kable` from the `knitr` package.

### With `knitr::kable`

``` r
library(knitr)
kable(iris[1:12,])
```

| Sepal.Length | Sepal.Width | Petal.Length | Petal.Width | Species |
| -----------: | ----------: | -----------: | ----------: | :------ |
|          5.1 |         3.5 |          1.4 |         0.2 | setosa  |
|          4.9 |         3.0 |          1.4 |         0.2 | setosa  |
|          4.7 |         3.2 |          1.3 |         0.2 | setosa  |
|          4.6 |         3.1 |          1.5 |         0.2 | setosa  |
|          5.0 |         3.6 |          1.4 |         0.2 | setosa  |
|          5.4 |         3.9 |          1.7 |         0.4 | setosa  |
|          4.6 |         3.4 |          1.4 |         0.3 | setosa  |
|          5.0 |         3.4 |          1.5 |         0.2 | setosa  |
|          4.4 |         2.9 |          1.4 |         0.2 | setosa  |
|          4.9 |         3.1 |          1.5 |         0.1 | setosa  |
|          5.4 |         3.7 |          1.5 |         0.2 | setosa  |
|          4.8 |         3.4 |          1.6 |         0.2 | setosa  |

A much more elegant and powerful solution is to create fully interactive tables with the [`DT` package](https://rstudio.github.io/DT/).
This JavaScirpt based environment provides a wrapper to the `DataTables` library using jQuery. The resulting tables can be sorted, queried and resized by the
user. Note, R Markdown source files containing JavaScript components can only be rendered into HTML and not PDF.

### With `DT::datatable`

``` r
library(DT)
datatable(iris)
```

<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-1" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-1">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120","121","122","123","124","125","126","127","128","129","130","131","132","133","134","135","136","137","138","139","140","141","142","143","144","145","146","147","148","149","150"],[5.1,4.9,4.7,4.6,5,5.4,4.6,5,4.4,4.9,5.4,4.8,4.8,4.3,5.8,5.7,5.4,5.1,5.7,5.1,5.4,5.1,4.6,5.1,4.8,5,5,5.2,5.2,4.7,4.8,5.4,5.2,5.5,4.9,5,5.5,4.9,4.4,5.1,5,4.5,4.4,5,5.1,4.8,5.1,4.6,5.3,5,7,6.4,6.9,5.5,6.5,5.7,6.3,4.9,6.6,5.2,5,5.9,6,6.1,5.6,6.7,5.6,5.8,6.2,5.6,5.9,6.1,6.3,6.1,6.4,6.6,6.8,6.7,6,5.7,5.5,5.5,5.8,6,5.4,6,6.7,6.3,5.6,5.5,5.5,6.1,5.8,5,5.6,5.7,5.7,6.2,5.1,5.7,6.3,5.8,7.1,6.3,6.5,7.6,4.9,7.3,6.7,7.2,6.5,6.4,6.8,5.7,5.8,6.4,6.5,7.7,7.7,6,6.9,5.6,7.7,6.3,6.7,7.2,6.2,6.1,6.4,7.2,7.4,7.9,6.4,6.3,6.1,7.7,6.3,6.4,6,6.9,6.7,6.9,5.8,6.8,6.7,6.7,6.3,6.5,6.2,5.9],[3.5,3,3.2,3.1,3.6,3.9,3.4,3.4,2.9,3.1,3.7,3.4,3,3,4,4.4,3.9,3.5,3.8,3.8,3.4,3.7,3.6,3.3,3.4,3,3.4,3.5,3.4,3.2,3.1,3.4,4.1,4.2,3.1,3.2,3.5,3.6,3,3.4,3.5,2.3,3.2,3.5,3.8,3,3.8,3.2,3.7,3.3,3.2,3.2,3.1,2.3,2.8,2.8,3.3,2.4,2.9,2.7,2,3,2.2,2.9,2.9,3.1,3,2.7,2.2,2.5,3.2,2.8,2.5,2.8,2.9,3,2.8,3,2.9,2.6,2.4,2.4,2.7,2.7,3,3.4,3.1,2.3,3,2.5,2.6,3,2.6,2.3,2.7,3,2.9,2.9,2.5,2.8,3.3,2.7,3,2.9,3,3,2.5,2.9,2.5,3.6,3.2,2.7,3,2.5,2.8,3.2,3,3.8,2.6,2.2,3.2,2.8,2.8,2.7,3.3,3.2,2.8,3,2.8,3,2.8,3.8,2.8,2.8,2.6,3,3.4,3.1,3,3.1,3.1,3.1,2.7,3.2,3.3,3,2.5,3,3.4,3],[1.4,1.4,1.3,1.5,1.4,1.7,1.4,1.5,1.4,1.5,1.5,1.6,1.4,1.1,1.2,1.5,1.3,1.4,1.7,1.5,1.7,1.5,1,1.7,1.9,1.6,1.6,1.5,1.4,1.6,1.6,1.5,1.5,1.4,1.5,1.2,1.3,1.4,1.3,1.5,1.3,1.3,1.3,1.6,1.9,1.4,1.6,1.4,1.5,1.4,4.7,4.5,4.9,4,4.6,4.5,4.7,3.3,4.6,3.9,3.5,4.2,4,4.7,3.6,4.4,4.5,4.1,4.5,3.9,4.8,4,4.9,4.7,4.3,4.4,4.8,5,4.5,3.5,3.8,3.7,3.9,5.1,4.5,4.5,4.7,4.4,4.1,4,4.4,4.6,4,3.3,4.2,4.2,4.2,4.3,3,4.1,6,5.1,5.9,5.6,5.8,6.6,4.5,6.3,5.8,6.1,5.1,5.3,5.5,5,5.1,5.3,5.5,6.7,6.9,5,5.7,4.9,6.7,4.9,5.7,6,4.8,4.9,5.6,5.8,6.1,6.4,5.6,5.1,5.6,6.1,5.6,5.5,4.8,5.4,5.6,5.1,5.1,5.9,5.7,5.2,5,5.2,5.4,5.1],[0.2,0.2,0.2,0.2,0.2,0.4,0.3,0.2,0.2,0.1,0.2,0.2,0.1,0.1,0.2,0.4,0.4,0.3,0.3,0.3,0.2,0.4,0.2,0.5,0.2,0.2,0.4,0.2,0.2,0.2,0.2,0.4,0.1,0.2,0.2,0.2,0.2,0.1,0.2,0.2,0.3,0.3,0.2,0.6,0.4,0.3,0.2,0.2,0.2,0.2,1.4,1.5,1.5,1.3,1.5,1.3,1.6,1,1.3,1.4,1,1.5,1,1.4,1.3,1.4,1.5,1,1.5,1.1,1.8,1.3,1.5,1.2,1.3,1.4,1.4,1.7,1.5,1,1.1,1,1.2,1.6,1.5,1.6,1.5,1.3,1.3,1.3,1.2,1.4,1.2,1,1.3,1.2,1.3,1.3,1.1,1.3,2.5,1.9,2.1,1.8,2.2,2.1,1.7,1.8,1.8,2.5,2,1.9,2.1,2,2.4,2.3,1.8,2.2,2.3,1.5,2.3,2,2,1.8,2.1,1.8,1.8,1.8,2.1,1.6,1.9,2,2.2,1.5,1.4,2.3,2.4,1.8,1.8,2.1,2.4,2.3,1.9,2.3,2.5,2.3,1.9,2,2.3,1.8],["setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Sepal.Length<\/th>\n      <th>Sepal.Width<\/th>\n      <th>Petal.Length<\/th>\n      <th>Petal.Width<\/th>\n      <th>Species<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

### Figures

Plots generated by the R code chunks in an R Markdown document can be automatically
inserted in the output file. The size of the figure can be controlled with the `fig.height`
and `fig.width` arguments.

``` r
library(ggplot2)
dsmall <- diamonds[sample(nrow(diamonds), 1000), ]
ggplot(dsmall, aes(color, price/carat)) + geom_jitter(alpha = I(1 / 2), aes(color=color))
```

<img src="/en/tutorials/rmarkdown/rmarkdown_files/figure-html/some_jitter_plot-1.png" width="672" />

Sometimes it can be useful to explicitly write an image to a file and then insert that
image into the final document by referencing its file name in the R Markdown source. For
instance, this can be useful for time consuming analyses. The following code will generate a
file named `myplot.png`. To insert the file in the final document, one can use standard
Markdown or HTML syntax, *e.g.*: `<img src="myplot.png"/>`.

``` r
png("myplot.png")
ggplot(dsmall, aes(color, price/carat)) + geom_jitter(alpha = I(1 / 2), aes(color=color))
dev.off()
```

    ## png 
    ##   2

<center>

<img title="some_title" src="../images/myplot.png"/>

</center>

### Inline R code

To evaluate R code inline, one can enclose an R expression with a single back-tick
followed by `r` and then the actual expression. For instance, the back-ticked version
of ‘r 1 + 1’ evaluates to 2 and ‘r pi’ evaluates to 3.1415927.

### Mathematical equations

To render mathematical equations, one can use standard Latex syntax. When expressions are
enclosed with single `$` signs then they will be shown inline, while
enclosing them with double `$$` signs will show them in display mode. For instance, the following
Latex syntax `d(X,Y) = \sqrt[]{ \sum_{i=1}^{n}{(x_{i}-y_{i})^2} }` renders in display mode as follows:

`$$d(X,Y) = \sqrt[]{ \sum_{i=1}^{n}{(x_{i}-y_{i})^2} }$$`

To learn LaTeX syntax for mathematical equations, one can consult various online manuals, such as
this [Wikibooks tutorial](https://en.wikibooks.org/wiki/LaTeX/Mathematics), or use an online
equation rendering and checking tool, such as this [one](https://arachnoid.com/latex/).

### Citations and bibliographies

Citations and bibliographies can be autogenerated in R Markdown in a similar
way as in Latex/Bibtex. Reference collections should be stored in a separate
file in Bibtex or other supported formats. To cite a publication in an R Markdown
script, one uses the syntax `[@<id1>]` where `<id1>` needs to be replaced with a
reference identifier present in the Bibtex database listed in the metadata section
of the R Markdown script (*e.g.* `bibtex.bib`). For instance, to cite Lawrence et al. 
(2013), one uses its reference identifier (*e.g.* `Lawrence2013-kt`) as `<id1>` (Lawrence et al. 2013).
This will place the citation inline in the text and add the corresponding
reference to a reference list at the end of the output document. For the latter a
special section called `References` needs to be specified at the end of the R Markdown script.
To fine control the formatting of citations and reference lists, users want to consult this
[R Markdown page](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html).
Also, for general reference management and obtaining references in Bibtex format [Paperpile](https://paperpile.com/features)
can be very helpful.

### Viewing R Markdown report on HPCC cluster

R Markdown reports located on UCR’s HPCC Cluster can be viewed locally in a web browser (without moving
the source HTML) by creating a symbolic link from a user’s `.html` directory. This way any updates to
the report will show up immediately without creating another copy of the HTML file. For instance, if user
`ttest` has generated an R Markdown report under `~/bigdata/today/rmarkdown/sample.html`, then the
symbolic link can be created as follows:

``` r
cd ~/.html
ln -s ~/bigdata/today/rmarkdown/sample.html sample.html
```

After this one can view the report in a web browser using this URL <https://cluster.hpcc.ucr.edu/~ttest/rmarkdown/sample.html>.
If necessary access to the URL can be restricted with a password following the instructions [here](http://hpcc.ucr.edu/manuals_linux-cluster_sharing.html#sharing-files-on-the-web).

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
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] ggplot2_3.4.2 DT_0.27       knitr_1.42   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.3      jsonlite_1.8.4    highr_0.10        dplyr_1.1.2      
    ##  [5] compiler_4.3.0    tidyselect_1.2.0  jquerylib_0.1.4   scales_1.2.1     
    ##  [9] yaml_2.3.7        fastmap_1.1.1     R6_2.5.1          labeling_0.4.2   
    ## [13] generics_0.1.3    htmlwidgets_1.6.2 tibble_3.2.1      bookdown_0.33    
    ## [17] munsell_0.5.0     bslib_0.4.2       pillar_1.9.0      rlang_1.1.1      
    ## [21] utf8_1.2.3        cachem_1.0.8      xfun_0.39         sass_0.4.6       
    ## [25] viridisLite_0.4.2 cli_3.6.1         withr_2.5.0       magrittr_2.0.3   
    ## [29] crosstalk_1.2.0   digest_0.6.31     grid_4.3.0        lifecycle_1.0.3  
    ## [33] vctrs_0.6.2       evaluate_0.21     glue_1.6.2        farver_2.1.1     
    ## [37] blogdown_1.16     fansi_1.0.4       colorspace_2.1-0  rmarkdown_2.21   
    ## [41] tools_4.3.0       pkgconfig_2.0.3   ellipsis_0.3.2    htmltools_0.5.5

## References

<div id="refs" class="references hanging-indent">

<div id="ref-Lawrence2013-kt">

Lawrence, Michael, Wolfgang Huber, Hervé Pagès, Patrick Aboyoun, Marc Carlson, Robert Gentleman, Martin T Morgan, and Vincent J Carey. 2013. “Software for Computing and Annotating Genomic Ranges.” *PLoS Comput. Biol.* 9 (8): e1003118. <https://doi.org/10.1371/journal.pcbi.1003118>.

</div>

</div>
