---
title: Overview of Course Projects 
linkTitle: "Overview of Projects"
description: >
type: docs
weight: 401
---

<br></br>

## Introduction

During the tutorial sessions of this class all students will perform the basic
data analysis of at least two NGS Workflows including RNA-Seq and VAR-Seq.
In addition, every student will work on a Challenge
Project addressing a specific data analysis task within one of the general NGS
Workflows. Students will also present a scientific paper closely related to
their challenge topic (see
[here](https://girke.bioinformatics.ucr.edu/GEN242/assignments/presentations/paper_presentations/)).
To facilitate teamwork and communication with instructors, each course project will be
assigned a private GitHub repository.

The results of the Challenge Projects will be presented by each student
during the last week of the course (see Slideshow Template
[here](https://bit.ly/3oMz9gb)).
In addition, each student will write a detailed analysis report for the assigned
course project. This report needs to include all analysis steps of the
corresponding NGS Workflow (_e.g._ full RNA-Seq analysis) as well as the
code and results of the Challenge Project. The final project reports will be written
in R Markdown. A basic tutorial on R Markdown is available [here](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/rmarkdown/rmarkdown/). 
Both the R Markdown script (`.Rmd`) along with the rendered HTML or PDF report will 
be submitted to each student's private project GitHub repository. All helper code used for 
the challenge project needs to be organized in well documented R functions of each 
project's `*_Fct.R` script. The custom functions defined in `*_Fct.R` need to be imported (sourced)
and used in the main Rmd project report. Other scripts used by the challenge projects need to be called from the `*_Fct.R` (_e.g._ via R's [system function](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/rprogramming/rprogramming/#calling-external-software)) and also uploaded to the project repos. The expected structure of the final project report is outlined below. 

The reports should be submitted to each student’s private project GitHub repository. For 
the report each student should create in this repository a new directory named after their
workflow project and include in it the following files: 

* `.Rmd` source script of project report 
* Report rendered from `.Rmd` source in HTML or PDF format
* `._Fct.R` file containing all helper functions written for challenge project
* __Submission Deadline__ for reports: 6:00 PM, June 11th, 2024


## Structure of final project report

1. Abstract
2. Introduction
3. Methods
    + Short description of methods used by NGS workflow
    + Detailed description of methods used for challenge project
4. Results and Discussion
    + Includes all components of NGS workflow as well as challenge project
5. Conclusions
6. Acknowledgments
7. References
8. Supplement (optional)


