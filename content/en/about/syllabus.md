---
title: "Syllabus - GEN242"
linkTitle: "Syllabus"
type: docs
description: >
weight: 2
---

## Course title

Data Analysis in Genome Biology <br>
GEN242 - Spring 2024

## Printable syllabus

See Google Doc version [here](https://bit.ly/42UEIbB).

## Instructor
Name: Thomas Girke <br> 
Email: thomas.girke@ucr.edu <br>
Office location: virtual via Zoom <br>
Office hour: Tue 4:30 - 5:30 PM & Fri 4:00 - 5:00 PM <br>
Zoom URL: privately shared

## Description

Introduction to algorithms, statistical methods and data analysis programming
routines relevant for genome biology. The class consists of three
main components: lectures, hands-on practicals and student course projects. The
lecture topics cover databases, sequence (NGS) analysis, phylogenetics,
comparative genomics, genome-wide profiling methods, network biology and more.
The hands-on practicals include homework assignments and course projects
focusing on data analysis programming of next generation genome data using
command-line tools on a computer cluster and the programming environment R.
Credit: 4 units (2x 1.5 hours lectures, 1 hour discussion)

## Objectives of course

- Acquire understanding of algorithms used in bioinformatics
- Obtain hands-on experience in large scale data analysis.

## Prerequisites
The main prerequisite for this course is a strong interest in acquiring the
skills required for mastering the computational aspects of modern genome
research.

## Structure of course
Two lectures per week (1.5 hours each) plus one discussion section (1 hour).
During the first weeks the discussion section will be used for data analysis
tutorials using Linux command-line tools and R. 

## Delivery mode
This class will be instructed online via Zoom in synchronous mode. Upon
request, office hours can be attended in person (hybrid mode).


## Time
Lecture: Tue/Thu 2:00-3:20 PM <br>
Discussion: Thu 3:30-4:20 PM

## Location
Online via video conferencing software

## Grading
1. Homework assignments: 40%
2. Scientific paper presentation: 20%
3. Course project presentations: 20%
4. Final project report: 20%

Additional details about the grading system are provided in this [table](https://bit.ly/3udPPMA) (see both tabs).

__Grading policy:__ Given the diverse educational background of the students in GEN242, all assignments are designed to be solvable by students from both experimental and quantitative disciplines, including those with no or only limited prior experience in programming and/or data modeling. The weight of each of the four gradable components in this class is given above in percent. 
__(1)__ The homeworks include 8-10 assignments throughout the class. They cover algorithms and data analysis programming problems using the R language. The grading of these assignments is mainly based on correctness, reproducibility and reusability of the analysis code. 
__(2-4)__ Students will work on a Challenge Project (individually or in group) addressing a specific data analysis problem in genome data sciences. As part of their project, students will present a scientific paper __(2)__ closely related to their project (see reading list for details). The results of the Challenge Projects __(3)__ will be presented and discussed by each student at the end of the course. In addition, each student will write a detailed analysis report __(4)__ of the assigned course project. The latter will be written in the style of a scientific publication and should include a detailed description of the results including all analysis code to fully reproduce the project results followed by a critical discussion of the outcomes. The grading of both the paper and project presentations __(2-3)__ includes anonymous feedback from all students as well as the instructor, where understanding of the material, clarity of the oral presentations and critical thinking are the main grading criteria. The final project reports __(4)__ will be graded by the instructor with an emphasis on scientific and coding accuracy, overall understanding of the topic, as well as reproducibility of the results.

## Materials needed
Students are expected to bring to each class meeting a laptop with a functional wireless
connection and a recent internet browser version (e.g. Firefox, Chrome or
Safari) preinstalled. Tablet computers with mobile operating systems are not
suitable for running the required software. User accounts on a research
computer cluster will be provided at the beginning of the course. To log in to
the cluster, students also need to install a terminal application for their
operating system (_e.g._ [iTerm2](http://www.iterm2.com/#/section/home) on OS X,
and [PuTTY](http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html) or 
[MobaXterm](http://mobaxterm.mobatek.net/) on Windows) as well as a file exchange software such as
[FileZilla](https://filezilla-project.org/download.php?type=client). In
addition, a recent version of [R](http://www.r-project.org) and
[RStudio](http://rstudio.org) should be installed.

If possible students may want to attend class sessions from a monitor setup
with either one large monitor (wide enough to display several windows) or two
separate monitors. This allows simultaneous viewing of presentations on one
screen and following along hands-on practicals on the other screen. 

## Schedule


|Week     |Topic                                                                |
|---------|---------------------------------------------------------------------|
| Week 1  | Course Introduction                                                 |
|         | Databases and Software for Genome Biology                           |
|         | Discussion: Introduction to Linux and HPC                           |
|         | Reading: A1, T1, T2                                                 |
| Week 2  | Sequencing Technologies                                             |
|         | Discussion: Introduction to R                                       |
|         | Reading: A2-A4, T3                                                  |
| Week 3  | Sequence Alignments and Searching                                   |
|         | Multiple Sequence Alignments                                        |
|         | Discussion: Programming in R and Parallel Evaluations               |
|         | Reading: A5-A6, T4-T5                                               |
| Week 4  | Short Read Alignment Algorithms                                     |
|         | Discussion: Basics of NGS Analysis                                  |
|         | Reading: A7-A10, T6                                                 |
| Week 5  | Gene Expression Analysis, Microarrays, bulk RNA-Seq and scRNA-Seq   |
|         | Discussion: NGS Workflow Overview; RNA-Seq Analysis                 |
|         | Reading: A11-A15, T7-T8                                             |
| Week 6  | Analysis of ChIP-Seq Experiments                                    |
|         | Discussion: ChIP-Seq Analysis                                       |
|         | Reading: A16-A18, T9-T10                                            |
| Week 7  | Students present publication related to their chosen course project |
|         | Discussion: Q&A about papers                                        |
|         | Reading: A19-A23                                                    |
| Week 8  | Clustering algorithms                                               |
|         | Pathway and GO annotation systems                                   |
|         | Discussion: Gene Set Enrichment Analysis                            |
|         | Reading: A24-A26, T7 (Sec 3.14-3.15), T11                           |
| Week 9  | Genome and Transcriptome Assembly Algorithms                        | 
|         | Profile HMMs for Protein Family Modeling                            |
|         | Introduction to Phylogenetics                                       |
|         | Discussion: Graphics and Data Visualization                         |
|         | Reading: A27-A29, T12                                               |
| Week 10 | Final presentations of student data analysis projects               |
|         | Discussion: Tips and tricks for efficient data analysis programming |
|         | Reading: A30-A31, T3 (Sec 12,13-17)                                 |


## Reading list


### Journal articles


A1. Huber W, Carey VJ, Gentleman R, Anders S, Carlson M, Carvalho BS, Bravo HC, Davis S, Gatto L, Girke T, et al (2015) Orchestrating high-throughput genomic analysis with Bioconductor. Nat Methods 12: 115–121

A2. Metzker, M. L., Jan 2010. Sequencing technologies - the next generation. Nat Rev Genet 11 (1), 31–46. 

A3. Needleman SB, Wunsch CD (1970) A general method applicable to the search for similarities in the amino acid sequence of two proteins. J Mol Biol 48, 443-453. 

A4. Smith TF, Waterman MS (1981) Identification of common molecular subsequences. J Mol Biol 147, 195-197. 

A5. Corpet F (1988) Multiple sequence alignment with hierarchical clustering. Nucleic Acids Res 16, 10881-90. 

A6. Altschul, S. F., Gish, W., Miller, W., Myers, E. W., Lipman, D. J., Oct 1990. Basic local alignment search tool. J Mol Biol 215 (3), 403–410.

A7. Li, H, Durbin, R (2009) Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics, 25: 1754-1760.

A8. Dobin, A., Davis, C.A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M., Gingeras, T.R., 2012. STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15–21. 

A9. Langmead, B, Salzberg, S L (2012) Fast gapped-read alignment with Bowtie 2. Nat Methods, 9: 357-359. 

A10. Kim D, Langmead B, Salzberg SL (2015) HISAT: a fast spliced aligner with low memory requirements. Nat Methods 12: 357–360

A11. Bray NL, Pimentel H, Melsted P, Pachter L (2016) Near-optimal probabilistic RNA-seq quantification. Nat Biotechnol. doi: 10.1038/nbt.3519

A12. Love MI, Huber W, Anders S (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol 15: 550

A13. Zhou X, Lindsay H, Robinson MD (2014) Robustly detecting differential expression in RNA sequencing data using observation weights. Nucleic Acids Res 42: e91

A14. Anders, S, Reyes, A, Huber, W (2012) Detecting differential usage of exons from RNA-seq data. Genome Res, 22: 2008-2017.

A15. Soneson, C, Delorenzi, M (2013) A comparison of methods for differential expression analysis of RNA-seq data. BMC Bioinformatics, 14: 91-91.

A16. Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, Nussbaum C, Myers RM, Brown M, Li W, et al (2008) Model-based analysis of ChIP-Seq (MACS). Genome Biol. doi: 10.1186/gb-2008-9-9-r137

A17. Wilbanks EG, Facciotti MT (2010) Evaluation of algorithm performance in ChIP-seq peak detection. PLoS One. doi: 10.1371/journal.pone.0011471.

A18. Landt et al. (2012) ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia. Genome Res, 22: 1813-1831. 

A19. McLeay, Robert C, and Timothy L Bailey. 2010. “Motif Enrichment Analysis: A Unified Framework and an Evaluation on ChIP Data.” BMC Bioinformatics 11: 165.

A20. Machanick, P, Bailey, T L (2011) MEME-ChIP: motif analysis of large DNA datasets. Bioinformatics, 27: 1696-1697.

A21. Tompa, M, N Li, T L Bailey, G M Church, B De Moor, E Eskin, A V Favorov, et al. 2005. “Assessing Computational Tools for the Discovery of Transcription Factor Binding Sites.” Nature Biotechnology 23 (1): 137–44.

A22. DePristo MA, Banks E, Poplin R, Garimella KV, Maguire JR, Hartl C, Philippakis AA, del Angel G, Rivas MA, Hanna M, et al (2011) A framework for variation discovery and genotyping using next-generation DNA sequencing data. Nat Genet 43: 491–498.

A23. Shihab HA, Rogers MF, Gough J, Mort M, Cooper DN, Day INM, Gaunt TR, Campbell C (2015) An integrative approach to predicting the functional effects of non-coding and coding sequence variation. Bioinformatics 31: 1536–1543.

A24. Raymond JW, Blankley CJ, Willett P (2003) Comparison of chemical clustering methods using graph- and fingerprint-based similarity measures. J Mol Graph Model 21: 421–433.

A25. Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, Paulovich A, Pomeroy SL, Golub TR, Lander ES, et al (2005) Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc Natl Acad Sci U S A 102: 15545–15550.

A26. Ashburner M, Ball CA, Blake JA, Botstein D, Butler H, Cherry JM, Davis AP, Dolinski K, Dwight SS, Eppig JT, et al (2000) Gene ontology: tool for the unification of biology. The Gene Ontology Consortium. Nat Genet 25: 25–29.

A27. Zyla J, Marczyk M, Domaszewska T, Kaufmann SHE, Polanska J, Weiner J (2019) Gene set enrichment for reproducible science: comparison of CERNO and eight other algorithms. Bioinformatics 35: 5146–5154

A28. Alkan, C, Sajjadian, S, Eichler, E E (2011) Limitations of next-generation genome sequence assembly. Nat Methods, 8: 61-65. 

A29. Eddy SR (1998) Profile hidden Markov models. Bioinformatics 14: 755–763.

A30. Grabherr, M G, Haas, B J, Yassour, M, Levin, J Z, Thompson, D A, Amit, I, Adiconis, X, Fan, L, Raychowdhury, R, Zeng, Q, Chen, Z, Mauceli, E, Hacohen, N, Gnirke, A, Rhind, N, di Palma, F, Birren, B W, Nusbaum, C, Lindblad-Toh, K, Friedman, N, Regev, A (2011) Full-length transcriptome assembly from RNA-Seq data without a reference genome. Nat Biotechnol, 29: 644-652. 

A31. Zeitouni, B, Boeva, V, Janoueix-Lerosey, I, Loeillet, S, Legoix-ne, P, Nicolas, A, Delattre, O, Barillot, E (2010) SVDetect: a tool to identify genomic structural variations from paired-end and mate-pair sequencing data. Bioinformatics, 26: 1895-1896.

A32. Ronquist F, Teslenko M, van der Mark P, Ayres DL, Darling A, Höhna S, Larget B, Liu L, Suchard MA, Huelsenbeck JP (2012) MrBayes 3.2: efficient Bayesian phylogenetic inference and model choice across a large model space. Syst Biol 61: 539–542.

A33. Law CW, Zeglinski K, Dong X, Alhamdoosh M, Smyth GK, Ritchie ME (2020) A guide to creating design matrices for gene expression experiments. F1000Res 9: 1444

A34. Langfelder P, Luo R, Oldham MC, Horvath S (2011) Is my network module preserved and reproducible? PLoS Comput Biol 7: e1001057

A35. McInnes L, Healy J, Melville J (2018) UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. arXiv [stat.ML]

A36. Rodriguez MZ, Comin CH, Casanova D, Bruno OM, Amancio DR, Costa L da F, Rodrigues FA (2019) Clustering algorithms: A comparative approach. PLoS One 14: e0210236

A37. Shulse CN, Cole BJ, Ciobanu D, Lin J, Yoshinaga Y, Gouran M, Turco GM, Zhu Y, O’Malley RC, Brady SM, et al (2019) High-Throughput Single-Cell Transcriptome Profiling of Plant Cell Types. Cell Rep 27: 2241–2247.e4

A38. Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, Paulovich A, Pomeroy SL, Golub TR, Lander ES, et al (2005) Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc Natl Acad Sci U S A 102: 15545–15550

A39. Sun S, Zhu J, Ma Y, Zhou X (2019) Accuracy, robustness and scalability of dimensionality reduction methods for single-cell RNA-seq analysis. Genome Biol 20: 269

### Tutorials

T1. [GitHub Introduction](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/github/github/)

T2. [Introduction to Computer Clusters and Linux](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/linux/linux/)

T3. [Introduction to R](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/rbasics/rbasics/) 

T4. [Programming in R](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/rprogramming/rprogramming/)

T5. [Parallel R](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/rparallel/rparallel/)

T6. [NGS Analysis Basics](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/rsequences/rsequences/)

T7. [NGS Workflows](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/systempiper/systempiper/)
 
T8. [RNA-Seq Workflow](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/sprnaseq/sprnaseq/) 

T9. [scRNA-Seq Embedding Methods](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/scrnaseq/scrnaseq/)

T10. [ChIP-Seq Workflow](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/spchipseq/spchipseq/)

T11. [R Markdown](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/rmarkdown/rmarkdown/)

T12. [Functional Enrichment Analysis](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/rfea/rfea/) 

T13. [Clustering and Network Analysis](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/rclustering/rclustering/)

T14. [Project Data](https://girke.bioinformatics.ucr.edu/GEN242/assignments/projects/project_data/)

T15. [Data Visualization](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/rgraphics/rgraphics/)

T16. [Shiny Apps](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/shinyapps/shinyapps/)

T17. [Building R Packages](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/rpackages/rpackages/)

T18. [dplyr, tidyr and some SQLite](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/dplyr/dplyr/)

T19. [Advanced: Common Workflow Language (CWL)](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/cmdtocwl/cmdtocwl/)


### Books 

Note: there is no need to purchase any books for this course as most reading material will be based on journal articles!

General
Jonathan Pevsner (2009) Bioinformatics and Functional Genomics. Wiley-Blackwell; 2nd Edition, 992 pages.

Algorithms
Jones N and Pevzner P (2004) An Introduction to Bioinformatics Algorithms. MIT Press, Massachusetts, 435 pages.

Sequence Analysis
Durbin, R, Eddy, S, Krogh, A, Mitchison, G. (1998) Biological Sequence Analysis: Probabilistic Models of Proteins and Nucleic Acids. Cambridge University Press, UK, 356 pages.

Parida L (2008) Pattern Discovery in Bioinformatics: Theory & Algorithms. CRC Press, London, 526 pages.

Profiling Bioinformatics
Gentleman, R, Carey, V, Dudoit, S, Irizarry, R, Huber, W (2005) Bioinformatics and Computational Biology Solutions Using R and Bioconductor. Springer, New York, 473 pages.

Phylogenetics 
Felsenstein, J (2004) Inferring Phylogenies. Sinauer, Massachusetts, 664 pages.

Paradis (2006) Analysis of Phylogenetics and Evolution with R. Springer, New York, 211 pages.

