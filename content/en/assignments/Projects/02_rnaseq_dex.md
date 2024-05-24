---
title: "RNA-Seq - Differentially Expressed Transcript (DET) Analysis"
linkTitle: "DEG Methods"
description: >
type: docs
weight: 405
---

<br></br>

## RNA-Seq Workflow  

1. Read quality assessment, filtering and trimming 
2. Map reads against reference genome 
3. Perform read counting for required ranges (_e.g._ exonic gene ranges)
4. Normalization of read counts
5. Identification of differentially expressed genes (DEGs)
6. Clustering of gene expression profiles 
7. Gene set enrichment analysis

## Challenge Projects

### Analysis of Differentially Expressed Exons and Transcripts

+ Run the workflow from start to finish (steps 1-7) on the full RNA-Seq data set from Howard et al. (2013).
+ Challenge project tasks
    + Group 1: Perform differential exon analysis with [DEXseq](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html). Assess the results as follows:
        + Identify genes that show differential exon usage according to DEXseq. Optionally, perform functional gene set enrichment analysis on the obained gene set.  
        + Compare the results with the findings of the splice variant analysis reported by Howard et al (2013).
        + Optional: compare the performance of DEXseq and Kallisto/Sleuth (see below) with the results from the Howard et al (2013) paper in the form of ROC plots. As 'pseudo' ground truth the consensus DET set or similar could be used.
    + Group 2: Same as above but with [Kallisto](https://pachterlab.github.io/kallisto/download.html)/[Sleuth](https://github.com/pachterlab/sleuth) (Pimentel et al, 2017) or [DTUrtle](https://tobitekath.github.io/DTUrtle/) (Tekath and Dugas, 2021).

## References

+ Anders S, Reyes A, Huber W (2012) Detecting differential usage of exons from RNA-seq data. Genome Res 22: 2008–2017 [PubMed](https://pubmed.ncbi.nlm.nih.gov/22722343/)
+ Howard, B.E. et al., 2013. High-throughput RNA sequencing of pseudomonas-infected Arabidopsis reveals hidden transcriptome complexity and novel splice variants. PloS one, 8(10), p.e74183. [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/24098335)
+ Guo Y, Li C-I, Ye F, Shyr Y (2013) Evaluation of read count based RNAseq analysis methods. BMC Genomics 14 Suppl 8: S2 [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/24564449)
+ Hardcastle TJ, Kelly KA (2010) baySeq: empirical Bayesian methods for identifying differential expression in sequence count data. BMC Bioinformatics 11: 422 [PubMed](https://pubmed.ncbi.nlm.nih.gov/20698981/)
+ Liu R, Holik AZ, Su S, Jansz N, Chen K, Leong HS, Blewitt ME, Asselin-Labat M-L, Smyth GK, Ritchie ME (2015) Why weight? Modelling sample and observational level variability improves power in RNA-seq analyses. Nucleic Acids Res. doi: 10.1093/nar/gkv412. [PubMed](https://pubmed.ncbi.nlm.nih.gov/25925576/)
+ Love MI, Huber W, Anders S (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol 15: 550 [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/25516281)
+ Pimentel H, Bray NL, Puente S, Melsted P, Pachter L (2017) Differential analysis of RNA-seq incorporating quantification uncertainty. Nat Methods 14: 687–690. [PubMed](https://pubmed.ncbi.nlm.nih.gov/28581496/)
+ Tekath T, Dugas M (2021) Differential transcript usage analysis of bulk and single-cell RNA-seq data with DTUrtle. Bioinformatics 37: 3781–3787. [PubMed](https://pubmed.ncbi.nlm.nih.gov/34469510/)
+ Zhou X, Lindsay H, Robinson MD (2014) Robustly detecting differential expression in RNA sequencing data using observation weights. Nucleic Acids Res 42: e91 [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/24753412)




