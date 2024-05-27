---
title: "Compare Performance of Variant Callers"
linkTitle: "Variant Callers"
description: >
type: docs
weight: 410
---

<br></br>

## VAR-Seq Workflow  

1. Read preprocessing: filtering, quality trimming
2. Alignments
3. Alignment statistics
4. Variant calling: focus of challenge project
5. Variant filtering
6. Variant annotation
7. Combine results from many samples
8. Summary statistics of samples


## Challenge Project: Performance Comparisons of Variant Callers

+ Run the workflow from start to finish (steps 1-8) on the VAR-Seq data set from 
on the data set from Lu _et al_ (2012).
+ Challenge project tasks
    + Compare the performance of at least 2 variant callers, e.g. GATK, BCFtools, Octopus and DeepVariant. Include in your comparisons the following analysis/visualization steps
        + Report unique and common variants identified by tested variant callers.
        + Compare the results from (a) with the variants identified by Lu et al, 2012
        + Plot results from a-b as venn diagrams or similar (_e.g._ upset plots)
        + If there is enough time and interest, plot the performance of the variant callers in the form of ROC plots and calculate AUC values. As pseudo ground truth, one can either use the published variants or the union of the variants identified by all methods. 

## References

+ Barbitoff YA, Abasov R, Tvorogova VE, Glotov AS, Predeus AV (2022) Systematic benchmark of state-of-the-art variant calling pipelines identifies major factors affecting accuracy of coding sequence variant discovery. BMC Genomics 23: 155. [PubMed](https://pubmed.ncbi.nlm.nih.gov/35193511/)
+ Cooke DP, Wedge DC, Lunter G (2021) A unified haplotype-based method for accurate and comprehensive variant calling. Nat Biotechnol 39: 885–892. [PubMed](https://pubmed.ncbi.nlm.nih.gov/33782612/)
+ DePristo MA, Banks E, Poplin R, Garimella KV, Maguire JR, Hartl C, Philippakis AA, del Angel G, Rivas MA, Hanna M, et al (2011) A framework for variation discovery and genotyping using next-generation DNA sequencing data. Nat Genet 43: 491–498. [PubMed](https://pubmed.ncbi.nlm.nih.gov/21478889/) 
+ Li H (2011) A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics 27: 2987–2993. [PubMed](https://pubmed.ncbi.nlm.nih.gov/21903627/)
+ Lu P, Han X, Qi J, Yang J, Wijeratne AJ, Li T, Ma H (2012) Analysis of Arabidopsis genome-wide variations before and after meiosis and meiotic recombination by resequencing Landsberg erecta and all four products of a single meiosis. Genome Res 22: 508–518. [PubMed](https://pubmed.ncbi.nlm.nih.gov/22106370/)
+ Poplin R, Chang P-C, Alexander D, Schwartz S, Colthurst T, Ku A, Newburger D, Dijamco J, Nguyen N, Afshar PT, et al (2018) A universal SNP and small-indel variant caller using deep neural networks. Nat Biotechnol 36: 983–987. [PubMed](https://pubmed.ncbi.nlm.nih.gov/30247488/)





