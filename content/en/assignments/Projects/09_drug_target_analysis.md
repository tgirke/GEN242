---
title: "Drug-target analysis"
linkTitle: "Drug-target analysis"
description: >
type: docs
weight: 409
---

<br></br>

## ChIP-Seq Workflow  

1. Read quality assessment, filtering and trimming
2. Align reads to reference genome
3. Compute read coverage across genome
4. Peak calling with different methods and consensus peak identification
5. Annotate peaks
6. Differential binding analysis
7. Gene set enrichment analysis
8. Motif prediction to identify putative TF binding sites

## Challenge Project: Functional enrichment analysis (FEA)

+ Run workflow from start to finish (steps 1-8) on ChIP-Seq data set from Kaufman et al. (2010)
+ Challenge project tasks
    + Identify protein coding genes in peak regions 
    + Identify corresponding human orthologs 
    + Perform drug-target annotation analysis, e.g. with [drugTargetInteractions](https://bioconductor.org/packages/release/bioc/html/drugTargetInteractions.html) package
    + Identify similar drugs with two different structural similarity search algorithms (e.g. 2 fingerprint methods) 
    + Challenge question: which of the two structural similarity search tools identifies more similar small molecules that have annotated protein targets 
      in ChEMBL (DrugBank). Explore options on how to visualize the performance results. 

## References

+ Chen X, Reynolds CH (2002) Performance of similarity measures in 2D fragment-based similarity searching: comparison of structural descriptors and similarity coefficients. J Chem Inf Comput Sci 42: 1407–1414 [PubMed](https://pubmed.ncbi.nlm.nih.gov/12444738/)
+ Kaufmann, K, F Wellmer, J M Muiño, T Ferrier, S E Wuest, V Kumar, A Serrano-Mislata, et al. 2010. “Orchestration of Floral Initiation by APETALA1.” Science 328 (5974): 85–89. [PubMed](https://pubmed.ncbi.nlm.nih.gov/20360106/)



