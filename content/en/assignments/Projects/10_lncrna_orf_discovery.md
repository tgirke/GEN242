---
title: "lncRNAs and other features"
linkTitle: "lncRNAs and ORFs"
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
    + Parses DNA sequences of identified peak footprints
    + Identify in the identified peak sequences 1-2 of the following feature types: 
        + Long non-coding RNAs (lncRNAs; Han et al., 2019; Hu et al., 2017)
        + Open reading frames (ORFs)
        + miRNAs
        + Repeats

## References

+ Kaufmann, K, F Wellmer, J M Muiño, T Ferrier, S E Wuest, V Kumar, A Serrano-Mislata, et al. 2010. “Orchestration of Floral Initiation by APETALA1.” Science 328 (5974): 85–89. [PubMed](https://pubmed.ncbi.nlm.nih.gov/20360106/)
+ Han S, Liang Y, Ma Q, Xu Y, Zhang Y, Du W, Wang C, Li Y (2019) LncFinder: an integrated platform for long non-coding RNA identification utilizing sequence intrinsic composition, structural information and physicochemical property. Brief Bioinform 20: 2009–2027 [PubMed](https://pubmed.ncbi.nlm.nih.gov/30084867/)
+ Hu L, Xu Z, Hu B, Lu ZJ (2017) COME: a robust coding potential calculation tool for lncRNA identification and characterization based on multiple features. Nucleic Acids Res 45: e2 [PubMed](https://pubmed.ncbi.nlm.nih.gov/27608726/)



