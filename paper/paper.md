---
title: 'genefindr: Rapid Gene Characterization Using Public Genomic Databases in R'
tags:
  - R
  - bioinformatics
  - genomics
  - gene characterization
  - cancer biology
authors:
  - name: Cydnie Martin
    orcid: 0009-0004-8479-6356
    affiliation: 1
affiliations:
  - name: Wayne State University School of Medicine, Detroit, USA
    index: 1
date: 20 June 2026
bibliography: paper.bib
---

# Summary

In modern genomics research, identifying and characterizing candidate genes is a critical but time-consuming process that commonly requires querying several disparate databases. genefindr is an R package that addresses this challenge by enabling rapid, comprehensive gene characterization through a single function call. Querying eight publicly available databases simultaneously -- MyGene.info [@mygene2019], Open Targets [@opentargets2023], Human Protein Atlas [@hpa2023], UniProt [@uniprot2023], GTEx [@gtex2020], cBioPortal [@cbioportal2012], PubMed [@pubmed], and ClinVar [@clinvar2014] -- genefindr provides a brief, structured overview of target genes including functional summaries, disease association scores, somatic mutation frequencies, protein evidence, antibody availability, and recommended positive control tissues for experimental validation.

# Statement of Need

Characterizing candidate genes during bioinformatic pipelines is a lengthy process that requires searching multiple databases to understand gene function within specific disease contexts. To assess mutation frequency, researchers commonly query TCGA via cBioPortal or COSMIC [@forbes2017]. To evaluate clinical relevance, ClinVar and literature searches are typically required. If a gene is considered a candidate for functional validation using techniques such as western blotting, immunofluorescence, or RT-PCR, researchers must additionally identify commercially available antibodies and appropriate positive control tissues. When no suitable antibody exists, the entire process must recommence with a different candidate gene. This workflow is time-consuming and fragmented, particularly during preliminary stages of bioinformatic pipelines.

genefindr streamlines this initial gene characterization process and is designed for basic scientists who use bioinformatic analyses to inform wet lab validation. It is not intended to replace formal differential expression analyses such as DESeq2 [@love2014] or edgeR [@robinson2010], but rather to complement them by providing a rapid preliminary characterization of candidate genes. genefindr enables simultaneous assessment of multiple genes within specific disease contexts, supports multi-site comparisons across ten cancer types, and generates structured output tables for organized result management, all without leaving the R environment.

# Functionality

genefindr's core function, `findr()`, accepts a gene symbol and optional disease context to return a structured summary rendered from eight databases. Users may input a vector of disease contexts or tissue sites for cross- comparisons of a single gene, enabling identification of context-specific associations. The `findr_multi()` function extends this by accepting vectors of gene symbols for batch characterization. Results can be returned as a formatted data frame for downstream use or export.

Supported cancer sites include breast, prostate, lung, colon, ovarian, liver, brain, pancreatic, skin, and blood. Disease associations are derived from Open Targets [@opentargets2023] and mutation frequencies from TCGA PanCancer Atlas 2018 via cBioPortal [@cbioportal2012]. Normal tissue expression data from GTEx [@gtex2020] provides suggested positive control tissues for experimental design. Variant counts from ClinVar [@clinvar2014] are used to contextualize the clinical significance and mutational landscape of each gene.

```r
library(genefindr)

# Characterize a single gene in disease context
findr("TP53", site = "breast")

# Batch characterization of multiple genes
results <- findr_multi(c("TP53", "BRCA1", "MYC"),
                       site = "breast",
                       output = "table")

# Multi-site comparison
findr("TP53", site = c("breast", "lung", "colon"))
```

# Acknowledgements

The author thanks Dr. Moray J. Campbell for support and guidance during the process of developing genefindr. The author also acknowledges the developers and maintainers of the open source databases that make genefindr possible.

# References
