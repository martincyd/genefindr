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

In modern genomics research, bioinformatic analyses such as RNA-seq often identify genes of interest that require downstream experimental validation. However, characterization of these candidate genes is a time-consuming process that commonly requires querying several disparate databases independently. genefindr is an R package that addresses this challenge by enabling rapid, comprehensive gene characterization through a single function call. Aggregating information from eight publicly available databases, MyGene.info [@mygene2016], Open Targets [@opentargets2023], Human Protein Atlas [@hpa2017], UniProt [@uniprot2023], GTEx [@gtex2020], cBioPortal [@cbioportal2012], PubMed [@pubmed], and ClinVar [@clinvar2014], genefindr provides a structured overview of target genes including functional summaries, disease association scores, somatic mutation frequencies, protein evidence, antibody availability, and recommended positive-control tissues for experimental validation.

# Statement of Need

Characterizing candidate genes identified in bioinformatic pipelines is a lengthy process that requires searching multiple databases independently to understand gene function within specific disease contexts. To assess mutation frequency, researchers commonly query TCGA via cBioPortal or COSMIC [@forbes2017]. To evaluate clinical relevance, ClinVar and literature searches are typically required. If a gene is considered a candidate for functional validation using techniques such as western blotting, immunofluorescence, or RT-PCR, researchers must also identify commercially available antibodies or appropriate positive-control tissues. In the absence of suitable antibodies for these assays, additional evaluation of alternative genes may be required. As such, this workflow is both time-consuming and fragmented, particularly during preliminary stages of bioinformatic analyses, as it often requires several web portals and interfaces.

genefindr streamlines this initial gene characterization process and is designed for researchers who use bioinformatic analyses to inform downstream wet-lab experimentation. It is not intended to replace comprehensive literature searches or formal differential expression analyses such as DESeq2 [@love2014] or edgeR [@robinson2010], but rather to complement them by providing a rapid preliminary characterization of candidate genes. genefindr enables simultaneous assessment of multiple genes within diverse disease contexts and generates structured output tables for organized result management while remaining in the R environment.

# Functionality

genefindr's core function, `findr()`, accepts a gene symbol and optional disease context to return a structured summary generated from eight databases. Users may input a vector of disease contexts or cancer tissue sites for cross-comparisons of a single gene, enabling identification of context-specific associations. The `findr_multi()` function extends this by accepting vectors of gene symbols for batch characterization. Results are returned as structured data frames suitable for downstream analysis, visualization, or export.

Disease associations are derived from Open Targets [@opentargets2023] and mutation frequencies from TCGA PanCancer Atlas 2018 via cBioPortal [@cbioportal2012]. Normal tissue expression data from GTEx [@gtex2020] are included to highlight tissues with relatively high expression of target genes, suggesting suitable positive-control tissues for downstream experimental design. Variant counts from ClinVar [@clinvar2014] are used to contextualize the clinical significance and mutational landscape of each gene. Supported cancer sites include breast, prostate, lung, colon, ovarian, liver, brain, pancreatic, skin, and blood. Though mutation frequency data from cBioPortal and predefined site terms are currently optimized for oncological contexts, disease associations from Open Targets and the `disease` parameter support broader queries including neurological, autoimmune, and cardiovascular conditions, making genefindr applicable beyond cancer biology. 

Integrating information from these complementary resources into a single interface promotes reproducible gene characterization and prioritization of candidate genes for experimental validation while reducing reliance on manual reference tracking and repetitive data queries.



```r
library(genefindr)

# Characterize a single gene in disease context
findr("TP53", site = "breast")

findr("APOE", disease = "alzheimer")

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
