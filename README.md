# genefindr <img src="man/figures/logo.png" align="right" height="139" />

> Rapid gene characterization without leaving R.

**genefindr** is an R package that provides instant, comprehensive gene characterization by querying eight public databases simultaneously. Designed for researchers who need quick, reliable gene overviews to support experimental decisions — from western blot candidate selection to target prioritization.

## Why genefindr?

Characterizing a gene typically means manually searching GeneCards, then Open Targets, then Human Protein Atlas, then PubMed... **genefindr does all of this in one function call.**

## Installation

```r
# Install from GitHub
pak::pak("martincyd/genefindr")

# Or install from local folder
pak::pak("local::path/to/genefindr")
```

## Quick Start

```r
library(genefindr)

# Single gene, single context
findr("TP53", site = "breast")

# Multiple genes
findr_multi(c("TP53", "BRCA1", "MYC"), site = "breast")

# Multi-site comparison
findr("TP53", site = c("breast", "lung", "colon"))

# Export results as a table
results <- findr_multi(c("TP53", "BRCA1"), site = "breast", output = "table")
write.csv(results, "candidates.csv")

# Non-cancer disease context
findr("APOE", disease = "alzheimer")

# Works with non-coding RNAs too
findr("MALAT1", site = "lung")
```

## Example Output

```
============================
Gene: TP53
Context: breast
============================
Name: tumor protein p53
Type: protein-coding
Molecular weight: 43,653 Da
Subcellular location: Cytoplasm, Nucleus, PML body, Endoplasmic reticulum...
Isoforms: 9 (p53, p53gamma, Del40-p53, Del40-p53beta...)
Note: multiple isoforms detected — verify antibody targets correct isoform

PubMed publications: 38,434 total | 4,636 in context of breast

Summary: This gene encodes a tumor suppressor protein...

--- Disease Association (Open Targets) ---
Association score: 0.743
Matched disease: Hereditary breast cancer

--- Protein & Expression (Human Protein Atlas) ---
Protein evidence: Evidence at protein level
Antibody available: Yes

--- Normal Expression (GTEx) ---
Suggested positive control tissue: Skin Sun Exposed Lower leg (37 TPM)

--- Tumor Mutation Frequency (cBioPortal/TCGA) ---
Mutated in 32.6% of TCGA BREAST tumors (353/1084 samples, TCGA PanCancer Atlas 2018)

--- Clinical Variants (ClinVar) ---
Pathogenic: 1,734 | Benign: 1,962 | VUS: 2,252

============================
Data sources: MyGene.info, Open Targets, Human Protein Atlas, UniProt, GTEx, cBioPortal, PubMed, ClinVar
============================
```

## Data Sources

| Database | Data provided |
|----------|--------------|
| [MyGene.info](https://mygene.info) | Gene name, type, summary |
| [Open Targets](https://www.opentargets.org) | Disease association scores |
| [Human Protein Atlas](https://www.proteinatlas.org) | Protein evidence, antibody availability |
| [UniProt](https://www.uniprot.org) | Molecular weight, subcellular location, isoforms |
| [GTEx](https://gtexportal.org/home/index.html) | Normal tissue expression |
| [cBioPortal/TCGA](https://www.cbioportal.org) | Tumor mutation frequency |
| [PubMed](https://pubmed.ncbi.nlm.nih.gov) | Publication counts |
| [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar) | Clinical variant counts |

## Supported Sites

`breast`, `prostate`, `lung`, `colon`, `ovarian`, `liver`, `brain`, `pancreatic`, `skin`, `blood`

## Note on gene symbols

genefindr requires official HGNC gene symbols. Common aliases may return unexpected results:

| Common name | Official symbol |
|-------------|----------------|
| HER2 | ERBB2 |
| p53 | TP53 |
| cMYC | MYC |
| VEGF | VEGFA |
| HER1 | EGFR |

When in doubt, look up the official symbol at [genenames.org](https://www.genenames.org).

## Functions

- `findr()` — characterize a single gene
- `findr_multi()` — characterize multiple genes at once

## Citation

If you use genefindr in your research, please cite the data sources above and this package:
Martin, Cydnie (2026). genefindr: Rapid Gene Characterization Using Public Genomic Databases. R package version 0.0.0.9000. https://github.com/martincyd/genefindr

## License

GPL-3


