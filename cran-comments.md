## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes for CRAN reviewers

- This is a first submission
- All examples are wrapped in \dontrun{} as they require internet access to query external APIs
- The package queries eight public, free, open-access databases: MyGene.info, Open Targets, Human Protein Atlas, UniProt, GTEx, cBioPortal, PubMed, and ClinVar
- All API calls include timeout handling (15 seconds) to prevent hanging
- The Sys.sleep() calls between NCBI queries are intentional to respect NCBI's rate limiting guidelines (max 3 requests/second without API key)
