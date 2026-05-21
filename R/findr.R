#' Characterize a gene using public genomic databases
#'
#' @param gene Character. Gene symbol (e.g. "TP53")
#' @param disease Character. Disease context (e.g. "alzheimer")
#' @param site Character. Cancer site (e.g. "breast", "lung"). Can be a vector for multi-site comparison.
#' @return Invisibly returns a data frame of results
#' @export
#' @examples
#' \dontrun{
#' findr("TP53", site = "breast")
#' }

findr <- function(gene, disease = NULL, site = NULL) {

  site_terms <- list(
    breast     = c("breast cancer", "breast carcinoma", "breast adenocarcinoma", "hereditary breast"),
    prostate   = c("prostate cancer", "prostate carcinoma", "prostate adenocarcinoma"),
    lung       = c("lung cancer", "lung carcinoma", "lung adenocarcinoma", "lung squamous"),
    colon      = c("colorectal", "colon cancer", "colorectal carcinoma", "colorectal adenocarcinoma"),
    ovarian    = c("ovarian", "ovarian cancer", "ovarian carcinoma", "ovarian neoplasm"),
    liver      = c("hepatocellular", "liver cancer", "hepatocellular carcinoma"),
    brain      = c("glioma", "glioblastoma", "brain cancer", "alzheimer", "neuroblastoma"),
    pancreatic = c("pancreatic", "pancreatic carcinoma", "pancreatic neoplasm"),
    skin       = c("melanoma", "skin cancer", "basal cell carcinoma"),
    blood      = c("leukemia", "lymphoma", "myeloma", "acute myeloid", "chronic lymphocytic")
  )

  tcga_studies <- list(
    breast     = "brca_tcga_pan_can_atlas_2018",
    prostate   = "prad_tcga_pan_can_atlas_2018",
    lung       = "luad_tcga_pan_can_atlas_2018",
    colon      = "coadread_tcga_pan_can_atlas_2018",
    ovarian    = "ov_tcga_pan_can_atlas_2018",
    liver      = "lihc_tcga_pan_can_atlas_2018",
    brain      = "gbm_tcga_pan_can_atlas_2018",
    pancreatic = "paad_tcga_pan_can_atlas_2018",
    skin       = "skcm_tcga_pan_can_atlas_2018",
    blood      = "laml_tcga_pan_can_atlas_2018"
  )

  # handle multiple sites
  if (!is.null(site) && length(site) > 1) {
    results <- list()
    first <- findr(gene, site = site[1])
    results[[site[1]]] <- first
    for (s in site[-1]) {
      cat("\n--- Context:", s, "---\n")
      res <- findr_context_only(gene, site = s)
      results[[s]] <- res
      cat("\n")
    }
    return(invisible(NULL))
  }

  # resolve site or disease to search terms
  if (!is.null(site)) {
    site_clean <- tolower(trimws(site))
    search_terms <- site_terms[[site_clean]]
    if (is.null(search_terms)) {
      cat("Unknown site '", site, "'. Try: breast, prostate, lung, colon, ovarian, liver, brain, pancreatic, skin, blood\n")
      return(invisible(NULL))
    }
    context_label <- site
    tcga_study <- tcga_studies[[site_clean]]
  } else if (!is.null(disease)) {
    search_terms <- c(disease)
    context_label <- disease
    tcga_study <- NULL
    site_clean <- NULL
  } else {
    search_terms <- NULL
    context_label <- NULL
    tcga_study <- NULL
    site_clean <- NULL
  }

  # query MyGene.info
  url <- paste0("https://mygene.info/v3/query?q=", gene, "&fields=name,summary,type_of_gene,ensembl.gene,entrezgene,symbol&species=human")

  result <- tryCatch({
    response <- httr2::request(url) |> httr2::req_perform()
    httr2::resp_body_json(response)
  }, error = function(e) {
    cat("Error: Could not reach MyGene.info. Please check your internet connection.\n")
    return(NULL)
  })

  if (is.null(result)) return(invisible(NULL))

  if (length(result$hits) == 0) {
    cat("Error: Gene '", gene, "' not found. Please check the gene symbol and try again.\n", sep="")
    return(invisible(NULL))
  }

  hit <- result$hits[[1]]
  ensembl_id <- hit$ensembl$gene
  entrez_id <- as.integer(hit$entrezgene)
  official_symbol <- hit$symbol

  if (hit$type_of_gene == "biological-region") {
    cat("Warning: '", gene, "' returned a biological region, not a gene. Try the official gene symbol instead.\n", sep="")
    return(invisible(NULL))
  }

  # query Human Protein Atlas
  hpa_url <- paste0("https://www.proteinatlas.org/api/search_download.php?search=", gene, "&format=json&columns=g,eg,up,pe,ab,rnacancertype&compress=no")
  hpa_response <- httr2::request(hpa_url) |> httr2::req_perform()
  hpa_result <- httr2::resp_body_json(hpa_response)
  hpa <- if (length(hpa_result) > 0) hpa_result[[1]] else NULL

  # query UniProt
  mol_weight <- NULL
  subcell_loc <- NULL
  isoform_count <- NULL
  isoform_names <- NULL
  if (!is.null(hpa) && !is.null(hpa$Uniprot[[1]])) {
    uniprot_id <- hpa$Uniprot[[1]]
    uniprot_url <- paste0("https://rest.uniprot.org/uniprotkb/", uniprot_id, "?format=json")
    uniprot_response <- httr2::request(uniprot_url) |> httr2::req_perform()
    uniprot_result <- httr2::resp_body_json(uniprot_response)
    mol_weight <- uniprot_result$sequence$molWeight
    loc <- Filter(function(x) x$commentType == "SUBCELLULAR LOCATION", uniprot_result$comments)
    if (length(loc) > 0) {
      all_locs <- unlist(lapply(loc, function(l)
        sapply(l$subcellularLocations, function(x) x$location$value)
      ))
      base_locs <- unique(trimws(unlist(strsplit(all_locs, ","))))
      subcell_loc <- paste(unique(base_locs), collapse = ", ")
    }
    alt <- Filter(function(x) x$commentType == "ALTERNATIVE PRODUCTS", uniprot_result$comments)
    if (length(alt) > 0) {
      isoform_count <- length(alt[[1]]$isoforms)
      raw_names <- sapply(alt[[1]]$isoforms, function(x) {
        if (length(x$synonyms) > 0) x$synonyms[[1]]$value else x$name$value
      })
      isoform_names <- raw_names[!grepl("^[0-9]+$", raw_names) & !grepl("^[A-Z0-9]{5,}$", raw_names)]
    }
  }

  # query GTEx
  gtex_top_tissue <- NULL
  gtex_top_tpm <- NULL
  tryCatch({
    gene_info <- gtexr::get_genes(geneId = gene, .verbose = FALSE)
    if (nrow(gene_info) > 0) {
      gencode_id <- gene_info$gencodeId[1]
      gtex_result <- gtexr::get_median_gene_expression(gencodeId = gencode_id, datasetId = "gtex_v8", .verbose = FALSE)
      if (nrow(gtex_result) > 0) {
        real_tissues <- gtex_result[!grepl("Cells_", gtex_result$tissueSiteDetailId), ]
        top <- real_tissues[which.max(real_tissues$median), ]
        gtex_top_tissue <- gsub("_", " ", top$tissueSiteDetailId)
        gtex_top_tpm <- round(top$median, 1)
      }
    }
  }, error = function(e) NULL)

  # query cBioPortal
  mut_frequency <- NULL
  mut_count <- NULL
  total_samples <- NULL
  if (!is.null(tcga_study) && !is.null(entrez_id)) {
    tryCatch({
      sample_url <- paste0("https://www.cbioportal.org/api/sample-lists/", tcga_study, "_all")
      sample_resp <- httr2::request(sample_url) |> httr2::req_perform()
      total_samples <- httr2::resp_body_json(sample_resp)$sampleCount
      mut_url <- paste0("https://www.cbioportal.org/api/molecular-profiles/", tcga_study, "_mutations/mutations?sampleListId=", tcga_study, "_all&entrezGeneId=", entrez_id, "&pageSize=10000")
      mut_resp <- httr2::request(mut_url) |> httr2::req_perform()
      mut_count <- length(httr2::resp_body_json(mut_resp))
      mut_frequency <- round(mut_count / total_samples * 100, 1)
    }, error = function(e) NULL)
  }


  # query PubMed
  pubmed_total <- NULL
  pubmed_context <- NULL

  pubmed_total <- tryCatch({
    total_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=", gene, "&retmode=json&retmax=0")
    total_resp <- httr2::request(total_url) |> httr2::req_perform()
    as.integer(httr2::resp_body_json(total_resp)$esearchresult$count)
  }, error = function(e) NULL)

  Sys.sleep(0.5)

  if (!is.null(context_label)) {
    pubmed_context <- tryCatch({
      context_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=", gene, "+", gsub(" ", "+", context_label), "&retmode=json&retmax=0")
      context_resp <- httr2::request(context_url) |> httr2::req_perform()
      as.integer(httr2::resp_body_json(context_resp)$esearchresult$count)
    }, error = function(e) NULL)
  }

  # query ClinVar
  clinvar_pathogenic <- tryCatch({
    path_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=", official_symbol, "[gene]+pathogenic[clinical_significance]&retmode=json&retmax=0")
    as.integer(httr2::resp_body_json(httr2::req_perform(httr2::request(path_url)))$esearchresult$count)
  }, error = function(e) NULL)

  Sys.sleep(0.5)

  clinvar_benign <- tryCatch({
    ben_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=", official_symbol, "[gene]+benign[clinical_significance]&retmode=json&retmax=0")
    as.integer(httr2::resp_body_json(httr2::req_perform(httr2::request(ben_url)))$esearchresult$count)
  }, error = function(e) NULL)

  Sys.sleep(0.5)

  clinvar_vus <- tryCatch({
    vus_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=", official_symbol, "[gene]+%22uncertain+significance%22&retmode=json&retmax=0")
    as.integer(httr2::resp_body_json(httr2::req_perform(httr2::request(vus_url)))$esearchresult$count)
  }, error = function(e) NULL)

  # print basic info
  cat("============================\n")
  cat("Gene:", toupper(gene), "\n")
  if (!is.null(context_label)) cat("Context:", context_label, "\n")
  cat("============================\n")
  cat("Name:", hit$name, "\n")
  cat("Type:", hit$type_of_gene, "\n")
  if (!is.null(hit$type_of_gene) && hit$type_of_gene != "protein-coding") {
    cat("Note: non-protein-coding gene — protein-based fields (molecular weight, antibody, isoforms) not applicable\n")
  }
  if (!is.null(mol_weight)) cat("Molecular weight:", format(mol_weight, big.mark = ","), "Da\n")
  if (!is.null(subcell_loc)) cat("Subcellular location:", subcell_loc, "\n")
  if (!is.null(isoform_count)) {
    if (length(isoform_names) > 0) {
      cat("Isoforms:", isoform_count, "(", paste(isoform_names, collapse = ", "), ")\n")
    } else {
      cat("Isoforms:", isoform_count, "\n")
    }
    if (isoform_count > 3) cat("Note: multiple isoforms detected — verify antibody targets correct isoform\n")
  }
  if (!is.null(pubmed_total)) {
    cat("\n")
    cat("PubMed publications:", format(pubmed_total, big.mark = ","), "total")
    if (!is.null(pubmed_context)) cat(" |", format(pubmed_context, big.mark = ","), "in context of", context_label)
    cat("\n")
  }
  cat("\n")
  cat("Summary:", hit$summary, "\n")

  # disease association
  assoc_score <- NULL
  assoc_disease <- NULL
  if (!is.null(search_terms)) {
    cat("\n--- Disease Association (Open Targets) ---\n")
    query <- paste0('{
      target(ensemblId: "', ensembl_id, '") {
        associatedDiseases {
          rows { score disease { name } }
        }
      }
    }')
    ot_response <- httr2::request("https://api.platform.opentargets.org/api/v4/graphql") |>
      httr2::req_body_json(list(query = query)) |>
      httr2::req_perform()
    ot_result <- httr2::resp_body_json(ot_response)
    rows <- ot_result$data$target$associatedDiseases$rows
    matches <- Filter(function(r) {
      any(sapply(search_terms, function(term) grepl(term, r$disease$name, ignore.case = TRUE)))
    }, rows)
    if (length(matches) > 0) {
      assoc_score <- round(matches[[1]]$score, 3)
      assoc_disease <- matches[[1]]$disease$name
      cat("Association score:", assoc_score, "\n")
      cat("Matched disease:", assoc_disease, "\n")
    } else {
      cat("No direct association found for '", context_label, "' in top results\n")
    }
  }

  # HPA output
  cat("\n--- Protein & Expression (Human Protein Atlas) ---\n")
  protein_evidence <- NULL
  antibody_available <- FALSE
  if (!is.null(hpa)) {
    protein_evidence <- hpa$Evidence
    antibody <- hpa$Antibody
    antibody_available <- length(antibody) > 0 && any(antibody != "")
    cat("Protein evidence:", protein_evidence, "\n")
    cat("Antibody available:", ifelse(antibody_available, "Yes", "No"), "\n")
  } else {
    cat("No Human Protein Atlas data found\n")
  }

  # GTEx output
  cat("\n--- Normal Expression (GTEx) ---\n")
  if (!is.null(gtex_top_tissue)) {
    cat("Suggested positive control tissue:", gtex_top_tissue, "(", gtex_top_tpm, "TPM)\n")
  } else {
    cat("GTEx data not available\n")
  }

  # cBioPortal output
  if (!is.null(mut_frequency)) {
    cat("\n--- Tumor Mutation Frequency (cBioPortal/TCGA) ---\n")
    if (mut_frequency > 100) {
      cat("Mutation count:", mut_count, "across", total_samples, "TCGA", toupper(site_clean), "tumors — exceeds 100% suggesting multiple mutations per sample (common in very large genes)\n")
    } else {
      cat("Mutated in", mut_frequency, "% of TCGA", toupper(site_clean), "tumors (", mut_count, "/", total_samples, "samples, TCGA PanCancer Atlas 2018)\n")
    }
  }

  # ClinVar output
  if (!is.null(clinvar_pathogenic)) {
    cat("\n--- Clinical Variants (ClinVar) ---\n")
    cat("Pathogenic:", format(clinvar_pathogenic, big.mark = ","),
        "| Benign:", format(clinvar_benign, big.mark = ","),
        "| VUS:", format(clinvar_vus, big.mark = ","), "\n")
  }

  cat("\n============================\n")
  cat("Data sources: MyGene.info, Open Targets, Human Protein Atlas, UniProt, GTEx, cBioPortal, PubMed, ClinVar\n")
  cat("============================\n")

  invisible(data.frame(
    gene = toupper(gene),
    context = ifelse(is.null(context_label), NA, context_label),
    molecular_weight_da = ifelse(is.null(mol_weight), NA, mol_weight),
    subcellular_location = ifelse(is.null(subcell_loc), NA, subcell_loc),
    isoforms = ifelse(is.null(isoform_count), NA, isoform_count),
    pubmed_total = ifelse(is.null(pubmed_total), NA, pubmed_total),
    pubmed_in_context = ifelse(is.null(pubmed_context), NA, pubmed_context),
    association_score = ifelse(is.null(assoc_score), NA, assoc_score),
    matched_disease = ifelse(is.null(assoc_disease), NA, assoc_disease),
    protein_evidence = ifelse(is.null(protein_evidence), NA, protein_evidence),
    antibody_available = antibody_available,
    mutation_frequency_pct = ifelse(is.null(mut_frequency), NA, mut_frequency),
    mut_count = ifelse(is.null(mut_count), NA, mut_count),
    total_samples = ifelse(is.null(total_samples), NA, total_samples),
    expected_band_kda = ifelse(is.null(mol_weight), NA, round(mol_weight / 1000)),
    positive_control_tissue = ifelse(is.null(gtex_top_tissue), NA, gtex_top_tissue),
    gtex_tpm = ifelse(is.null(gtex_top_tpm), NA, gtex_top_tpm),
    clinvar_pathogenic = ifelse(is.null(clinvar_pathogenic), NA, clinvar_pathogenic),
    clinvar_benign = ifelse(is.null(clinvar_benign), NA, clinvar_benign),
    clinvar_vus = ifelse(is.null(clinvar_vus), NA, clinvar_vus),
    stringsAsFactors = FALSE
  ))
}

#' Characterize multiple genes using public genomic databases
#'
#' @param genes Character vector. Gene symbols (e.g. c("TP53", "BRCA1"))
#' @param disease Character. Disease context
#' @param site Character. Cancer site
#' @param output Character. "print" or "table"
#' @return Invisibly returns a data frame of results
#' @export
#' @examples
#' \dontrun{
#' findr_multi(c("TP53", "BRCA1"), site = "breast")
#' }

findr_multi <- function(genes, disease = NULL, site = NULL, output = "print") {
  results <- list()
  for (gene in genes) {
    result <- findr(gene, disease = disease, site = site)
    results[[gene]] <- result
    cat("\n")
  }
  if (output == "table") {
    return(do.call(rbind, results))
  }
  invisible(do.call(rbind, results))
}

findr_context_only <- function(gene, site) {

  site_terms <- list(
    breast     = c("breast cancer", "breast carcinoma", "breast adenocarcinoma", "hereditary breast"),
    prostate   = c("prostate cancer", "prostate carcinoma", "prostate adenocarcinoma"),
    lung       = c("lung cancer", "lung carcinoma", "lung adenocarcinoma", "lung squamous"),
    colon      = c("colorectal", "colon cancer", "colorectal carcinoma", "colorectal adenocarcinoma"),
    ovarian    = c("ovarian", "ovarian cancer", "ovarian carcinoma", "ovarian neoplasm"),
    liver      = c("hepatocellular", "liver cancer", "hepatocellular carcinoma"),
    brain      = c("glioma", "glioblastoma", "brain cancer", "alzheimer", "neuroblastoma"),
    pancreatic = c("pancreatic", "pancreatic carcinoma", "pancreatic neoplasm"),
    skin       = c("melanoma", "skin cancer", "basal cell carcinoma"),
    blood      = c("leukemia", "lymphoma", "myeloma", "acute myeloid", "chronic lymphocytic")
  )

  tcga_studies <- list(
    breast     = "brca_tcga_pan_can_atlas_2018",
    prostate   = "prad_tcga_pan_can_atlas_2018",
    lung       = "luad_tcga_pan_can_atlas_2018",
    colon      = "coadread_tcga_pan_can_atlas_2018",
    ovarian    = "ov_tcga_pan_can_atlas_2018",
    liver      = "lihc_tcga_pan_can_atlas_2018",
    brain      = "gbm_tcga_pan_can_atlas_2018",
    pancreatic = "paad_tcga_pan_can_atlas_2018",
    skin       = "skcm_tcga_pan_can_atlas_2018",
    blood      = "laml_tcga_pan_can_atlas_2018"
  )

  site_clean <- tolower(trimws(site))
  search_terms <- site_terms[[site_clean]]
  tcga_study <- tcga_studies[[site_clean]]

  # get entrez and ensembl
  url <- paste0("https://mygene.info/v3/query?q=", gene, "&fields=ensembl.gene,entrezgene&species=human")
  response <- httr2::request(url) |> httr2::req_perform()
  result <- httr2::resp_body_json(response)
  hit <- result$hits[[1]]
  ensembl_id <- hit$ensembl$gene
  entrez_id <- as.integer(hit$entrezgene)

  # Open Targets
  assoc_score <- NULL
  assoc_disease <- NULL
  cat("\n--- Disease Association (Open Targets) ---\n")
  query <- paste0('{
      target(ensemblId: "', ensembl_id, '") {
        associatedDiseases {
          rows { score disease { name } }
        }
      }
    }')
  ot_response <- httr2::request("https://api.platform.opentargets.org/api/v4/graphql") |>
    httr2::req_body_json(list(query = query)) |>
    httr2::req_perform()
  ot_result <- httr2::resp_body_json(ot_response)
  rows <- ot_result$data$target$associatedDiseases$rows
  matches <- Filter(function(r) {
    any(sapply(search_terms, function(term) grepl(term, r$disease$name, ignore.case = TRUE)))
  }, rows)
  if (length(matches) > 0) {
    assoc_score <- round(matches[[1]]$score, 3)
    assoc_disease <- matches[[1]]$disease$name
    cat("Association score:", assoc_score, "\n")
    cat("Matched disease:", assoc_disease, "\n")
  } else {
    cat("No direct association found for '", site, "' in top results\n")
  }

  # PubMed context count
  pubmed_context <- NULL
  tryCatch({
    context_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=", gene, "+", site, "&retmode=json&retmax=0")
    context_resp <- httr2::request(context_url) |> httr2::req_perform()
    pubmed_context <- as.integer(httr2::resp_body_json(context_resp)$esearchresult$count)
    cat("PubMed publications in context of", site, ":", format(pubmed_context, big.mark = ","), "\n")
  }, error = function(e) NULL)

  # cBioPortal
  mut_frequency <- NULL
  mut_count <- NULL
  total_samples <- NULL
  tryCatch({
    sample_url <- paste0("https://www.cbioportal.org/api/sample-lists/", tcga_study, "_all")
    sample_resp <- httr2::request(sample_url) |> httr2::req_perform()
    total_samples <- httr2::resp_body_json(sample_resp)$sampleCount
    mut_url <- paste0("https://www.cbioportal.org/api/molecular-profiles/", tcga_study, "_mutations/mutations?sampleListId=", tcga_study, "_all&entrezGeneId=", entrez_id, "&pageSize=10000")
    mut_resp <- httr2::request(mut_url) |> httr2::req_perform()
    mut_count <- length(httr2::resp_body_json(mut_resp))
    mut_frequency <- round(mut_count / total_samples * 100, 1)
  }, error = function(e) NULL)

  if (!is.null(mut_frequency)) {
    cat("\n--- Tumor Mutation Frequency (cBioPortal/TCGA) ---\n")
    cat("Mutated in", mut_frequency, "% of TCGA", toupper(site_clean), "tumors (", mut_count, "/", total_samples, "samples, TCGA PanCancer Atlas 2018)\n")
  }

  invisible(data.frame(
    gene = toupper(gene),
    context = site,
    association_score = ifelse(is.null(assoc_score), NA, assoc_score),
    matched_disease = ifelse(is.null(assoc_disease), NA, assoc_disease),
    pubmed_in_context = ifelse(is.null(pubmed_context), NA, pubmed_context),
    mutation_frequency_pct = ifelse(is.null(mut_frequency), NA, mut_frequency),
    mut_count = ifelse(is.null(mut_count), NA, mut_count),
    total_samples = ifelse(is.null(total_samples), NA, total_samples),
    stringsAsFactors = FALSE
  ))
}
