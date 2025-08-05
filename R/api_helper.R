# R/api_helpers.R
# for the comparison tab - mod_comparison.R


# UniProt
# Fetch UniProt JSON metadata
get_uniprot_protein_data <- function(pid) {
  url <- paste0("https://rest.uniprot.org/uniprotkb/", pid, ".json")
  res <- httr::GET(url)
  if (httr::status_code(res) == 200) {
    httr::content(res, as = "parsed", type = "application/json")
  } else {
    NULL
  }
}

# Extract sequence string from UniProt JSON
get_protein_sequence <- function(protein_data) {
  if (!is.null(protein_data$sequence$value)) {
    protein_data$sequence$value
  } else {
    NULL
  }
}

# Extract domain information from UniProt JSON
extract_uniprot_domains <- function(protein_data) {
  if (!is.null(protein_data$features)) {
    domains <- protein_data$features[sapply(protein_data$features, function(x) x$type == "Domain")]
    if (length(domains) > 0) {
      return(data.frame(
        start = as.integer(sapply(domains, function(x) x$location$start$value)),
        end   = as.integer(sapply(domains, function(x) x$location$end$value)),
        description = sapply(domains, function(x) x$description),
        stringsAsFactors = FALSE
      ))
    }
  }
  return(NULL)
}

# Ensembl
# Fetch transcript table for a UniProt ID
get_transcripts_by_uniprot <- function(uniprot_id) {
  ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  result <- biomaRt::getBM(
    attributes = c("uniprotswissprot", "ensembl_transcript_id", "external_gene_name", "transcript_biotype", "transcript_is_canonical"),
    filters = "uniprotswissprot",
    values = uniprot_id,
    mart = ensembl
  )
  subset(result, transcript_biotype == "protein_coding") |> unique()
}

# GnomAD
# Fetch gnomAD variants for a canonical transcript ID
get_gnomad_variants_for_transcript <- function(transcript_id) {
  url <- "https://gnomad.broadinstitute.org/api"
  query <- '
  query TranscriptVariants($transcriptID: String!) {
    transcript(transcript_id: $transcriptID, reference_genome: GRCh38) {
      variants(dataset: gnomad_r4) {
        variant_id
        pos
        consequence
        transcript_id
        hgvsp
        exome { af ac an }
      }
    }
  }'
  
  body <- list(
    query = query,
    variables = list(transcriptID = transcript_id)
  )
  
  res <- httr::POST(url, body = body, encode = "json", httr::content_type_json())
  if (httr::status_code(res) != 200) return(NULL)
  
  content <- httr::content(res, as = "parsed", simplifyVector = TRUE)
  vars <- tryCatch(content$data$transcript$variants, error = function(e) NULL)
  if (is.null(vars)) return(NULL)
  
  df <- data.frame(
    Variant_ID = vars$variant_id,
    Position = vars$pos,
    Consequence = vars$consequence,
    HGVSp = vars$hgvsp,
    Transcript_ID = vars$transcript_id,
    AF = vars$exome$af,
    AC = vars$exome$ac,
    AN = vars$exome$an,
    stringsAsFactors = FALSE
  )
  
  df$AA_Position <- sapply(df$HGVSp, function(hgvsp) {
    if (is.na(hgvsp) || hgvsp == "") return(NA)
    matches <- regmatches(hgvsp, regexec("\\d+", hgvsp))[[1]]
    if (length(matches) > 0) as.integer(matches[1]) else NA
  })
  
  return(df)
}

# PTM
# Extract PTM sites from UniProt JSON
get_ptm_data <- function(protein_data) {
  if (is.null(protein_data$features)) return(NULL)
  
  ptm_features <- protein_data$features[sapply(protein_data$features, function(x) x$type == "Modified residue")]
  if (length(ptm_features) == 0) return(NULL)
  
  df <- data.frame(
    Position = as.integer(sapply(ptm_features, function(x) x$location$start$value)),
    Type = sapply(ptm_features, function(x) x$description)
  )
  
  df$TypeCategory <- dplyr::case_when(
    grepl("phospho", df$Type, ignore.case = TRUE) ~ "Phosphorylation",
    grepl("acetyl", df$Type, ignore.case = TRUE) ~ "Acetylation",
    grepl("succinyl", df$Type, ignore.case = TRUE) ~ "Succinylation",
    grepl("methyl", df$Type, ignore.case = TRUE) ~ "Methylation",
    TRUE ~ "Other"
  )
  
  df <- df |>
    dplyr::group_by(Position, TypeCategory) |>
    dplyr::summarise(
      tooltip = paste(paste0(Position, "：", Type), collapse = "\n"), 
      .groups = "drop"
    )
  
  return(df)
}

# AlphaMissense
# Get AlphaMissense URL from AlphaFold API
fetch_alphamissense_csv_url <- function(uniprot_id) {
  url <- paste0("https://alphafold.ebi.ac.uk/api/prediction/", uniprot_id)
  res <- httr::GET(url)
  if (httr::status_code(res) == 200) {
    json <- httr::content(res, as = "parsed", type = "application/json")
    if (!is.null(json[[1]]$amAnnotationsUrl)) {
      return(json[[1]]$amAnnotationsUrl)
    }
  }
  return(NULL)
}

# Download and process AlphaMissense CSV
get_alphamissense_data <- function(uniprot_id) {
  url <- fetch_alphamissense_csv_url(uniprot_id)
  if (is.null(url)) return(NULL)
  
  temp_file <- tempfile(fileext = ".csv")
  tryCatch({
    download.file(url, destfile = temp_file, quiet = TRUE)
  }, error = function(e) {
    warning("AlphaMissense download failed.")
    return(NULL)
  })
  
  df <- readr::read_csv(temp_file, show_col_types = FALSE)
  
  df <- df %>%
    dplyr::mutate(
      aa_from = stringr::str_sub(protein_variant, 1, 1),
      position = as.integer(stringr::str_extract(protein_variant, "\\d+")),
      aa_to = stringr::str_sub(protein_variant, -1, -1),
      aa_to = factor(aa_to, levels = rev(c(
        "A", "C", "D", "E", "F", "G", "H", "I", "K", 
        "L", "M", "N", "P", "Q", "R", "S", "T", "V", 
        "W", "Y"
      ))),
      am_class = dplyr::recode(am_class,
                               "Amb" = "ambiguous",
                               "LBen" = "likely_benign",
                               "LPath" = "likely_pathogenic")
    ) %>%
    dplyr::rename(score = am_pathogenicity)
  
  return(df)
}