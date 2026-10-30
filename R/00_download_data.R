#!/usr/bin/env Rscript
###############################################################################
# 00_download_data.R
# Downloads all OPEN-ACCESS datasets and records release/date in
# metadata/data_sources.csv.
#
#   TCGA PanCancer Atlas (UCEC, COAD/READ)  : cBioPortal DataHub tarballs
#   CPTAC (UCEC 2020, COAD 2019)             : cBioPortal DataHub tarballs
#   gnomAD constraint + variants for target gene    : gnomAD GraphQL API
#   ClinVar records for target gene                 : NCBI E-utilities
#   DepMap                                   : MANUAL download (see below)
#
# Controlled-access TCGA germline data (dbGaP phs000178) is NOT downloaded here.
# If approved, place it in data/controlled/ (never commit or share that folder).
###############################################################################

source(file.path("R", "00_config.R")); source(file.path("R", "utils.R"))
load_packages()

download_datahub <- function(study_id, dest_dir) {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  study_dir <- file.path(dest_dir, study_id)
  if (dir.exists(study_dir) && length(list.files(study_dir))) {
    log_msg("Already present: ", study_id); return(invisible(study_dir))
  }
  url <- sprintf("%s/%s.tar.gz", DATAHUB_URL, study_id)
  tgz <- file.path(dest_dir, paste0(study_id, ".tar.gz"))
  log_msg("Downloading ", url)
  utils::download.file(url, tgz, mode = "wb", quiet = FALSE)
  utils::untar(tgz, exdir = dest_dir)
  unlink(tgz)
  record_source(study_id, "cBioPortal DataHub (download date)", url, "open")
  invisible(study_dir)
}

# ---- 1. TCGA and CPTAC from cBioPortal DataHub -----------------------------
for (s in TCGA_STUDIES)  download_datahub(s, file.path(DIR_RAW, "tcga"))
for (s in CPTAC_STUDIES) download_datahub(s, file.path(DIR_RAW, "cptac"))

# ---- 2. gnomAD constraint and variants for target gene ----------------------------
gnomad_query <- function(query) {
  r <- httr::POST("https://gnomad.broadinstitute.org/api",
                  body = list(query = query), encode = "json",
                  httr::timeout(120))
  httr::stop_for_status(r)
  jsonlite::fromJSON(httr::content(r, as = "text", encoding = "UTF-8"),
                     simplifyVector = TRUE)
}

gn_dir <- file.path(DIR_RAW, "gnomad_clinvar")
q_gene <- sprintf('{
  gene(gene_symbol: "%s", reference_genome: GRCh38) {
    gene_id symbol chrom start stop
    gnomad_constraint { exp_lof obs_lof oe_lof oe_lof_lower oe_lof_upper pli
                        exp_mis obs_mis oe_mis mis_z exp_syn obs_syn oe_syn }
    variants(dataset: gnomad_r4) {
      variant_id pos ref alt consequence hgvsc hgvsp lof lof_filter flags
      exome  { ac an homozygote_count }
      genome { ac an homozygote_count }
    }
  }
}', GENE_TARGET)
res_gn <- tryCatch(gnomad_query(q_gene), error = function(e) {
  log_msg("gnomAD query failed: ", conditionMessage(e)); NULL })
if (!is.null(res_gn)) {
  saveRDS(res_gn, file.path(gn_dir, paste0("gnomad_", GENE_TARGET, "_r4.rds")))
  jsonlite::write_json(res_gn, file.path(gn_dir, paste0("gnomad_", GENE_TARGET, "_r4.json")),
                       auto_unbox = TRUE, pretty = TRUE)
  record_source(paste0("gnomAD_", GENE_TARGET), "gnomad_r4 (GraphQL API)",
                "https://gnomad.broadinstitute.org/api", "open")
}

# ---- 3. ClinVar records for target gene via E-utilities ---------------------------
eutils <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
cv_ids <- tryCatch({
  r <- httr::GET(sprintf("%s/esearch.fcgi", eutils),
                 query = list(db = "clinvar", term = paste0(GENE_TARGET, "[gene]"),
                              retmax = 5000, retmode = "json"))
  httr::stop_for_status(r)
  jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"))$esearchresult$idlist
}, error = function(e) { log_msg("ClinVar search failed: ", conditionMessage(e)); NULL })

if (length(cv_ids)) {
  chunks <- split(cv_ids, ceiling(seq_along(cv_ids) / 200))
  recs <- list()
  for (ch in chunks) {
    r <- httr::GET(sprintf("%s/esummary.fcgi", eutils),
                   query = list(db = "clinvar", id = paste(ch, collapse = ","),
                                retmode = "json"))
    httr::stop_for_status(r)
    js <- jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"),
                             simplifyVector = FALSE)$result
    for (id in setdiff(names(js), "uids")) {
      x <- js[[id]]
      germ <- x$germline_classification
      genes <- vapply(x$genes, function(g) g$symbol, character(1))
      recs[[id]] <- data.frame(
        clinvar_id = id,
        title = x$title %||% NA,
        variant_type = x$obj_type %||% NA,
        molecular_consequence = paste(unlist(x$molecular_consequence_list), collapse = ";"),
        classification = germ$description %||% NA,
        review_status = germ$review_status %||% NA,
        conditions = paste(vapply(germ$trait_set %||% list(), function(t) t$trait_name %||% "",
                                  character(1)), collapse = ";"),
        genes = paste(genes, collapse = ";")
      )
    }
    Sys.sleep(0.4)  # NCBI rate limit without API key
  }
  clinvar <- do.call(rbind, recs)
  write.csv(clinvar, file.path(gn_dir, paste0("clinvar_", GENE_TARGET, "_esummary.csv")), row.names = FALSE)
  record_source(paste0("ClinVar_", GENE_TARGET), "NCBI E-utilities (download date)", eutils, "open")
}

# ---- 4. DepMap (manual) ----------------------------------------------------
dm_dir <- file.path(DIR_RAW, "depmap")
have <- basename(list.files(dm_dir))
has_model <- "Model.csv" %in% have
has_expr <- any(grepl("^OmicsExpression(TPMLogp1HumanProteinCodingGenes|ProteinCodingGenesTPMLogp1)", have))
needed <- c("Model.csv", "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv (26Q3+)")
if (!(has_model && has_expr)) {
  message(
    "\nDepMap files must be downloaded manually from https://depmap.org/portal/data_page/",
    "\n(choose the latest public release, e.g. 'DepMap Public 25Q2') and placed in:\n  ", dm_dir,
    "\nRequired: ", paste(needed, collapse = ", "),
    "\nOptional: OmicsGlobalSignatures.csv (MSI score), CRISPRGeneEffect.csv,",
    "\n          CCLE proteomics (Nusinow 2020) for protein-level correlation.",
    "\nThen add the release name to metadata/data_sources.csv.\n")
} else {
  log_msg("DepMap required files present.")
}

save_session("00_download_data")
log_msg("Download step finished.")
