###############################################################################
# 00_config.R
# Shared configuration. Every analysis script sources this file first.
#
# To analyse a different gene, cancer type or partner genes, edit the
# "USER SETTINGS" block below, or set environment variables before running:
#   TARGET_GENE=ATM PARTNER_GENES=MLH1,MSH2 CNA_PROXY_GENE=ATM Rscript R/run_all.R
###############################################################################

env_vec <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (nzchar(v)) trimws(strsplit(v, ",")[[1]]) else default
}

# ============================ USER SETTINGS =================================

# Gene whose association with mismatch repair (MMR) deficiency is tested
GENE_TARGET   <- env_vec("TARGET_GENE", "DRG2")

# Genes expected to change with the target (correlation tests, plots)
GENES_PARTNER <- env_vec("PARTNER_GENES", c("MLH1", "MLH3"))

# Gene on the same chromosome arm as the target, used as an arm-level
# copy-number proxy (DRG2 and TP53 are both on 17p). Use NA to skip.
CNA_PROXY_GENE <- env_vec("CNA_PROXY_GENE", "TP53")
if (identical(CNA_PROXY_GENE, "NA")) CNA_PROXY_GENE <- NA_character_

# cBioPortal DataHub study IDs (https://www.cbioportal.org/datasets).
# Names are short labels used in tables and plots.
TCGA_STUDIES <- c(
  UCEC = "ucec_tcga_pan_can_atlas_2018",
  CRC  = "coadread_tcga_pan_can_atlas_2018"
)
CPTAC_STUDIES <- c(
  UCEC = "ucec_cptac_2020",
  COAD = "coad_cptac_2019"
)
# DepMap lineages treated as subgroups (regular expressions on OncotreeLineage)
DEPMAP_SUBGROUPS <- c(
  Colorectal  = "Bowel|Colorectal|Colon",
  Endometrial = "Uterus|Endometri"
)

# ========================== END USER SETTINGS ===============================

# ---- Project root ----------------------------------------------------------
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = getwd()), mustWork = TRUE)
if (basename(PROJECT_DIR) == "R") PROJECT_DIR <- dirname(PROJECT_DIR)

DIR_RAW       <- file.path(PROJECT_DIR, "data", "raw")
DIR_PROCESSED <- file.path(PROJECT_DIR, "data", "processed")
DIR_CONTROL   <- file.path(PROJECT_DIR, "data", "controlled")
# Results are kept per target gene so several genes can be analysed side by side
DIR_RESULTS   <- file.path(PROJECT_DIR, "results", GENE_TARGET)
DIR_TABLES    <- file.path(DIR_RESULTS, "tables")
DIR_FIGURES   <- file.path(DIR_RESULTS, "figures")
DIR_QC        <- file.path(DIR_RESULTS, "qc")
DIR_SESSION   <- file.path(DIR_RESULTS, "sessionInfo")
DIR_META      <- file.path(PROJECT_DIR, "metadata")
for (d in c(DIR_RAW, DIR_PROCESSED, DIR_CONTROL, DIR_TABLES, DIR_FIGURES,
            DIR_QC, DIR_SESSION, DIR_META)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- Fixed gene sets -------------------------------------------------------
GENES_MMR_CORE <- c("MLH1", "MSH2", "MSH6", "PMS2", "EPCAM")   # Lynch syndrome genes
GENES_MMR_EXT  <- c("MLH3", "MSH3", "PMS1", "EXO1", "POLE", "POLD1", "MUTYH")
GENES_FOCUS    <- unique(c(GENE_TARGET, GENES_PARTNER, "MLH1", "MLH3", "MSH2",
                           "MSH6", "PMS2", "MSH3", "PMS1", "EXO1"))

DATAHUB_URL <- "https://datahub.assets.cbioportal.org"

# ---- Cohort definition thresholds (pre-specified; see metadata/cohort_definitions.md)
MSI_MANTIS_CUTOFF <- 0.4   # MANTIS score > 0.4 = MSI-H
MSI_SENSOR_CUTOFF <- 10    # MSIsensor score >= 10 = MSI-H
MLH1_METH_CUTOFF  <- 0.3   # mean beta of MLH1 TSS200 probes > 0.3 = methylated
DEPMAP_MSI_CUTOFF <- 20    # DepMap MSIsensor2 score >= 20 = MSI (DepMap definition)
TRUNCATING_CLASSES <- c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
                        "Splice_Site", "Splice_Region", "Translation_Start_Site",
                        "Nonstop_Mutation")

# ---- Statistics ------------------------------------------------------------
FDR_CUTOFF <- 0.05
COR_METHOD <- "spearman"
SEED <- 20261029
set.seed(SEED)

# ---- Packages --------------------------------------------------------------
REQUIRED_PKGS <- c("data.table", "ggplot2", "jsonlite", "httr")
load_packages <- function(pkgs = REQUIRED_PKGS) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing packages: ", paste(missing, collapse = ", "),
         "\nInstall with install.packages(c(", paste0('"', missing, '"', collapse = ", "), "))",
         "\nor restore the pinned environment with renv::restore().")
  }
  suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))
}

options(stringsAsFactors = FALSE, scipen = 999, timeout = 3600)
message(sprintf("[config] target=%s | partners=%s | arm proxy=%s | results -> %s",
                GENE_TARGET, paste(GENES_PARTNER, collapse = ","), CNA_PROXY_GENE, DIR_RESULTS))
