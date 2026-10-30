#!/usr/bin/env Rscript
###############################################################################
# 05_drg2_population_variants.R
# Questions:
#   Q1. How tolerant is target gene to loss-of-function in the general population
#       (gnomAD pLI, LOEUF)? A dominant Lynch-like predisposition from
#       heterozygous target gene loss is less plausible if LoF carriers are common.
#   Q2. How many high-confidence LoF carriers are there (allele counts)?
#   Q3. What does ClinVar record for target gene (and any cancer-related conditions)?
# Note: large multi-gene CNVs in ClinVar (e.g. the 17p11.2 Smith-Magenis region
# for DRG2) are flagged separately because
# they span many genes.
###############################################################################

source(file.path("R", "00_config.R")); source(file.path("R", "utils.R"))
load_packages()

gdir <- file.path(DIR_RAW, "gnomad_clinvar")

# ---- Q1/Q2 gnomAD ----------------------------------------------------------
f_gn <- file.path(gdir, paste0("gnomad_", GENE_TARGET, "_r4.rds"))
if (file.exists(f_gn)) {
  gn <- readRDS(f_gn)$data$gene
  cons <- as.data.frame(gn$gnomad_constraint)
  cons$gene <- GENE_TARGET
  loeuf <- if ("oe_lof_upper" %in% names(cons)) cons$oe_lof_upper else NA
  cons$interpretation <- ifelse(is.na(loeuf), "LOEUF unavailable",
    ifelse(loeuf < 0.6, "LoF-intolerant (LOEUF < 0.6)",
    ifelse(loeuf < 1, "Moderately constrained (LOEUF 0.6 to 1)", "LoF-tolerant (LOEUF >= 1)")))
  write.csv(cons, file.path(DIR_TABLES, paste0("05_gnomAD_", GENE_TARGET, "_constraint.csv")), row.names = FALSE)

  v <- as.data.frame(gn$variants)
  if (nrow(v)) {
    ac <- function(x, f) { y <- x[[f]]; if (is.null(y)) rep(0, nrow(x)) else ifelse(is.na(y), 0, y) }
    v$exome_ac  <- ac(v$exome, "ac");  v$exome_an  <- ac(v$exome, "an")
    v$genome_ac <- ac(v$genome, "ac"); v$genome_an <- ac(v$genome, "an")
    v$exome_hom <- ac(v$exome, "homozygote_count"); v$genome_hom <- ac(v$genome, "homozygote_count")
    v$total_ac <- v$exome_ac + v$genome_ac
    v$total_an <- v$exome_an + v$genome_an
    v$af <- ifelse(v$total_an > 0, v$total_ac / v$total_an, NA)
    v$hom_total <- v$exome_hom + v$genome_hom
    v$flags <- vapply(v$flags, function(z) paste(unlist(z), collapse = ";"), character(1))
    v$hc_lof <- !is.na(v$lof) & v$lof == "HC" & (is.na(v$lof_filter) | v$lof_filter == "")
    v <- v[, setdiff(names(v), c("exome", "genome"))]
    write.csv(v, file.path(DIR_TABLES, paste0("05_gnomAD_", GENE_TARGET, "_variants.csv")), row.names = FALSE)

    lof <- v[v$hc_lof, ]
    summ <- data.frame(
      n_variants_total = nrow(v),
      n_hc_lof_variants = nrow(lof),
      hc_lof_allele_count = sum(lof$total_ac),
      hc_lof_homozygotes = sum(lof$hom_total),
      cumulative_hc_lof_af = sum(lof$af, na.rm = TRUE),
      approx_carrier_frequency = 2 * sum(lof$af, na.rm = TRUE),
      n_missense = sum(grepl("missense", v$consequence)))
    write.csv(summ, file.path(DIR_TABLES, paste0("05_gnomAD_", GENE_TARGET, "_LoF_summary.csv")), row.names = FALSE)
    print(cons); print(summ)
  }
} else log_msg("gnomAD file missing; run 00_download_data.R")

# ---- Q3 ClinVar ------------------------------------------------------------
f_cv <- file.path(gdir, paste0("clinvar_", GENE_TARGET, "_esummary.csv"))
if (file.exists(f_cv)) {
  cv <- read.csv(f_cv, stringsAsFactors = FALSE)
  cv$n_genes <- lengths(strsplit(cv$genes, ";"))
  cv$scope <- ifelse(cv$n_genes > 1, "Multi-gene (large CNV)", paste0(GENE_TARGET, "-specific"))
  cv$cancer_related <- grepl("cancer|carcinoma|neoplas|tumou?r|Lynch|polyposis",
                             cv$conditions, ignore.case = TRUE)
  write.csv(cv, file.path(DIR_TABLES, paste0("05_ClinVar_", GENE_TARGET, "_annotated.csv")), row.names = FALSE)
  tab <- as.data.frame(table(scope = cv$scope, classification = cv$classification))
  tab <- tab[tab$Freq > 0, ]
  write.csv(tab, file.path(DIR_TABLES, paste0("05_ClinVar_", GENE_TARGET, "_summary.csv")), row.names = FALSE)
  print(tab)
  log_msg(paste0(GENE_TARGET, "-specific ClinVar records mentioning cancer: "),
          sum(cv$cancer_related & cv$scope == paste0(GENE_TARGET, "-specific")))
} else log_msg("ClinVar file missing; run 00_download_data.R")

save_session("05_drg2_population_variants")
log_msg("Population variant analysis complete.")
