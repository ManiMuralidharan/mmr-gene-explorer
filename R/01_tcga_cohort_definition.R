#!/usr/bin/env Rscript
###############################################################################
# 01_tcga_cohort_definition.R
# Classifies TCGA UCEC and COAD/READ primary tumours into MMR groups using
# PRE-SPECIFIED rules (metadata/cohort_definitions.md):
#
#   MSS                      : not MSI-H
#   MSI_MLH1_methylated      : MSI-H and MLH1 promoter methylated (sporadic)
#   MSI_MMR_mutated          : MSI-H, MLH1 unmethylated, truncating somatic
#                              mutation in a core MMR gene
#   MSI_unexplained          : MSI-H, MLH1 unmethylated, no truncating core MMR
#                              mutation  ("Lynch-like" proxy; germline unknown)
#   MSI_indeterminate        : MSI-H, MLH1 methylation unavailable
#
# NOTE: Open TCGA data contain somatic calls only. "MSI_unexplained" is a
# tumour-level proxy for Lynch-like syndrome; germline confirmation requires
# dbGaP-controlled data (data/controlled/).
###############################################################################

source(file.path("R", "00_config.R")); source(file.path("R", "utils.R"))
load_packages()

classify_study <- function(cancer, study_id) {
  sdir <- file.path(DIR_RAW, "tcga", study_id)
  stopifnot(dir.exists(sdir))
  log_msg("Classifying ", cancer, " (", study_id, ")")

  # ---- Clinical / MSI scores ----
  cs <- read_cbio_clinical(find_file(sdir, "^data_clinical_sample\\.txt$"))
  cp <- read_cbio_clinical(find_file(sdir, "^data_clinical_patient\\.txt$"))
  cs$SAMPLE_ID <- tcga_sample(cs$SAMPLE_ID)
  cs <- cs[substr(cs$SAMPLE_ID, 14, 15) == "01", , drop = FALSE]   # primary tumours

  # PanCancer molecular subtype lives in the PATIENT file (e.g. UCEC_MSI, COAD_POLE)
  cp_sub <- cp; cp_sub$PATIENT_ID <- tcga_patient(cp_sub$PATIENT_ID)
  if ("SUBTYPE" %in% names(cp_sub) && !"SUBTYPE" %in% names(cs))
    cs$SUBTYPE <- cp_sub$SUBTYPE[match(tcga_patient(cs$SAMPLE_ID), cp_sub$PATIENT_ID)]

  num <- function(x) suppressWarnings(as.numeric(x))
  mantis <- if ("MSI_SCORE_MANTIS" %in% names(cs)) num(cs$MSI_SCORE_MANTIS) else NA
  sensor <- if ("MSI_SENSOR_SCORE" %in% names(cs)) num(cs$MSI_SENSOR_SCORE) else NA
  msi_h <- ifelse(!is.na(mantis), mantis > MSI_MANTIS_CUTOFF,
                  ifelse(!is.na(sensor), sensor >= MSI_SENSOR_CUTOFF, NA))
  # Fallback to molecular subtype when scores are missing
  if ("SUBTYPE" %in% names(cs)) {
    # PanCancer subtypes: UCEC_MSI (endometrial), GI.HM-indel (MSI colorectal)
    sub_msi <- grepl("MSI|HM-indel", cs$SUBTYPE, ignore.case = TRUE)
    msi_h[is.na(msi_h) & !is.na(cs$SUBTYPE)] <- sub_msi[is.na(msi_h) & !is.na(cs$SUBTYPE)]
  }

  # ---- POLE-ultramutated tumours (excluded from MSI groups; separate label) ----
  pole_ultra <- if ("SUBTYPE" %in% names(cs)) (!is.na(cs$SUBTYPE) & grepl("POLE|HM-SNV", cs$SUBTYPE, ignore.case = TRUE)) else rep(FALSE, nrow(cs))

  # ---- MLH1 promoter methylation ----
  # cBioPortal methylation files are PROBE-level. MLH1 promoter score = mean beta of all
  # HM450 probes annotated to MLH1 with TSS200 in their annotation (core promoter).
  # Selection uses annotation only (not MSI or target gene). Upstream TSS1500-only probes
  # (e.g. cg10990993, cg18320188) are constitutively methylated and were excluded.
  meth_file <- find_file(sdir, "^data_methylation_hm450\\.txt$", required = FALSE)
  if (is.na(meth_file)) meth_file <- find_file(sdir, "^data_methylation.*\\.txt$", required = FALSE)
  mlh1_beta <- setNames(rep(NA_real_, nrow(cs)), cs$SAMPLE_ID)
  if (!is.na(meth_file)) {
    md <- data.table::fread(meth_file, sep = "\t", data.table = FALSE, check.names = FALSE)
    is_mlh1 <- vapply(strsplit(as.character(md$NAME), ";"), function(g) "MLH1" %in% g, logical(1))
    is_core <- grepl("TSS200", md$DESCRIPTION)
    md <- md[is_mlh1 & is_core, , drop = FALSE]
    log_msg(study_id, ": using ", nrow(md), " MLH1 TSS200 probes: ", paste(md$ENTITY_STABLE_ID, collapse = ","))
    if (nrow(md)) {
      m <- as.matrix(md[, setdiff(names(md), c("ENTITY_STABLE_ID", "NAME", "DESCRIPTION", "TRANSCRIPT_ID")), drop = FALSE])
      suppressWarnings(storage.mode(m) <- "numeric")
      colnames(m) <- tcga_sample(colnames(m))
      score <- colMeans(m, na.rm = TRUE); score[is.nan(score)] <- NA
      score <- score[!duplicated(names(score))]
      common <- intersect(names(mlh1_beta), names(score))
      mlh1_beta[common] <- score[common]
      log_msg(study_id, ": MLH1 promoter score for ", sum(!is.na(mlh1_beta)), " of ", length(mlh1_beta), " tumours")
      write.csv(data.frame(probe = md$ENTITY_STABLE_ID, annotation = md$DESCRIPTION),
                file.path(DIR_QC, paste0("01_MLH1_probes_used_", cancer, ".csv")), row.names = FALSE)
    }
  } else log_msg("No methylation file for ", study_id)

  # Fallback for HM27-only samples: TSS200 MLH1 probe(s) present in the merged HM27/HM450 file
  merged_file <- find_file(sdir, "^data_methylation_hm27_hm450_merged\\.txt$", required = FALSE)
  missing_ids <- names(mlh1_beta)[is.na(mlh1_beta)]
  if (!is.na(merged_file) && length(missing_ids)) {
    mg <- data.table::fread(merged_file, sep = "\t", data.table = FALSE, check.names = FALSE)
    is_mlh1 <- vapply(strsplit(as.character(mg$NAME), ";"), function(g) "MLH1" %in% g, logical(1))
    mg <- mg[is_mlh1 & grepl("TSS200", mg$DESCRIPTION), , drop = FALSE]
    if (nrow(mg)) {
      m <- as.matrix(mg[, setdiff(names(mg), c("ENTITY_STABLE_ID", "NAME", "DESCRIPTION", "TRANSCRIPT_ID")), drop = FALSE])
      suppressWarnings(storage.mode(m) <- "numeric")
      colnames(m) <- tcga_sample(colnames(m))
      sc <- colMeans(m, na.rm = TRUE); sc[is.nan(sc)] <- NA; sc <- sc[!duplicated(names(sc))]
      fill <- intersect(missing_ids, names(sc)[!is.na(sc)])
      mlh1_beta[fill] <- sc[fill]
      log_msg(study_id, ": fallback probe(s) ", paste(mg$ENTITY_STABLE_ID, collapse = ","),
              " filled ", length(fill), " HM27-only tumours")
    }
  }

  # ---- Somatic mutations in core MMR genes and target gene ----
  mut <- data.table::fread(find_file(sdir, "^data_mutations.*\\.txt$"), sep = "\t",
                           data.table = FALSE, quote = "",
                           select = c("Hugo_Symbol", "Tumor_Sample_Barcode",
                                      "Variant_Classification", "HGVSp_Short"))
  mut$SAMPLE_ID <- tcga_sample(mut$Tumor_Sample_Barcode)
  mmr_trunc <- unique(mut$SAMPLE_ID[mut$Hugo_Symbol %in% GENES_MMR_CORE &
                                      mut$Variant_Classification %in% TRUNCATING_CLASSES])
  mmr_any <- unique(mut$SAMPLE_ID[mut$Hugo_Symbol %in% GENES_MMR_CORE &
                                    mut$Variant_Classification != "Silent"])
  sequenced <- unique(mut$SAMPLE_ID)

  # ---- Assign groups ----
  meth <- mlh1_beta[cs$SAMPLE_ID] > MLH1_METH_CUTOFF
  grp <- rep(NA_character_, nrow(cs))
  grp[pole_ultra] <- "POLE_ultramutated"
  idx <- !pole_ultra
  grp[idx & msi_h %in% FALSE] <- "MSS"
  grp[idx & msi_h %in% TRUE & meth %in% TRUE] <- "MSI_MLH1_methylated"
  grp[idx & msi_h %in% TRUE & is.na(meth)] <- "MSI_indeterminate"
  grp[idx & msi_h %in% TRUE & meth %in% FALSE &
        cs$SAMPLE_ID %in% mmr_trunc] <- "MSI_MMR_mutated"
  grp[idx & msi_h %in% TRUE & meth %in% FALSE &
        !cs$SAMPLE_ID %in% mmr_trunc & cs$SAMPLE_ID %in% sequenced] <- "MSI_unexplained"

  out <- data.frame(
    cancer = cancer, study_id = study_id,
    SAMPLE_ID = cs$SAMPLE_ID, PATIENT_ID = tcga_patient(cs$SAMPLE_ID),
    msi_mantis = mantis, msi_sensor = sensor, msi_high = msi_h,
    subtype = if ("SUBTYPE" %in% names(cs)) cs$SUBTYPE else NA,
    mlh1_promoter_beta = unname(mlh1_beta[cs$SAMPLE_ID]),
    mmr_truncating_somatic = cs$SAMPLE_ID %in% mmr_trunc,
    mmr_any_nonsilent_somatic = cs$SAMPLE_ID %in% mmr_any,
    mutation_profiled = cs$SAMPLE_ID %in% sequenced,
    mmr_group = grp
  )

  # Patient covariates for later adjustment
  cp$PATIENT_ID <- tcga_patient(cp$PATIENT_ID)
  keep_cp <- intersect(c("PATIENT_ID", "AGE", "SEX", "AJCC_PATHOLOGIC_TUMOR_STAGE",
                         "RACE", "OS_MONTHS", "OS_STATUS"), names(cp))
  out <- merge(out, cp[, keep_cp, drop = FALSE], by = "PATIENT_ID", all.x = TRUE)
  out
}

cohort <- do.call(rbind, Map(classify_study, names(TCGA_STUDIES), TCGA_STUDIES))
cohort$mmr_group <- factor(cohort$mmr_group,
  levels = c("MSS", "MSI_MLH1_methylated", "MSI_MMR_mutated", "MSI_unexplained",
             "MSI_indeterminate", "POLE_ultramutated"))

saveRDS(cohort, file.path(DIR_PROCESSED, "tcga_cohort.rds"))
write.csv(cohort, file.path(DIR_TABLES, "01_tcga_cohort_mmr_groups.csv"), row.names = FALSE)

summ <- as.data.frame.matrix(table(cohort$cancer, cohort$mmr_group, useNA = "ifany"))
summ$cancer <- rownames(summ)
write.csv(summ, file.path(DIR_TABLES, "01_tcga_mmr_group_counts.csv"), row.names = FALSE)
print(summ)

# QC: MSI score vs MLH1 methylation
qc <- cohort[!is.na(cohort$msi_mantis) & !is.na(cohort$mlh1_promoter_beta), ]
if (nrow(qc)) {
  p <- ggplot(qc, aes(msi_mantis, mlh1_promoter_beta, colour = mmr_group)) +
    geom_point(alpha = 0.7, size = 1.6) +
    geom_vline(xintercept = MSI_MANTIS_CUTOFF, linetype = 2) +
    geom_hline(yintercept = MLH1_METH_CUTOFF, linetype = 2) +
    facet_wrap(~ cancer) +
    labs(x = "MSI score (MANTIS)", y = "MLH1 promoter methylation (beta)",
         colour = "MMR group", title = "TCGA MMR group assignment") +
    theme_mmr()
  ggsave(file.path(DIR_QC, "01_msi_vs_mlh1_methylation.png"), p, width = 10, height = 4.5, dpi = 300)
}

save_session("01_tcga_cohort_definition")
log_msg("Cohort definition complete.")
