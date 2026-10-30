#!/usr/bin/env Rscript
###############################################################################
# 03_cptac_protein_correlations.R
# Questions (CPTAC UCEC 2020, COAD 2019):
#   Q1. Does target gene protein correlate with partner genes (default partner genes) / other MMR proteins?
#   Q2. Is target gene protein lower in MSI tumours, and in MSI tumours without MLH1
#       methylation (when annotated)?
#   Q3. Does the mRNA-level relationship hold at protein level (mRNA vs protein)?
# MLH3 is often poorly detected by mass spectrometry; missing values are
# reported, not imputed.
###############################################################################

source(file.path("R", "00_config.R")); source(file.path("R", "utils.R"))
load_packages()

prot_genes <- unique(c(GENES_FOCUS, GENES_MMR_EXT))
res_cor <- list(); res_grp <- list(); res_rna_prot <- list(); detect <- list()

for (cancer in names(CPTAC_STUDIES)) {
  sid <- CPTAC_STUDIES[[cancer]]
  sdir <- file.path(DIR_RAW, "cptac", sid)
  if (!dir.exists(sdir)) { log_msg("Missing ", sid, "; skipping"); next }
  pf <- find_file(sdir, "^data_protein_quantification.*\\.txt$", required = FALSE)
  if (is.na(pf)) { log_msg("No protein quantification in ", sid); next }
  pm <- read_cbio_matrix(pf)
  # Some CPTAC files use "SYMBOL|ID"; read_cbio_matrix already strips after '|'
  log_msg(cancer, ": protein matrix ", nrow(pm), " proteins x ", ncol(pm), " samples")

  detect[[cancer]] <- data.frame(
    cancer = cancer, gene = prot_genes,
    detected = prot_genes %in% rownames(pm),
    n_nonmissing = vapply(prot_genes, function(g)
      if (g %in% rownames(pm)) sum(is.finite(pm[g, ])) else 0L, integer(1)))

  if (!GENE_TARGET %in% rownames(pm)) { log_msg(paste0(GENE_TARGET, " protein not detected in "), sid); next }

  # Q1 protein correlations: focus genes + genome-wide rank of partner genes
  res_cor[[cancer]] <- cor_table(pm, GENE_TARGET, setdiff(prot_genes, GENE_TARGET),
                                 label = paste(cancer, "CPTAC protein"))
  gw <- cor_table(pm, GENE_TARGET, setdiff(rownames(pm), GENE_TARGET),
                  label = paste(cancer, "CPTAC protein genome-wide"))
  if (!is.null(gw)) {
    gw <- gw[order(-gw$rho), ]; gw$rank <- seq_len(nrow(gw))
    write.csv(gw, file.path(DIR_TABLES, paste0("03_", cancer, paste0("_", GENE_TARGET, "_protein_genomewide_correlations.csv"))),
              row.names = FALSE)
  }

  # Q2 MSI annotation from clinical sample file
  cs <- read_cbio_clinical(find_file(sdir, "^data_clinical_sample\\.txt$"))
  msi_col <- grep("MSI", names(cs), value = TRUE, ignore.case = TRUE)
  meth_col <- grep("MLH1.*METH|METH.*MLH1", names(cs), value = TRUE, ignore.case = TRUE)
  log_msg(cancer, " MSI columns: ", paste(msi_col, collapse = ", "),
          " | MLH1 methylation columns: ", paste(meth_col, collapse = ", "))
  if (length(msi_col)) {
    v <- as.character(cs[[msi_col[1]]])
    num <- suppressWarnings(as.numeric(v))
    msi <- if (sum(!is.na(num)) > 0.5 * length(v)) ifelse(num > MSI_MANTIS_CUTOFF, "MSI-H", "MSS")
           else ifelse(grepl("MSI-H|^MSI$|high", v, ignore.case = TRUE), "MSI-H",
                       ifelse(grepl("MSS|MSI-L|low|stable", v, ignore.case = TRUE), "MSS", NA))
    grp <- data.frame(SAMPLE_ID = cs$SAMPLE_ID, msi = msi)
    if (length(meth_col)) {
      mv <- as.character(cs[[meth_col[1]]])
      grp$mlh1_meth <- grepl("^(yes|true|methylated|1)$", mv, ignore.case = TRUE)
      grp$msi_group <- ifelse(grp$msi == "MSI-H" & grp$mlh1_meth, "MSI_MLH1_methylated",
                       ifelse(grp$msi == "MSI-H", "MSI_unmethylated", grp$msi))
    } else grp$msi_group <- grp$msi
    grp <- grp[grp$SAMPLE_ID %in% colnames(pm) & !is.na(grp$msi_group), ]
    grp$TARGET_protein <- pm[GENE_TARGET, grp$SAMPLE_ID]
    for (g in GENES_PARTNER) if (g %in% rownames(pm)) grp[[paste0(g, "_protein")]] <- pm[g, grp$SAMPLE_ID]
    grp$cancer <- cancer
    res_grp[[cancer]] <- grp

    tests <- lapply(setdiff(unique(grp$msi_group), "MSS"), function(g) {
      a <- grp$TARGET_protein[grp$msi_group == g]; b <- grp$TARGET_protein[grp$msi_group == "MSS"]
      a <- a[is.finite(a)]; b <- b[is.finite(b)]
      if (length(a) < 3 || length(b) < 3) return(NULL)
      data.frame(cancer = cancer, comparison = paste(g, "vs MSS"), n_group = length(a),
                 n_mss = length(b), median_diff = median(a) - median(b),
                 wilcox_p = wilcox.test(a, b, exact = FALSE)$p.value)
    })
    tt <- do.call(rbind, tests)
    if (!is.null(tt)) write.csv(tt, file.path(DIR_TABLES, paste0("03_", cancer, paste0("_", GENE_TARGET, "_protein_by_MSI.csv"))),
                                row.names = FALSE)
  }

  # Q3 mRNA vs protein for target gene, MLH1, MLH3
  rf <- find_file(sdir, "^data_mrna_seq.*\\.txt$", required = FALSE)
  if (!is.na(rf)) {
    rm <- read_cbio_matrix(rf, genes = c(GENE_TARGET, GENES_PARTNER))
    if (max(rm, na.rm = TRUE) > 50) rm <- log2(rm + 1)
    common <- intersect(colnames(rm), colnames(pm))
    res_rna_prot[[cancer]] <- do.call(rbind, lapply(intersect(rownames(rm), rownames(pm)), function(g) {
      x <- rm[g, common]; y <- pm[g, common]; ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < 10) return(NULL)
      ct <- cor.test(x[ok], y[ok], method = COR_METHOD, exact = FALSE)
      data.frame(cancer = cancer, gene = g, n = sum(ok), rho_mRNA_protein = unname(ct$estimate),
                 p_value = ct$p.value)
    }))
  }
}

write_if <- function(x, f) if (length(x) && !is.null(do.call(rbind, x)))
  write.csv(do.call(rbind, x), file.path(DIR_TABLES, f), row.names = FALSE)
write_if(detect, "03_CPTAC_protein_detection.csv")
write_if(res_cor, paste0("03_CPTAC_", GENE_TARGET, "_protein_correlations.csv"))
write_if(res_rna_prot, "03_CPTAC_mRNA_protein_concordance.csv")

grp_all <- do.call(rbind, lapply(res_grp, function(x) x[, intersect(c("cancer", "SAMPLE_ID", "msi_group",
                                   "TARGET_protein", paste0(GENES_PARTNER, "_protein")), names(x))]))
if (!is.null(grp_all) && nrow(grp_all)) {
  write.csv(grp_all, file.path(DIR_TABLES, paste0("03_CPTAC_", GENE_TARGET, "_protein_by_sample.csv")), row.names = FALSE)
  p <- ggplot(grp_all, aes(msi_group, TARGET_protein, fill = msi_group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) + geom_jitter(width = 0.15, size = 0.9) +
    facet_wrap(~ cancer, scales = "free_y") +
    labs(x = NULL, y = paste0(GENE_TARGET, " protein (log2 ratio)"), title = paste0("CPTAC: ", GENE_TARGET, " protein by MSI status")) +
    theme_mmr() + theme(legend.position = "none")
  ggsave(file.path(DIR_FIGURES, paste0("03_CPTAC_", GENE_TARGET, "_protein_by_MSI.png")), p, width = 9, height = 4.5, dpi = 300)
  p_col <- paste0(GENES_PARTNER[1], "_protein")
  if (p_col %in% names(grp_all)) {
    p2 <- ggplot(grp_all, aes(TARGET_protein, .data[[p_col]], colour = msi_group)) +
      geom_point(size = 1.4, alpha = 0.8) + geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "black") +
      facet_wrap(~ cancer, scales = "free") +
      labs(x = paste0(GENE_TARGET, " protein"), y = paste(GENES_PARTNER[1], "protein"), title = paste0("CPTAC: ", GENE_TARGET, " vs ", GENES_PARTNER[1], " protein")) + theme_mmr()
    ggsave(file.path(DIR_FIGURES, paste0("03_CPTAC_", GENE_TARGET, "_vs_", GENES_PARTNER[1], "_protein.png")), p2, width = 9, height = 4.5, dpi = 300)
  }
}

save_session("03_cptac_protein_correlations")
log_msg("CPTAC analysis complete.")
