#!/usr/bin/env Rscript
###############################################################################
# 04_depmap_expression_msi.R
# Quickest first check (run this first).
# Questions (DepMap / CCLE cell lines):
#   Q1. Does target gene mRNA correlate with partner genes (default partner genes) mRNA across all lines, and
#       within colorectal and endometrial lines?
#   Q2. Where do MLH1 and MLH3 rank among all genes correlated with target gene?
#   Q3. Is target gene lower in MSI cell lines (MSIScore, if OmicsSignatures.csv given)?
#   Q4. Optional: CRISPR co-dependency of target gene with MMR genes.
#   Q5. Optional: CCLE proteomics target gene vs partner genes (default partner genes) protein.
# Files are read by name from data/raw/depmap/ (manual download).
###############################################################################

source(file.path("R", "00_config.R")); source(file.path("R", "utils.R"))
load_packages()

dm <- file.path(DIR_RAW, "depmap")
f_model <- find_file(dm, "^Model\\.csv$")
# 26Q3+: OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv ; older: OmicsExpressionProteinCodingGenesTPMLogp1.csv
f_expr  <- find_file(dm, "^OmicsExpression(TPMLogp1HumanProteinCodingGenes|ProteinCodingGenesTPMLogp1).*\\.csv$")
# 26Q3+: OmicsGlobalSignatures.csv ; older: OmicsSignatures.csv
f_sig   <- find_file(dm, "^Omics(Global)?Signatures.*\\.csv$", required = FALSE)
f_crispr <- find_file(dm, "^CRISPRGeneEffect.*\\.csv$", required = FALSE)
f_prot  <- find_file(dm, "MS_CCLE|proteomics", required = FALSE)

# DepMap matrices: rows = ModelID, columns = "SYMBOL (EntrezID)"
read_depmap_matrix <- function(path) {
  dt <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  id_col <- names(dt)[1]
  # Newer releases may carry extra annotation columns; keep numeric gene columns only
  meta_cols <- intersect(c("SequencingID", "ProfileID", "ModelID", "IsDefaultEntryForModel",
                           "is_default_entry", "ModelConditionID", "IsDefaultEntryForMC"), names(dt))
  def_col <- intersect(c("IsDefaultEntryForModel", "is_default_entry"), names(dt))[1]
  if (!is.na(def_col))
    dt <- dt[dt[[def_col]] %in% c(TRUE, "Yes", "True", "TRUE", 1, "1"), , drop = FALSE]
  ids <- if ("ModelID" %in% names(dt)) dt$ModelID else dt[[id_col]]
  gene_cols <- setdiff(names(dt), c(id_col, meta_cols))
  m <- as.matrix(dt[, gene_cols, drop = FALSE]); suppressWarnings(storage.mode(m) <- "numeric")
  colnames(m) <- sub(" \\(.*\\)$", "", colnames(m))
  rownames(m) <- ids
  m <- m[!duplicated(rownames(m)), !duplicated(colnames(m)), drop = FALSE]
  t(m)   # genes x models
}

model <- data.table::fread(f_model, data.table = FALSE)
lineage_col <- intersect(c("OncotreeLineage", "lineage"), names(model))[1]
subtype_col <- intersect(c("OncotreePrimaryDisease", "OncotreeSubtype", "primary_disease"), names(model))[1]
model$lineage <- model[[lineage_col]]
model$disease <- model[[subtype_col]]
model$group <- "Other"
for (nm in names(DEPMAP_SUBGROUPS)) {
  hit <- grepl(DEPMAP_SUBGROUPS[[nm]], model$lineage, ignore.case = TRUE) |
         grepl(DEPMAP_SUBGROUPS[[nm]], model$disease, ignore.case = TRUE)
  model$group[hit & model$group == "Other"] <- nm
}

expr_ids <- tryCatch(data.table::fread(f_expr, select = intersect(c("ProfileID", "ModelID"),
                       names(data.table::fread(f_expr, nrows = 0))), data.table = FALSE),
                     error = function(e) data.frame())
expr <- read_depmap_matrix(f_expr)
log_msg("Expression: ", nrow(expr), " genes x ", ncol(expr), " models")
stopifnot(GENE_TARGET %in% rownames(expr))

# ---- Q1. Focus correlations ------------------------------------------------
others <- setdiff(unique(c(GENES_FOCUS, GENES_MMR_EXT)), GENE_TARGET)
q1 <- list(cor_table(expr, GENE_TARGET, others, label = "All cell lines"))
for (g in names(DEPMAP_SUBGROUPS)) {
  ids <- intersect(model$ModelID[model$group == g], colnames(expr))
  if (length(ids) >= 10)
    q1[[g]] <- cor_table(expr[, ids, drop = FALSE], GENE_TARGET, others, label = paste(g, "lines"))
}
q1 <- do.call(rbind, q1)
write.csv(q1, file.path(DIR_TABLES, paste0("04_DepMap_", GENE_TARGET, "_mRNA_correlations.csv")), row.names = FALSE)
print(q1[q1$gene_y %in% GENES_PARTNER, ])

# ---- Q2. Genome-wide rank of partner genes (default partner genes) -----------------------------------
x <- expr[GENE_TARGET, ]
keep <- rowSums(is.finite(expr)) >= 100 & apply(expr, 1, sd, na.rm = TRUE) > 0
rho_all <- apply(expr[keep, , drop = FALSE], 1, function(y)
  suppressWarnings(cor(x, y, method = COR_METHOD, use = "pairwise.complete.obs")))
rho_all <- sort(rho_all[names(rho_all) != GENE_TARGET], decreasing = TRUE)
gw <- data.frame(gene = names(rho_all), rho = unname(rho_all), rank = seq_along(rho_all),
                 percentile = round(100 * seq_along(rho_all) / length(rho_all), 2))
write.csv(gw, file.path(DIR_TABLES, paste0("04_DepMap_", GENE_TARGET, "_genomewide_correlations.csv")), row.names = FALSE)
print(gw[gw$gene %in% unique(c(GENES_PARTNER, GENES_MMR_CORE)), ])

# Null comparison: is |rho(target gene, MLH1)| unusual vs. random genes?
set.seed(SEED)
null_rho <- sample(rho_all, min(5000, length(rho_all)))
null_summary <- data.frame(
  gene = GENES_PARTNER,
  rho = unname(rho_all[GENES_PARTNER]),
  empirical_percentile_abs = sapply(GENES_PARTNER, function(g)
    if (g %in% names(rho_all)) mean(abs(null_rho) <= abs(rho_all[g])) else NA))
write.csv(null_summary, file.path(DIR_TABLES, "04_DepMap_partner_rank_vs_random_genes.csv"), row.names = FALSE)

# ---- Q3. MSI association ---------------------------------------------------
if (!is.na(f_sig)) {
  sig <- data.table::fread(f_sig, data.table = FALSE)
  id <- intersect(c("ModelID", names(sig)[1]), names(sig))[1]
  msi_col <- grep("^MSI", names(sig), value = TRUE)[1]
  if (!is.na(msi_col)) {
    ids <- sig[[id]]
    # Some releases index signatures by ProfileID (PR-*); map to ModelID via the expression file
    if (any(grepl("^PR-", ids)) && "ProfileID" %in% names(expr_ids)) {
      ids <- expr_ids$ModelID[match(ids, expr_ids$ProfileID)]
    }
    d <- data.frame(ModelID = ids, MSIScore = suppressWarnings(as.numeric(sig[[msi_col]])))
    d <- d[!is.na(d$ModelID) & !duplicated(d$ModelID), ]
    d <- d[d$ModelID %in% colnames(expr) & is.finite(d$MSIScore), ]
    d$TARGET <- expr[GENE_TARGET, d$ModelID]
    for (g in GENES_PARTNER) if (g %in% rownames(expr)) d[[g]] <- expr[g, d$ModelID]
    d <- merge(d, model[, c("ModelID", "group")], by = "ModelID")
    # MSI-high = top of MSIScore distribution (DepMap MSIScore ~ % unstable sites);
    # threshold recorded in output for transparency
    # DepMap defines MSI as MSIsensor2 score >= 20 (OmicsInferredMolecularSubtypes)
    thr <- DEPMAP_MSI_CUTOFF
    d$msi_status <- ifelse(d$MSIScore >= thr, "MSI-H", "MSS")
    tests <- do.call(rbind, lapply(c("All", names(DEPMAP_SUBGROUPS)), function(g) {
      dd <- if (g == "All") d else d[d$group == g, ]
      a <- dd$TARGET[dd$msi_status == "MSI-H"]; b <- dd$TARGET[dd$msi_status == "MSS"]
      if (length(a) < 3 || length(b) < 3) return(NULL)
      ct <- cor.test(dd$TARGET, dd$MSIScore, method = COR_METHOD, exact = FALSE)
      data.frame(subset = g, msi_threshold = thr, n_msi = length(a), n_mss = length(b),
                 median_diff_target = median(a) - median(b),
                 wilcox_p = wilcox.test(a, b, exact = FALSE)$p.value,
                 rho_target_MSIScore = unname(ct$estimate), rho_p = ct$p.value)
    }))
    write.csv(tests, file.path(DIR_TABLES, paste0("04_DepMap_", GENE_TARGET, "_vs_MSI.csv")), row.names = FALSE)
    write.csv(d, file.path(DIR_PROCESSED, paste0("depmap_", GENE_TARGET, "_MSI_per_model.csv")), row.names = FALSE)
    p <- ggplot(d, aes(msi_status, TARGET, fill = msi_status)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) + geom_jitter(width = 0.15, size = 0.6, alpha = 0.5) +
      facet_wrap(~ group) + labs(x = NULL, y = paste0(GENE_TARGET, " mRNA log2(TPM + 1)"),
                                 title = paste0("DepMap: ", GENE_TARGET, " by MSI status")) +
      theme_mmr() + theme(legend.position = "none")
    ggsave(file.path(DIR_FIGURES, paste0("04_DepMap_", GENE_TARGET, "_by_MSI.png")), p, width = 9, height = 4.5, dpi = 300)
  }
}

# Scatter target gene vs partner genes (default partner genes)
sc <- data.frame(ModelID = colnames(expr), TARGET = expr[GENE_TARGET, ])
for (g in GENES_PARTNER) if (g %in% rownames(expr)) sc[[g]] <- expr[g, ]
sc <- merge(sc, model[, c("ModelID", "group")], by = "ModelID")
for (g in intersect(GENES_PARTNER, names(sc))) {
  p <- ggplot(sc, aes(TARGET, .data[[g]])) + geom_point(aes(colour = group), size = 0.8, alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, colour = "black") +
    labs(x = paste0(GENE_TARGET, " mRNA log2(TPM+1)"), y = paste(g, "mRNA log2(TPM+1)"),
         title = paste(paste0("DepMap: ", GENE_TARGET, " vs"), g)) + theme_mmr()
  ggsave(file.path(DIR_FIGURES, paste0(paste0("04_DepMap_", GENE_TARGET, "_vs_"), g, ".png")), p, width = 6.5, height = 5, dpi = 300)
}

# ---- Q4. CRISPR co-dependency (optional) -----------------------------------
if (!is.na(f_crispr)) {
  cr <- read_depmap_matrix(f_crispr)
  if (GENE_TARGET %in% rownames(cr)) {
    q4 <- cor_table(cr, GENE_TARGET, others, method = "pearson", label = "CRISPR gene effect")
    write.csv(q4, file.path(DIR_TABLES, paste0("04_DepMap_", GENE_TARGET, "_CRISPR_codependency.csv")), row.names = FALSE)
  }
}

# ---- Q5. CCLE proteomics (optional) ----------------------------------------
if (!is.na(f_prot)) {
  pr <- tryCatch({
    dt <- data.table::fread(f_prot, data.table = FALSE, check.names = FALSE)
    # Gygi harmonized file: rows = models, columns = proteins (symbol in name)
    ids <- dt[[1]]; m <- as.matrix(dt[, -1]); suppressWarnings(storage.mode(m) <- "numeric")
    colnames(m) <- sub("^([A-Za-z0-9-]+).*$", "\\1", colnames(m)); rownames(m) <- ids
    t(m[, !duplicated(colnames(m)), drop = FALSE])
  }, error = function(e) { log_msg("Could not parse proteomics file: ", conditionMessage(e)); NULL })
  if (!is.null(pr) && GENE_TARGET %in% rownames(pr)) {
    q5 <- cor_table(pr, GENE_TARGET, others, label = "CCLE proteomics")
    write.csv(q5, file.path(DIR_TABLES, paste0("04_DepMap_", GENE_TARGET, "_protein_correlations.csv")), row.names = FALSE)
  }
}

save_session("04_depmap_expression_msi")
log_msg("DepMap analysis complete.")
