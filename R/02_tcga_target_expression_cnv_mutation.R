#!/usr/bin/env Rscript
###############################################################################
# 02_tcga_drg2_expression_cnv_mutation.R
# Questions (TCGA UCEC and COAD/READ):
#   Q1. Is target gene mRNA lower in MSI_unexplained tumours than in other groups?
#   Q2. Does target gene mRNA correlate with partner genes (default partner genes) mRNA (all tumours; MSS only)?
#   Q3. Are target gene deletions or mutations enriched in MSI_unexplained tumours?
#   Q4. Is the association independent of arm-level copy-number loss? Uses
#       CNA_PROXY_GENE (default TP53, on the same arm as DRG2: 17p) as a proxy.
#       Set CNA_PROXY_GENE to a gene on your target's chromosome arm.
###############################################################################

source(file.path("R", "00_config.R")); source(file.path("R", "utils.R"))
load_packages()

cohort <- readRDS(file.path(DIR_PROCESSED, "tcga_cohort.rds"))
genes_needed <- unique(c(GENES_FOCUS, GENES_PARTNER, GENES_MMR_EXT, CNA_PROXY_GENE))

expr_all <- list(); cna_all <- list()
for (cancer in names(TCGA_STUDIES)) {
  sdir <- file.path(DIR_RAW, "tcga", TCGA_STUDIES[[cancer]])
  # RSEM (batch-normalised) expression; log2(x + 1)
  ef <- find_file(sdir, "^data_mrna_seq_v2_rsem\\.txt$")
  em <- read_cbio_matrix(ef, genes = genes_needed)
  colnames(em) <- tcga_sample(colnames(em)); em <- em[, !duplicated(colnames(em)), drop = FALSE]
  expr_all[[cancer]] <- log2(em + 1)
  # GISTIC discrete CNA (-2 deep del, -1 shallow del, 0, 1 gain, 2 amp)
  cf <- find_file(sdir, "^data_cna\\.txt$", required = FALSE)
  if (!is.na(cf)) {
    cm <- read_cbio_matrix(cf, genes = genes_needed)
    colnames(cm) <- tcga_sample(colnames(cm)); cm <- cm[, !duplicated(colnames(cm)), drop = FALSE]
    cna_all[[cancer]] <- cm
  }
}

# ---- Long table of focus genes per sample ----------------------------------
long <- do.call(rbind, lapply(names(expr_all), function(cancer) {
  em <- expr_all[[cancer]]
  df <- data.frame(SAMPLE_ID = colnames(em), cancer = cancer, check.names = FALSE)
  for (g in intersect(genes_needed, rownames(em))) df[[paste0("expr_", g)]] <- em[g, ]
  if (!is.null(cna_all[[cancer]])) {
    cm <- cna_all[[cancer]]
    for (g in intersect(c(GENE_TARGET, CNA_PROXY_GENE, GENES_PARTNER), rownames(cm)))
      df[[paste0("cna_", g)]] <- cm[g, match(df$SAMPLE_ID, colnames(cm))]
  }
  df
}))
dat <- merge(cohort, long, by = c("SAMPLE_ID", "cancer"))
# Generic aliases for the target gene's columns
dat$expr_TARGET <- if (paste0("expr_", GENE_TARGET) %in% names(dat)) dat[[paste0("expr_", GENE_TARGET)]] else NA_real_
if (paste0("cna_", GENE_TARGET) %in% names(dat)) dat$cna_TARGET <- dat[[paste0("cna_", GENE_TARGET)]]
if (all(is.na(dat$expr_TARGET))) stop("Target gene ", GENE_TARGET, " not found in TCGA expression data.")
saveRDS(dat, file.path(DIR_PROCESSED, "tcga_focus_gene_matrix.rds"))
log_msg("Samples with expression and MMR group: ", nrow(dat))

# ---- Q1. target gene expression by MMR group --------------------------------------
main_groups <- c("MSS", "MSI_MLH1_methylated", "MSI_MMR_mutated", "MSI_unexplained")
d1 <- dat[dat$mmr_group %in% main_groups & is.finite(dat$expr_TARGET), ]
d1$mmr_group <- droplevels(d1$mmr_group)

q1 <- do.call(rbind, lapply(split(d1, d1$cancer), function(x) {
  kw <- kruskal.test(expr_TARGET ~ mmr_group, data = x)
  ref <- x$expr_TARGET[x$mmr_group == "MSI_unexplained"]
  do.call(rbind, lapply(setdiff(levels(x$mmr_group), "MSI_unexplained"), function(g) {
    other <- x$expr_TARGET[x$mmr_group == g]
    if (length(ref) < 3 || length(other) < 3) return(NULL)
    wt <- wilcox.test(ref, other, exact = FALSE)
    data.frame(cancer = x$cancer[1], comparison = paste("MSI_unexplained vs", g),
               n_unexplained = length(ref), n_other = length(other),
               median_unexplained = median(ref), median_other = median(other),
               log2_diff_median = median(ref) - median(other),
               wilcox_p = wt$p.value, kruskal_p_all_groups = kw$p.value)
  }))
}))
if (!is.null(q1)) q1$fdr <- p.adjust(q1$wilcox_p, "BH")
write.csv(q1, file.path(DIR_TABLES, paste0("02_", GENE_TARGET, "_expression_by_MMR_group.csv")), row.names = FALSE)

# Adjusted model: target gene ~ group + age (+ sex for CRC); linear model on log2 RSEM
q1_lm <- do.call(rbind, lapply(split(d1, d1$cancer), function(x) {
  x$mmr_group <- relevel(x$mmr_group, ref = "MSS")
  covs <- c("mmr_group", intersect(c("AGE"), names(x)))
  if (x$cancer[1] == "CRC" && "SEX" %in% names(x)) covs <- c(covs, "SEX")
  x$AGE <- suppressWarnings(as.numeric(x$AGE))
  fit <- lm(reformulate(covs, "expr_TARGET"), data = x)
  co <- summary(fit)$coefficients
  co <- co[grepl("^mmr_group", rownames(co)), , drop = FALSE]
  data.frame(cancer = x$cancer[1], term = rownames(co), estimate = co[, 1],
             se = co[, 2], p_value = co[, 4], covariates = paste(covs, collapse = "+"))
}))
write.csv(q1_lm, file.path(DIR_TABLES, paste0("02_", GENE_TARGET, "_expression_by_MMR_group_adjusted_lm.csv")),
          row.names = FALSE)

p1 <- ggplot(d1, aes(mmr_group, expr_TARGET, fill = mmr_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  facet_wrap(~ cancer, scales = "free_y") +
  labs(x = NULL, y = paste0(GENE_TARGET, " mRNA, log2(RSEM + 1)"), title = paste0(GENE_TARGET, " expression by MMR group")) +
  theme_mmr() + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
ggsave(file.path(DIR_FIGURES, paste0("02_", GENE_TARGET, "_expression_by_MMR_group.png")), p1, width = 10, height = 5, dpi = 300)

# ---- Q2. target gene vs partner genes (default partner genes) correlation -----------------------------------
q2 <- list()
for (cancer in names(expr_all)) {
  em <- expr_all[[cancer]]
  q2[[paste(cancer, "all")]] <- cor_table(em, GENE_TARGET, setdiff(genes_needed, GENE_TARGET),
                                          label = paste(cancer, "all tumours"))
  mss <- dat$SAMPLE_ID[dat$cancer == cancer & dat$mmr_group %in% "MSS"]
  q2[[paste(cancer, "MSS")]] <- cor_table(em[, intersect(mss, colnames(em)), drop = FALSE],
                                          GENE_TARGET, setdiff(genes_needed, GENE_TARGET),
                                          label = paste(cancer, "MSS only"))
}
q2 <- do.call(rbind, q2)
write.csv(q2, file.path(DIR_TABLES, paste0("02_", GENE_TARGET, "_mRNA_correlations.csv")), row.names = FALSE)

d2 <- dat[is.finite(dat$expr_TARGET), ]
for (g in GENES_PARTNER) {
  col <- paste0("expr_", g)
  if (!col %in% names(d2)) next
  p <- ggplot(d2, aes(expr_TARGET, .data[[col]], colour = mmr_group)) +
    geom_point(size = 1, alpha = 0.6) + geom_smooth(aes(group = 1), method = "lm", colour = "black", se = FALSE) +
    facet_wrap(~ cancer, scales = "free") +
    labs(x = paste0(GENE_TARGET, " mRNA (log2)"), y = paste(g, "mRNA (log2)"), colour = "MMR group",
         title = paste(paste0(GENE_TARGET, " vs"), g, "in TCGA tumours")) + theme_mmr()
  ggsave(file.path(DIR_FIGURES, paste0(paste0("02_", GENE_TARGET, "_vs_"), g, "_scatter.png")), p, width = 10, height = 4.5, dpi = 300)
}

# ---- Q3. target gene deletion / mutation enrichment --------------------------------
mut_rows <- list()
for (cancer in names(TCGA_STUDIES)) {
  sdir <- file.path(DIR_RAW, "tcga", TCGA_STUDIES[[cancer]])
  mut <- data.table::fread(find_file(sdir, "^data_mutations.*\\.txt$"), sep = "\t",
                           data.table = FALSE, quote = "",
                           select = c("Hugo_Symbol", "Tumor_Sample_Barcode",
                                      "Variant_Classification", "HGVSp_Short"))
  mut <- mut[mut$Hugo_Symbol == GENE_TARGET & mut$Variant_Classification != "Silent", ]
  if (nrow(mut)) { mut$SAMPLE_ID <- tcga_sample(mut$Tumor_Sample_Barcode); mut$cancer <- cancer
                   mut_rows[[cancer]] <- mut }
}
target_mut <- do.call(rbind, mut_rows)
if (!is.null(target_mut)) {
  target_mut <- merge(target_mut, cohort[, c("SAMPLE_ID", "mmr_group", "msi_mantis")],
                    by = "SAMPLE_ID", all.x = TRUE)
  write.csv(target_mut, file.path(DIR_TABLES, paste0("02_", GENE_TARGET, "_somatic_mutations.csv")), row.names = FALSE)
}

dat$target_mutated <- dat$SAMPLE_ID %in% (target_mut$SAMPLE_ID %||% character())
dat$target_deleted <- if ("cna_TARGET" %in% names(dat)) dat$cna_TARGET <= -1 else NA
dat$target_deep_del <- if ("cna_TARGET" %in% names(dat)) dat$cna_TARGET <= -2 else NA

fisher_row <- function(x, var, grp_a = "MSI_unexplained") {
  x <- x[x$mmr_group %in% main_groups & !is.na(x[[var]]), ]
  if (!nrow(x)) return(NULL)
  a <- x$mmr_group == grp_a
  tab <- table(factor(a, c(TRUE, FALSE)), factor(x[[var]], c(TRUE, FALSE)))
  ft <- fisher.test(tab)
  data.frame(cancer = x$cancer[1], feature = var,
             n_unexplained = sum(a), altered_unexplained = tab["TRUE", "TRUE"],
             n_other = sum(!a), altered_other = tab["FALSE", "TRUE"],
             odds_ratio = unname(ft$estimate), fisher_p = ft$p.value)
}
q3 <- do.call(rbind, lapply(split(dat, dat$cancer), function(x)
  do.call(rbind, lapply(c("target_mutated", "target_deleted", "target_deep_del"),
                        function(v) fisher_row(x, v)))))
write.csv(q3, file.path(DIR_TABLES, paste0("02_", GENE_TARGET, "_alteration_enrichment.csv")), row.names = FALSE)

# Expected direction check: MSI tumours are copy-number quiet, so target gene deletion
# is likely to be LESS frequent in MSI groups for reasons unrelated to target gene.
cna_by_group <- if ("cna_TARGET" %in% names(dat))
  as.data.frame(table(cancer = dat$cancer, group = dat$mmr_group, TARGET_CNA = dat$cna_TARGET)) else NULL
if (!is.null(cna_by_group))
  write.csv(cna_by_group, file.path(DIR_QC, paste0("02_", GENE_TARGET, "_CNA_by_group.csv")), row.names = FALSE)

# ---- Q4. target gene vs partner genes adjusted for arm-proxy copy number -------
q4 <- do.call(rbind, lapply(split(dat, dat$cancer), function(x) {
  proxy <- paste0("cna_", CNA_PROXY_GENE)
  outcomes <- intersect(paste0("expr_", GENES_PARTNER), names(x))
  if (is.na(CNA_PROXY_GENE) || !proxy %in% names(x) || !length(outcomes)) return(NULL)
  do.call(rbind, lapply(outcomes, function(y) {
    f <- as.formula(paste0("`", y, "` ~ expr_TARGET + factor(`", proxy, "`) + factor(mmr_group)"))
    xx <- x[x$mmr_group %in% main_groups, ]
    fit <- tryCatch(lm(f, data = xx), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    co <- summary(fit)$coefficients["expr_TARGET", ]
    data.frame(cancer = x$cancer[1], outcome = sub("expr_", "", y),
               model = paste0("outcome ~ ", GENE_TARGET, " + ", CNA_PROXY_GENE, " CNA (arm proxy) + MMR group"),
               beta_target = co[1], se = co[2], p_value = co[4], n = nobs(fit))
  }))
}))
if (!is.null(q4)) write.csv(q4, file.path(DIR_TABLES, paste0("02_", GENE_TARGET, "_partners_adjusted_for_arm_proxy.csv")), row.names = FALSE)

save_session("02_tcga_drg2_expression_cnv_mutation")
log_msg(paste0("TCGA ", GENE_TARGET, " analysis complete."))
