#!/usr/bin/env Rscript
###############################################################################
# 06_integrate_results.R
# Collects the key results of scripts 01 to 05 into:
#   results/tables/06_evidence_summary.csv   one row per question
#   reports/target gene_database_summary.md          short summary for collaborators
# Each line of evidence is labelled Supportive / Not supportive / Inconclusive
# using pre-specified rules (direction + FDR), so the summary cannot be
# selectively worded after seeing the data.
###############################################################################

source(file.path("R", "00_config.R")); source(file.path("R", "utils.R"))
load_packages()

rd <- function(f) { p <- file.path(DIR_TABLES, f); if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else NULL }
fmt <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, digits = d, format = "g"))
call_ev <- function(effect, p, want_sign, alpha = FDR_CUTOFF) {
  if (is.null(effect) || length(effect) == 0 || is.na(effect) || is.na(p)) return("Not available")
  if (p < alpha && sign(effect) == want_sign) return("Supportive")
  if (p < alpha && sign(effect) != want_sign) return("Contradicts")
  "Inconclusive"
}

ev <- list()
add <- function(source, question, estimate, p, verdict, note = "")
  ev[[length(ev) + 1]] <<- data.frame(source = source, question = question,
                                      estimate = estimate, p_or_fdr = p,
                                      verdict = verdict, note = note)

# ---- DepMap ----------------------------------------------------------------
dm <- rd(paste0("04_DepMap_", GENE_TARGET, "_mRNA_correlations.csv"))
if (!is.null(dm)) for (i in which(dm$gene_y %in% GENES_PARTNER)) {
  r <- dm[i, ]
  add("DepMap", paste0(paste0(GENE_TARGET, " vs "), r$gene_y, " mRNA (", r$dataset, ")"),
      paste0("rho=", fmt(r$rho), ", n=", r$n), fmt(r$fdr),
      call_ev(r$rho, r$fdr, +1))
}
dmr <- rd("04_DepMap_partner_rank_vs_random_genes.csv")
if (!is.null(dmr)) for (i in seq_len(nrow(dmr)))
  add("DepMap", paste("Genome-wide rank of", dmr$gene[i]),
      paste0("rho=", fmt(dmr$rho[i]), ", |rho| exceeds ", round(100 * dmr$empirical_percentile_abs[i]), "% of genes"),
      NA, "Descriptive")
dmsi <- rd(paste0("04_DepMap_", GENE_TARGET, "_vs_MSI.csv"))
if (!is.null(dmsi)) for (i in seq_len(nrow(dmsi)))
  add("DepMap", paste0(paste0(GENE_TARGET, " lower in MSI-H lines ("), dmsi$subset[i], ")"),
      paste0("median diff=", fmt(dmsi$median_diff_target[i]), "; n=", dmsi$n_msi[i], "/", dmsi$n_mss[i]),
      fmt(dmsi$wilcox_p[i]), call_ev(dmsi$median_diff_target[i], dmsi$wilcox_p[i], -1),
      "nominal P (3 subsets)")

# ---- TCGA ------------------------------------------------------------------
t1 <- rd(paste0("02_", GENE_TARGET, "_expression_by_MMR_group.csv"))
if (!is.null(t1)) for (i in seq_len(nrow(t1)))
  add("TCGA", paste0(t1$cancer[i], paste0(": ", GENE_TARGET, " mRNA, "), t1$comparison[i]),
      paste0("log2 diff=", fmt(t1$log2_diff_median[i]), "; n=", t1$n_unexplained[i], "/", t1$n_other[i]),
      fmt(t1$fdr[i]), call_ev(t1$log2_diff_median[i], t1$fdr[i], -1))
t2 <- rd(paste0("02_", GENE_TARGET, "_mRNA_correlations.csv"))
if (!is.null(t2)) for (i in which(t2$gene_y %in% GENES_PARTNER))
  add("TCGA", paste0(paste0(GENE_TARGET, " vs "), t2$gene_y[i], " mRNA (", t2$dataset[i], ")"),
      paste0("rho=", fmt(t2$rho[i]), ", n=", t2$n[i]), fmt(t2$fdr[i]),
      call_ev(t2$rho[i], t2$fdr[i], +1))
t3 <- rd(paste0("02_", GENE_TARGET, "_alteration_enrichment.csv"))
if (!is.null(t3)) for (i in seq_len(nrow(t3)))
  add("TCGA", paste0(t3$cancer[i], ": ", t3$feature[i], " enriched in MSI_unexplained"),
      paste0("OR=", fmt(t3$odds_ratio[i]), "; ", t3$altered_unexplained[i], "/", t3$n_unexplained[i],
             " vs ", t3$altered_other[i], "/", t3$n_other[i]),
      fmt(t3$fisher_p[i]), call_ev(log(t3$odds_ratio[i]), t3$fisher_p[i], +1),
      "Deletions expected to be rarer in MSI tumours (copy-number quiet)")
t4 <- rd(paste0("02_", GENE_TARGET, "_partners_adjusted_for_arm_proxy.csv"))
if (!is.null(t4)) for (i in seq_len(nrow(t4)))
  add("TCGA", paste0(t4$cancer[i], paste0(": ", GENE_TARGET, " to "), t4$outcome[i], paste0(" adjusted for ", CNA_PROXY_GENE, " CNA (arm proxy) + MMR group")),
      paste0("beta=", fmt(t4$beta_target[i]), ", n=", t4$n[i]), fmt(t4$p_value[i]),
      call_ev(t4$beta_target[i], t4$p_value[i], +1), "nominal P")

# ---- CPTAC -----------------------------------------------------------------
c1 <- rd(paste0("03_CPTAC_", GENE_TARGET, "_protein_correlations.csv"))
if (!is.null(c1)) for (i in which(c1$gene_y %in% GENES_PARTNER))
  add("CPTAC", paste0(paste0(GENE_TARGET, " vs "), c1$gene_y[i], " protein (", c1$dataset[i], ")"),
      paste0("rho=", fmt(c1$rho[i]), ", n=", c1$n[i]), fmt(c1$fdr[i]),
      call_ev(c1$rho[i], c1$fdr[i], +1))
for (f in list.files(DIR_TABLES, paste0("^03_.*_", GENE_TARGET, "_protein_by_MSI\\.csv$"))) {
  x <- rd(f)
  for (i in seq_len(nrow(x)))
    add("CPTAC", paste0(x$cancer[i], paste0(": ", GENE_TARGET, " protein, "), x$comparison[i]),
        paste0("median diff=", fmt(x$median_diff[i]), "; n=", x$n_group[i], "/", x$n_mss[i]),
        fmt(x$wilcox_p[i]), call_ev(x$median_diff[i], x$wilcox_p[i], -1), "nominal P")
}

# ---- gnomAD / ClinVar ------------------------------------------------------
g1 <- rd(paste0("05_gnomAD_", GENE_TARGET, "_constraint.csv"))
if (!is.null(g1))
  add("gnomAD", paste0(GENE_TARGET, " loss-of-function constraint"),
      paste0("pLI=", fmt(g1$pli[1]), ", LOEUF=", fmt(g1$oe_lof_upper[1])), NA, g1$interpretation[1],
      "Constraint argues for importance; tolerance argues against dominant predisposition")
g2 <- rd(paste0("05_gnomAD_", GENE_TARGET, "_LoF_summary.csv"))
if (!is.null(g2))
  add("gnomAD", "High-confidence LoF carriers in the population",
      paste0(g2$hc_lof_allele_count[1], " alleles, ", g2$hc_lof_homozygotes[1],
             " homozygotes; carrier freq ~", fmt(g2$approx_carrier_frequency[1])), NA, "Descriptive")
cv <- rd(paste0("05_ClinVar_", GENE_TARGET, "_annotated.csv"))
if (!is.null(cv))
  add("ClinVar", paste0(GENE_TARGET, "-specific records linked to cancer"),
      paste0(sum(cv$cancer_related & cv$scope == paste0(GENE_TARGET, "-specific")), " of ",
             sum(cv$scope == paste0(GENE_TARGET, "-specific")), paste0(" ", GENE_TARGET, "-specific records")), NA, "Descriptive")

evidence <- do.call(rbind, ev)
if (is.null(evidence)) stop("No result tables found. Run scripts 01 to 05 first.")
write.csv(evidence, file.path(DIR_TABLES, "06_evidence_summary.csv"), row.names = FALSE)

# ---- Markdown summary ------------------------------------------------------
cohort_counts <- rd("01_tcga_mmr_group_counts.csv")
md <- c(
  paste0("# ", GENE_TARGET, " and mismatch repair in human cancer: database summary"),
  "",
  paste("Generated:", Sys.Date()),
  "",
  "## Hypothesis tested",
  paste0("Lower ", GENE_TARGET, " is associated with lower ", paste(GENES_PARTNER, collapse = "/"),
         " and with mismatch repair-deficient tumours that are not explained by MLH1 promoter"),
  "methylation or truncating MMR mutations (a tumour-level proxy for Lynch-like syndrome).",
  "",
  "## TCGA cohort (pre-specified MMR groups)",
  if (!is.null(cohort_counts)) c("", paste(capture.output(print(cohort_counts, row.names = FALSE)), collapse = "\n")) else "Not run.",
  "",
  "## Evidence table",
  "",
  "| Source | Question | Estimate | P / FDR | Verdict |",
  "|---|---|---|---|---|",
  sprintf("| %s | %s | %s | %s | %s |", evidence$source, evidence$question,
          evidence$estimate, evidence$p_or_fdr, evidence$verdict),
  "",
  "## Verdict counts",
  paste0("- ", names(table(evidence$verdict)), ": ", as.integer(table(evidence$verdict))),
  "",
  "## Limitations",
  "- Open TCGA data are somatic only; 'MSI_unexplained' is a tumour-level proxy for Lynch-like",
  "  syndrome. Germline confirmation needs dbGaP-controlled data or clinical cohorts.",
  paste0("- Cross-sectional correlations cannot establish that ", GENE_TARGET, " regulates ",
         paste(GENES_PARTNER, collapse = "/"), "."),
  paste0("- Copy number of ", GENE_TARGET, " can be confounded by arm-level loss; adjusted using ",
         CNA_PROXY_GENE, " copy number as an arm proxy (set CNA_PROXY_GENE in 00_config.R)."),
  "- Proteomics may not detect low-abundance MMR proteins (see 03_CPTAC_protein_detection.csv).",
  ""
)
dir.create(file.path(PROJECT_DIR, "reports"), showWarnings = FALSE)
writeLines(md, file.path(PROJECT_DIR, "reports", paste0(GENE_TARGET, "_database_summary.md")))

print(table(evidence$verdict))
save_session("06_integrate_results")
log_msg(paste0("Integration complete: reports/", GENE_TARGET, "_database_summary.md"))
