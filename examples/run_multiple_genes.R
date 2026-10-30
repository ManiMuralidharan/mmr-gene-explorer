#!/usr/bin/env Rscript
###############################################################################
# examples/run_multiple_genes.R
# Runs the analysis steps for several target genes in turn, reusing downloaded
# data. Run from the repository root after R/00_download_data.R:
#   Rscript examples/run_multiple_genes.R
# gnomAD/ClinVar records are gene-specific: run 00_download_data.R once per gene
# (with TARGET_GENE set) if you need script 05 for that gene.
###############################################################################

targets <- list(
  # target         partners              arm-proxy gene
  list(gene = "DRG2", partners = "MLH1,MLH3", proxy = "TP53"),
  list(gene = "ATM",  partners = "MLH1,MSH2", proxy = "ATM")
)

steps <- c("R/04_depmap_expression_msi.R",
           "R/01_tcga_cohort_definition.R",
           "R/02_tcga_target_expression_cnv_mutation.R",
           "R/03_cptac_protein_correlations.R",
           "R/05_target_population_variants.R",
           "R/06_integrate_results.R")

for (t in targets) {
  Sys.setenv(TARGET_GENE = t$gene, PARTNER_GENES = t$partners, CNA_PROXY_GENE = t$proxy)
  message("\n######## ", t$gene, " ########")
  for (s in steps) {
    res <- tryCatch({ source(s, local = new.env()); "OK" },
                    error = function(e) paste("FAILED:", conditionMessage(e)))
    message(sprintf("%-45s %s", s, res))
  }
}
