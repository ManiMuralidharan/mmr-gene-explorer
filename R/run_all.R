#!/usr/bin/env Rscript
###############################################################################
# run_all.R : runs the full pipeline from the project root.
#   Rscript R/run_all.R
# Order: DepMap first (fastest sanity check), then TCGA, CPTAC, population
# variants, integration. Each step is isolated so one failure does not stop
# the others; failures are listed at the end.
###############################################################################

steps <- c(
  "R/00_download_data.R",
  "R/04_depmap_expression_msi.R",
  "R/01_tcga_cohort_definition.R",
  "R/02_tcga_target_expression_cnv_mutation.R",
  "R/03_cptac_protein_correlations.R",
  "R/05_target_population_variants.R",
  "R/06_integrate_results.R"
)

status <- character()
for (s in steps) {
  message("\n==================== ", s, " ====================")
  ok <- tryCatch({ source(s, local = new.env()); "OK" },
                 error = function(e) paste("FAILED:", conditionMessage(e)))
  status[s] <- ok
}
message("\n==================== Summary ====================")
print(data.frame(step = names(status), status = unname(status)), row.names = FALSE)
