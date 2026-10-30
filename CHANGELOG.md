# Changelog

## 1.1.0 (2026-10-30)
- Target gene, partner genes, arm-proxy gene, studies and DepMap subgroups are configurable
  (`R/00_config.R` or environment variables). Results are written per target gene.
- MLH1 promoter methylation read from probe-level cBioPortal files (HM450 TSS200 probes, with HM27 fallback).
- PanCancer SUBTYPE read from the patient file.
- cBioPortal DataHub URL updated to `datahub.assets.cbioportal.org`.
- DepMap 26Q3 file names supported; DepMap MSI threshold set to the DepMap definition (>= 20).

## 1.0.0 (2026-10-29)
- Initial pipeline: TCGA MMR grouping, TCGA/CPTAC/DepMap correlations, gnomAD/ClinVar, integration.
