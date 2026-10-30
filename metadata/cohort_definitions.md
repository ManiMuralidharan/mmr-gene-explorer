# Pre-specified cohort definitions and analysis rules

Version 1.1 (2026-10-30). Group rules were fixed before any target-gene result was examined.
Log every change at the bottom of this file with a date and reason.

## Hypothesis template

A target gene (`GENE_TARGET`, default DRG2) is required for normal expression of partner genes
(`GENES_PARTNER`, default MLH1 and MLH3). Loss of the target lowers mismatch repair (MMR) capacity and
contributes to MMR-deficient tumours, including Lynch-like syndrome (LLS): MMR-deficient tumours without a
germline MMR pathogenic variant and without somatic MLH1 promoter hypermethylation.

Expected direction if the hypothesis holds:
- Target mRNA/protein correlates positively with the partner genes.
- The target is lower, deleted or mutated more often in MSI_unexplained tumours than in MSS tumours and
  MLH1-methylated MSI tumours.

When you analyse a new target gene, add its hypothesis and the date to the change log before running
scripts 02 to 06.

## TCGA MMR groups (script 01)

Primary tumours only (sample type 01), PanCancer Atlas studies set in `TCGA_STUDIES`.

| Group | Rule |
|---|---|
| POLE_ultramutated | PanCancer SUBTYPE (patient file) contains POLE or HM-SNV; excluded from MSI groups |
| MSS | MANTIS score <= 0.4 (MSIsensor < 10 if MANTIS missing; SUBTYPE if both missing) |
| MSI_MLH1_methylated | MSI-H and MLH1 promoter score > 0.3 (sporadic MSI) |
| MSI_MMR_mutated | MSI-H, MLH1 score <= 0.3, somatic truncating mutation in MLH1, MSH2, MSH6, PMS2 or EPCAM |
| MSI_unexplained | MSI-H, MLH1 score <= 0.3, no truncating core MMR mutation, mutation-profiled (Lynch-like proxy) |
| MSI_indeterminate | MSI-H, MLH1 methylation not available; excluded from tests |

**MLH1 promoter score:** mean beta of all HM450 probes annotated to MLH1 with TSS200 in their
annotation (core promoter). HM27-only tumours use the TSS200 MLH1 probe(s) in the merged HM27/HM450 file
(cg00893636). Probe selection uses annotation only.

**Truncating:** Nonsense, Frame_Shift_Del/Ins, Splice_Site, Splice_Region, Translation_Start_Site, Nonstop.

Caveats:
- Open TCGA data are somatic. MSI_unexplained mixes true LLS, undetected germline Lynch syndrome, and
  tumours with biallelic somatic MMR inactivation by missense or copy-number events.
- Germline confirmation requires dbGaP phs000178 (controlled; `data/controlled/`, never shared).

## Primary tests

1. TCGA: target mRNA, MSI_unexplained vs each other group (Wilcoxon, BH-FDR); linear model adjusted
   for age (and sex for CRC).
2. TCGA: Spearman correlation of target vs partner genes, all tumours and MSS only.
3. CPTAC: Spearman correlation of target vs partner proteins.
4. DepMap: Spearman correlation of target vs partner genes (all lines; subgroups in `DEPMAP_SUBGROUPS`),
   and genome-wide rank of the partner genes among target-correlated genes.

## Secondary / descriptive

- Target deletion and mutation enrichment (Fisher). MSI tumours are copy-number quiet, so fewer
  deletions in MSI groups are expected regardless of target biology.
- Target to partner association adjusted for `CNA_PROXY_GENE` copy number (arm-level proxy) and MMR group.
- gnomAD constraint (pLI, LOEUF) and high-confidence LoF carrier frequency.
- ClinVar target-specific records; multi-gene CNVs reported separately.

## Decision rules (script 06)

- Supportive: FDR < 0.05 in the hypothesised direction.
- Contradicts: FDR < 0.05 in the opposite direction.
- Inconclusive: FDR >= 0.05.
- Nominal P values are labelled as such.

## Change log

- 2026-10-29: v1.0 created (target DRG2; partners MLH1, MLH3).
- 2026-10-30: v1.1 implementation fixes, group rules unchanged. MLH1 methylation read from probe-level
  files using TSS200 probes (upstream probes cg10990993/cg18320188 failed QC: methylated in most MSS
  tumours). PanCancer SUBTYPE read from the patient file. DepMap 26Q3 file names supported.
