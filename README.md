# mmr-gene-explorer

Test whether a gene of interest is linked to **mismatch repair (MMR) deficiency** in human cancer,
using only public data: TCGA, CPTAC, DepMap, gnomAD and ClinVar.

Given a target gene (for example `DRG2`) and partner genes (for example `MLH1`, `MLH3`), the pipeline:

1. Classifies TCGA tumours into MMR groups: MSS, MSI with MLH1 promoter methylation (sporadic),
   MSI with a truncating MMR mutation, MSI unexplained (a tumour-level proxy for Lynch-like syndrome)
   and POLE-ultramutated.
2. Tests whether the target is lower, deleted or mutated in MSI-unexplained tumours, and whether it
   correlates with the partner genes, with adjustment for arm-level copy-number loss.
3. Repeats the correlations at the protein level (CPTAC) and across cancer cell lines (DepMap),
   including a genome-wide rank and MSI status.
4. Summarises population constraint and loss-of-function carriers (gnomAD) and clinical variant
   records (ClinVar).
5. Labels every result **Supportive / Contradicts / Inconclusive** using rules fixed in advance, and
   writes a Markdown summary.

Defaults cover endometrial and colorectal cancer, the main Lynch syndrome tumour types. Any TCGA
PanCancer study on cBioPortal can be used instead.

## Quick start

```bash
git clone https://github.com/<ManiMuralidharan>/mmr-gene-explorer.git
cd mmr-gene-explorer
```

In R (4.3 or later), from the repository folder:

```r
install.packages(c("data.table", "ggplot2", "jsonlite", "httr"))
source("R/00_download_data.R")   # TCGA + CPTAC (about 1 GB), gnomAD, ClinVar
```

Download DepMap files manually from <https://depmap.org/portal/data_page/?tab=allData> into
`data/raw/depmap/`:

| File | Required |
|---|---|
| `Model.csv` | yes |
| `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` (26Q3+; older: `OmicsExpressionProteinCodingGenesTPMLogp1.csv`) | yes |
| `OmicsGlobalSignatures.csv` (MSI score; older: `OmicsSignatures.csv`) | recommended |
| `CRISPRGeneEffect.csv` | optional |

Then run everything:

```bash
Rscript R/run_all.R
```

Results go to `results/<TARGET_GENE>/` and the summary to `reports/<TARGET_GENE>_database_summary.md`.

## Analysing a different gene or cancer

Edit the **USER SETTINGS** block at the top of `R/00_config.R`, or use environment variables:

```bash
TARGET_GENE=ATM PARTNER_GENES=MLH1,MSH2 CNA_PROXY_GENE=ATM Rscript R/run_all.R
```

| Setting | Meaning | Default |
|---|---|---|
| `GENE_TARGET` / `TARGET_GENE` | Gene being tested | `DRG2` |
| `GENES_PARTNER` / `PARTNER_GENES` | Genes expected to change with the target | `MLH1, MLH3` |
| `CNA_PROXY_GENE` | Gene on the same chromosome arm, used to adjust for arm-level loss (`NA` to skip) | `TP53` |
| `TCGA_STUDIES` | cBioPortal PanCancer study IDs | UCEC, COADREAD |
| `CPTAC_STUDIES` | cBioPortal CPTAC study IDs | UCEC 2020, COAD 2019 |
| `DEPMAP_SUBGROUPS` | DepMap lineage subgroups (regular expressions) | Colorectal, Endometrial |

Study IDs are listed at <https://www.cbioportal.org/datasets>. Downloaded data are shared between runs;
only results are separated per target gene, so several genes can be compared side by side.

Before running a new gene, add its hypothesis to `metadata/cohort_definitions.md` so the
interpretation is fixed in advance.

### Other questions the pipeline can answer with small changes

- **A different DNA repair pathway:** change `GENES_PARTNER` (for example homologous recombination
  genes) and the grouping logic in script 01 (for example HRD scores instead of MSI).
- **Another cancer type:** add its PanCancer study ID to `TCGA_STUDIES`. MSI grouping is most
  meaningful in endometrial, colorectal and stomach cancer.
- **Expression-only screening of many genes:** loop `run_all.R` over a gene list with `TARGET_GENE`.

## Pipeline

| Script | Purpose |
|---|---|
| `R/00_config.R` | User settings, paths, thresholds |
| `R/utils.R` | Shared readers and helpers |
| `R/00_download_data.R` | Downloads open data and logs sources in `metadata/data_sources.csv` |
| `R/01_tcga_cohort_definition.R` | MMR groups from MSI score, MLH1 promoter methylation and MMR mutations |
| `R/02_tcga_target_expression_cnv_mutation.R` | Target expression by group, correlations, CNA and mutation enrichment |
| `R/03_cptac_protein_correlations.R` | Protein-level correlations and target protein by MSI status |
| `R/04_depmap_expression_msi.R` | Cell-line correlations, genome-wide rank, MSI, CRISPR co-dependency |
| `R/05_target_population_variants.R` | gnomAD constraint, LoF carriers, ClinVar |
| `R/06_integrate_results.R` | Evidence table with pre-specified verdicts and Markdown summary |
| `R/run_all.R` | Runs all steps and reports failures per step |

Group definitions and decision rules: [`metadata/cohort_definitions.md`](metadata/cohort_definitions.md).

## Interpretation cautions

- Correlations across tumours or cell lines cannot show that one gene regulates another.
- "MSI unexplained" is based on somatic data only. It approximates, but is not, a clinical diagnosis of
  Lynch-like syndrome. Germline analysis needs controlled-access data (dbGaP phs000178).
- Copy-number results can be confounded by arm-level loss and by the copy-number-quiet nature of MSI
  tumours.
- Low-abundance proteins (for example MLH3) are often missing from proteomics data.
- Public data formats change. If a file or column is renamed, scripts stop with a message naming what
  was not found.

## Data sources and terms

All data are downloaded from their providers and are not redistributed here. Please cite the original
sources and follow their terms of use: TCGA PanCancer Atlas and CPTAC (via cBioPortal), DepMap,
gnomAD and ClinVar.

## Reproducibility

```r
install.packages("renv")
renv::init()      # first time, records package versions in renv.lock
renv::restore()   # on another machine
```

Each script writes its `sessionInfo()` to `results/<TARGET_GENE>/sessionInfo/`.

## Citation

If you use this code, please cite it using the metadata in [`CITATION.cff`](CITATION.cff).

## License

MIT. See [`LICENSE`](LICENSE).
