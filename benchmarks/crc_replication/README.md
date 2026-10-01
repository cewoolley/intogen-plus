# Replication of the co-regulator findings in an independent colorectal cancer cohort

Pre-registered test, in a new cohort, of two findings from TCGA microsatellite-stable (MSS) colorectal cancer
(`benchmarks/tcga_coadread`):

- **H1**: truncating mutations in nine transcriptional co-regulators (EP300, CREBBP, KMT2C, KMT2B, SIN3A, ASXL1,
  ASXL2, USP9X, TNRC6B) are under positive selection in MSS tumours.
- **H2**: MSS tumours carrying one are more often stage IV at diagnosis.

The hypotheses, definitions and analysis plan are in [PREREGISTRATION.md](PREREGISTRATION.md). They are frozen in code
in `spec.py` and `gene_sets.tsv`; `replicate.py` records the commit and the SHA-256 of both files with the results.
The pipeline was written for the 100,000 Genomes Project colorectal cancers, and works with any cohort given in the
input format below. `results_tcga/` holds its output on the discovery cohort (TCGA, dry run).

## Inputs

`mutations.tsv.gz`, somatic mutations (PASS only), one row per mutation:

| Column | |
|---|---|
| `SAMPLE` | tumour identifier, as in `samples.tsv` |
| `CHROM`, `POS`, `REF`, `ALT` | GRCh38 or GRCh37; indels VCF-style or MAF-style |
| `T_ALT`, `T_DEPTH` | tumour alternate and total read counts (optional; allele fractions) |

It can be written from an IntOGen-style cohort table (`prepare_intogen.py`) or from somatic VCFs (`prepare_vcf.py`,
Strelka2 tier-1 counts or `AD`):

```bash
python prepare_intogen.py input <cohort table.tsv.gz> [more tables]
```

`prepare_intogen.py` reads `CHROMOSOME`, `POSITION`, `REF`, `ALT`, `SAMPLE` (and `DONOR` and read counts if present:
`T_ALT`/`T_DEPTH`, `t_alt_count` with `t_depth` or `t_ref_count`, or `VAF`), keeps one sample per donor as IntOGen
does, and writes indels MAF-style. Give it the table **before** IntOGen's `parse-variants` step: that step removes
hypermutated samples (WGS: more than 10,000 SNVs and above Q3 + 1.5 IQR of the cohort), which drops the MSI and POLE
tumours and can drop high-burden MSS tumours. The script warns when the table looks filtered. Without read counts,
allele fractions are unavailable: give the cohort's `PURITY` in `samples.tsv`, otherwise purity is left out of the H2
model (`spec.PURITY_MIN_AVAILABLE`) and the clonality analysis is skipped.

`samples.tsv`, one row per tumour (one tumour per patient):

| Column | |
|---|---|
| `SAMPLE` | tumour identifier |
| `MSI_STATUS` | `MSS`, `MSI` or `POLE` (MSI-H, MSI-L understood). If missing for any tumour, groups are derived from the mutations as in TCGA |
| `SITE_ICD10` or `LOCATION` | ICD-10 site (C18.0 to C20) or `proximal` / `distal` / `rectum` |
| `STAGE` | stage at diagnosis, I to IV (e.g. `Stage IIIB`, `3`, `IV`: only the Roman numeral is read, so map Arabic stage codes to I-IV) |
| `AGE`, `SEX` | at diagnosis |
| `PURITY` | tumour purity (optional; otherwise the median allele fraction of the tumour is used) |
| `OS`, `OS_TIME`, `PFI`, `PFI_TIME` | event (0/1) and time in days (optional; secondary analyses) |

## Run

```bash
# 1. blinded feasibility (cohort size, neutral expectation and stage mix; record it with the registration)
FEASIBILITY=1 ./run.sh input hg38 results
# 2. the registered analysis
./run.sh input hg38 results [cornish_drivers.tsv]
```

`run.sh` keeps the coding and splice-site mutations, assigns the groups, runs dNdScv (`dnds.R`) in the MSS and MSI
tumours, then `replicate.py`. The feasibility step reads no observed co-regulator mutation and no carrier-by-stage
table, so it can be run and reported before the registered analysis.

Outputs: `replication.json` and `summary.md` (or `feasibility.json`), aggregate results only: no sample identifiers,
and counts of tumours below 5 (`spec.MIN_CELL`) are suppressed. Check them against the current Genomics England
airlock rules before export. `work/` holds individual-level intermediate files and stays in the Research Environment.

## In the Genomics England Research Environment

- Import this repository through the airlock (the pipeline uses `core/intogen_core` of this fork).
- Software: R with `dndscv` (its package data include the GRCh38 RefCDS and covariates) and `data.table`;
  Python 3 with numpy, pandas, scipy and statsmodels.
- Tumours: the colorectal primary tumours of the cancer programme, one per participant (the one analysed by
  Cornish et al., Nature 2024, where it applies). MSI and POLE status: the cohort's calls if available.
- Somatic small variants: an IntOGen-style cohort table prepared in the Research Environment (`prepare_intogen.py`,
  before `parse-variants`), or the per-tumour somatic VCFs of the cancer analysis (`prepare_vcf.py`, manifest of
  SAMPLE and VCF path; run `prepare_vcf.py --check 20` on one VCF first: Strelka2 sample columns and FORMAT fields are
  assumed).
- Clinical data: stage, site, age and sex at diagnosis from the cancer registry tables; death from the mortality
  data. Table and column names differ between data releases, so they are mapped by hand into `samples.tsv`.
- Cornish et al.'s driver list (supplementary table) as a TSV/CSV with a `gene` column enables the gene-level
  comparison (secondary analysis S8).

## Files

| File | |
|---|---|
| `PREREGISTRATION.md` | hypotheses and analysis plan |
| `spec.py`, `gene_sets.tsv` | frozen definitions and the nine Reactome sets (MSigDB 7.5.1) |
| `run.sh` | the pipeline |
| `prepare_intogen.py`, `prepare_vcf.py`, `prepare_tcga.py` | inputs from an IntOGen cohort table, somatic VCFs, or the TCGA benchmark |
| `coding_filter.py`, `classify.py`, `dnds.R` | coding mutations, groups and location, dNdScv |
| `replicate.py` | primary, secondary and feasibility analyses |
| `power.py` | power of H1 and H2 (`results/power.tsv`) |
