# Driver networks in TCGA colorectal cancer (COADREAD)

Mining of TCGA COADREAD (MC3 somatic mutations, 553 tumours) with the pathway analysis of this fork:
gene-level drivers (dNdScv), selection of the long tail of gene sets, and co-occurrence / mutual exclusivity
of drivers and selected long tails with the burden-aware test, in microsatellite-stable (MSS) and
hypermutated (MSI) tumours separately. The results are summarised in `results/crc_summary.json` and, for the
transcriptional co-regulator signal of MSS tumours, `results/coregulators.json`.

Methylation, RNA-seq and copy number of these tumours could not be obtained in the environment where this
was run, so only mutations are analysed.

## Data

- Somatic mutations and clinical data: TCGA MC3 COAD and READ as distributed by
  [TCGAmutations](https://github.com/PoisonAlien/TCGAmutations) (`inst/extdata/MC3`).
- Gene sets: MSigDB 7.5.1 hallmarks and Reactome (msigdbr 7.5.1) and seven curated CRC pathways
  (`gene_sets.py`: WNT, TGF-β/BMP, PI3K, RTK/RAS, TP53/DNA damage, SWI/SNF and immune escape).
- dNdScv ([im3sanger/dndscv](https://github.com/im3sanger/dndscv), commit 43c5e2f) with its default hg19
  reference and epigenomic covariates.

## Groups

Tumours are grouped from their mutations (`classify.py`): POLE (missense mutation in the exonuclease domain,
codons 268-471; 11 tumours), MSI (more than 500 coding mutations, the gap of the bimodal burden distribution;
73 tumours, median indel fraction 0.15, 60 of them proximal) and MSS (469). The clinical field
`loss_expression_of_mismatch_repair_proteins_by_ihc` is not consistent with the mutation data and is not used.

## Run

```bash
export TCGAMUTATIONS=/path/to/TCGAmutations MSIGDBR_SYSDATA=/path/to/msigdbr/R/sysdata.rda
./run.sh      # about 15 minutes
```

| Step | Script | Output |
|---|---|---|
| mutations, read counts, clinical data and TCGA CDR survival | `mutations.R` | `coadread_mutations.tsv.gz`, `mutations_reads.tsv.gz`, `clinical.tsv` |
| MSS / MSI / POLE groups and location | `classify.py` | `tumour_groups.tsv` |
| dNdScv per stratum (ALL, MSS, MSI), gene sets, covariates | `dnds.R` | `data/` |
| curated CRC pathways | `gene_sets.py` | `gene_sets_crc.tsv.gz` |
| pathway analysis of the fork per stratum | `crc_pathways.py` | `results/<STRATUM>/` |
| exhaustive pairs, location and MSI associations, coherence | `crc_networks.py` | `results/<STRATUM>/network_pairs.tsv.gz` |
| location-stratified tests, composition of the co-regulator signal | `crc_robust.py` | `results/robustness.json` |
| the co-regulator long tail of MSS tumours: genes, overlap of the selected sets, MC3 filters, allele fractions, tumours | `crc_coregulators.py` | `results/coregulators.json` |
| stage IV and survival of co-regulator carriers: adjusted model, matched random genes, TCGA CDR endpoints | `crc_coregulators_clinical.py` | `results/coregulators_clinical.json` |
| immune escape in hypermutated tumours | `crc_immune.py` | `results/immune_escape.json` |
| detectable effect sizes (469 vs 2,023 tumours) | `crc_power.py` | `results/power.tsv` |
| summary | `crc_summary.py` | `results/crc_summary.json` |
