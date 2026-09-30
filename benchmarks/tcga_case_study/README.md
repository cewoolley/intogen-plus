# Case study: pathway analysis on TCGA exomes

Evaluation of the pathway analysis of this fork (`pathway-analysis`) with real
somatic mutations, against what the original IntOGen provides (gene-level
drivers) and the usual follow-up analyses:

- gene-set selection: over-representation of the drivers (ORA) and dNdScv
  `genesetdnds` (Martincorena et al. 2017), with random gene sets
  (calibration) and selection spiked into the long tail of gene sets (power);
- co-occurrence: pairwise Fisher's exact tests on the most mutated genes (as in
  maftools `somaticInteractions` or cBioPortal) and the additive model of
  DISCOVER (Canisius et al. 2016), on the real cohorts and with co-occurrence
  planted in real tumours.

The results are summarised in `docs/source/pathways.rst` (section *Evaluation
on TCGA exomes*) and in `results/summary.json`.

## Data

- Somatic mutations: TCGA MC3 (KIRC, LUAD, BRCA, COAD) as distributed by
  [TCGAmutations](https://github.com/PoisonAlien/TCGAmutations) (`inst/extdata/MC3`).
- Gene sets: MSigDB 7.5.1 hallmarks and Reactome, as bundled in
  [msigdbr 7.5.1](https://github.com/cran/msigdbr/tree/7.5.1) (`R/sysdata.rda`).
- dNdScv ([im3sanger/dndscv](https://github.com/im3sanger/dndscv), commit 43c5e2f),
  with its default hg19 reference and epigenomic covariates.

The original IntOGen is represented by its gene-level output: the genes with
dNdScv q-value < 0.1 are the drivers (and the genes excluded from the long tail).
A sensitivity analysis also excludes all the genes of the Cancer Gene Census
(v81, bundled with dNdScv).

## Requirements

R with dndscv, data.table (and the dependencies of dndscv), and the Python
environment of `intogen-core` (the scripts import it from this repository).

## Run

```bash
export TCGAMUTATIONS=/path/to/TCGAmutations MSIGDBR_SYSDATA=/path/to/msigdbr/R/sysdata.rda
./run.sh      # a few hours, most of it genesetdnds (about 35 CPU minutes per cohort)
```

| Step | Script | Output |
|---|---|---|
| dNdScv, gene sets, covariates | `prepare.R` | `data/` |
| inputs of the fork, tests against the neutral model, units for genesetdnds | `inputs.py` | `results/<COHORT>/fork.*` |
| dNdScv genesetdnds | `gsd.R` | `results/<COHORT>/gsd.tsv` |
| pathway analysis of the fork (default: background dN/dS) | `corrected.py` | `results/<COHORT>/final.*` |
| calibration with random gene sets | `calibration.py` | `results/<COHORT>/calibration_final.tsv` |
| power with selection spiked into the long tail | `spikein.py` | `results/<COHORT>/spikein_final.tsv` |
| gene sets called by each approach | `compare_sets.py` | `results/sets_compare.json` |
| co-occurrence: Fisher, DISCOVER, burden elasticity | `cooccurrence.py` | `results/cooccurrence.json` |
| summary | `consolidate.py` | `results/summary.json` |
