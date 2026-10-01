# Pre-registration: truncating mutations in transcriptional co-regulators in microsatellite-stable colorectal cancer

Replication in the 100,000 Genomes Project colorectal cancers of two findings from TCGA.

Format: OSF Preregistration, registration of existing data prior to analysis. The analysis code is part of this
registration: `benchmarks/crc_replication` of github.com/cewoolley/intogen-plus at the commit recorded on
registration, with SHA-256 of `spec.py` and `gene_sets.tsv` (printed by `replicate.py` in every output). At the
time of writing: `spec.py` 86b25d40ff3cb1eb8d3244a9bc800911bda7345a6c55d4b464fc5de77459d479, `gene_sets.tsv` dc65b83ab1a65785a41ac7efe2e1ae7e2b2a9f7f345470899b054fed340deaa7 (commit 3cc0ca1).

Items marked **[confirm]** must be checked and completed by the investigators before registering.

## 1. Study information

**Title.** Selection of truncating mutations in enhancer co-regulators and their association with stage IV disease
in microsatellite-stable colorectal cancer: a pre-registered replication in the 100,000 Genomes Project.

**Authors.** **[confirm]**

**Background.** In 469 microsatellite-stable (MSS) colorectal cancers from TCGA (MC3 exome calls), a pathway-level
analysis of genes too rarely mutated to be called drivers individually found nine Reactome gene sets selected through
truncating mutations (q = 0.005 to 0.096). They are one signal carried by nine transcriptional co-regulators: EP300,
CREBBP, KMT2C, KMT2B, SIN3A, ASXL1, ASXL2, USP9X and TNRC6B. In TCGA:

- 38 truncating substitutions in the nine genes against 5.1 expected under the neutral model. The genes were chosen
  after inspecting the data, so this ratio is inflated and is not itself the evidence.
- Missense mutations near neutral (67 vs 51.4) and synonymous at expectation (19 vs 18.5); allele fractions as clonal
  as APC/TP53 truncations; no excess of filtered calls or of the oxidative-damage substitution.
- 44 of 469 MSS tumours (9%) carry a truncating mutation (nonsense, essential splice or frameshift) in one of the
  genes. 13 of 44 carriers were stage IV at diagnosis against 59 of 407 other MSS tumours (14%): adjusted odds ratio
  2.8 (95% CI 1.4-5.8). Random gene sets matched for their number of carriers: 14% stage IV (empirical p = 0.0006).
  No difference in progression-free interval within stage (hazard ratio 1.05, 13 events among carriers).

The stage IV association was found among several clinical features examined after the gene list was chosen, so it
is a hypothesis.

**Hypotheses.**

- **H1.** In MSS colorectal cancers, truncating substitutions (nonsense and essential splice) in the nine genes
  exceed the neutral expectation of the fork's null model (one-sided).
- **H2.** MSS colorectal cancers carrying a truncating mutation (nonsense, essential splice or frameshift) in any of
  the nine genes are more often stage IV at diagnosis, adjusting for age, sex, rectum, proximal colon, coding
  mutation burden and purity (odds ratio > 1, one-sided).

## 2. Design

Observational analysis of existing whole-genome sequencing and registry data. No intervention, no randomisation.

**Blinding.** Before registration, the investigators may run the blinded feasibility step only
(`FEASIBILITY=1 ./run.sh`): the number of MSS tumours, the number with stage, the overall proportion at stage IV,
and the neutral expectation of H1. It reads no observed mutation in the nine genes and no carrier-by-stage table. Its
output is attached to the registration.

## 3. Sampling plan

**Existing data.** Registration prior to analysis of the data. **[confirm]** The investigators have not analysed
mutations in the nine genes, or their association with stage, in the 100,000 Genomes Project data. Prior knowledge:
Cornish et al. (Nature 2024) report driver genes in these tumours; whether any of the nine genes are among them is
**[confirm: state what the investigators know]**. H1 is a pathway-level test of the nine genes together and is not
equivalent to their individual driver calls; H2 has not been reported for this cohort to our knowledge.

**Data.** Colorectal adenocarcinomas of the 100,000 Genomes Project cancer programme with tumour and matched normal
whole-genome sequencing (the cohort of Cornish et al., about 2,000 tumours). One primary tumour per participant (the
one analysed by Cornish et al. where it applies, else the first in alphabetical order of sample identifier, as
IntOGen does). Somatic small variants: the programme's somatic calls (PASS), from the investigators' IntOGen-style
cohort table taken before IntOGen's hypermutator filter (`parse-variants`), so that no tumour is removed for its
mutation burden. Stage at diagnosis and site: cancer registry. Age and sex at diagnosis. Death: mortality data.
Purity: the cohort's estimates if available **[confirm]**.

**Groups.** MSI and POLE status from the cohort's existing calls (Cornish et al.) when available for every tumour;
otherwise from the mutations as in TCGA (`spec.py`: POLE = exonuclease-domain missense in a hypermutated tumour,
hypermutated = more than 500 coding mutations).

**Sample size.** Fixed by the available data; expected about 1,400 to 1,700 MSS tumours with stage.

**Power** (`power.py`, `results/power.tsv`; one-sided α = 0.05):

| | Effect | 1,000 MSS | 1,400 MSS | 1,700 MSS |
|---|---|---|---|---|
| H1 | observed / expected = 2 | 0.82 | 0.93 | 0.96 |
| H1 | observed / expected = 3 | 1.00 | 1.00 | 1.00 |
| H2 (9% carriers, 8% stage IV) | odds ratio 2.0 | 0.70 | 0.82 | 0.88 |
| H2 (9% carriers, 8% stage IV) | odds ratio 2.8 | 0.97 | 0.99 | 1.00 |
| H2 (6% carriers, 5% stage IV) | odds ratio 2.0 | 0.44 | 0.55 | 0.62 |

The TCGA estimates are expected to shrink in a new cohort; the planning effect for H2 is an odds ratio of 2.0. Power
for the observed cohort size and stage mix is computed by the feasibility step and reported with the registration.
The analysis is run whatever the power.

## 4. Variables

All definitions are in `spec.py` and the code; the text below summarises them.

- **Truncating substitution (H1):** dNdScv impact `Nonsense` or `Essential_Splice`.
- **Carrier (H2):** at least one truncating substitution or frameshift indel (dNdScv `no-SNV` with a length change not
  divisible by 3) in any of the nine genes, among the mutations dNdScv annotates in the MSS tumours.
- **Stage IV:** stage group IV at diagnosis (registry). Tumours without stage are excluded from H2.
- **Covariates (H2):** age / 10, male, rectum (C19-C20), proximal colon (C18.0, C18.2-C18.4), log coding mutations,
  purity. Purity is the cohort's estimate when available for at least 90% of tumours, else the median allele fraction
  of the tumour; if neither is available for at least half of the tumours (e.g. a mutation table without read
  counts), purity is left out of the model. Complete cases.
- **Drivers of a stratum:** dNdScv `qglobal_cv` < 0.1 in that stratum; excluded from the background fit of the null
  model and from the long tails of the secondary set tests (but the nine genes are always tested in H1).

## 5. Analysis plan

**Pipeline.** `run.sh`: coding and splice-site mutations (`dnds.R regions`), groups (`classify.py`), dNdScv per
stratum on GRCh38 (`dnds.R strata`, package defaults), analyses (`replicate.py`).

**H1.** The fork's null model of truncating substitutions in MSS tumours (`MutationLayer('mutation_truncating')` after
`background_omega`, which corrects the neutral expectation for genomic context with a regression fitted on non-driver
genes and the 20 epigenomic covariate components of dNdScv). Statistic: total truncating substitutions in the nine
genes. p: upper tail of the sum of the per-gene negative binomial posterior predictive distributions. Estimate:
observed / expected with an exact 95% interval.

**H2.** Logistic regression of stage IV on carrier status and the covariates. Estimate: odds ratio with Wald 95% CI.
p: one-sided Wald test of a positive coefficient.

**Multiplicity.** Fixed sequence: H1 at one-sided α = 0.05; H2 at one-sided α = 0.05 only if H1 is supported. This
keeps the family-wise error rate at 0.05. If H1 is not supported, H2 is reported as a secondary result.

**Secondary analyses** (estimates with nominal p; no confirmatory claim):

- S1. The nine Reactome sets: long-tail truncating test, with and without the nine genes (BH q over the nine sets).
- S2. Stage IV among carriers against 5,000 sets of random non-driver genes matched for their number of carriers.
- S3. Per gene: carriers, stage IV carriers, and H2 leaving each gene out.
- S4. Clonality: allele fraction relative to the tumour median for co-regulator truncations, APC/TP53 truncations and
  synonymous mutations.
- S5. Overall survival (and progression-free interval if available): Cox models of carriers, unadjusted, stratified
  by stage, and in stages I-III.
- S6. H1 in MSI tumours.
- S7. Sensitivity of H2: unadjusted; without the location covariates; carriers of substitutions only.
- S8. Gene-level overlap of the nine genes with the drivers of Cornish et al.

Any other analysis is exploratory and labelled as such.

**Inference criteria.** H1 is replicated if its one-sided p ≤ 0.05. H2 is replicated if H1 is replicated and its
one-sided p ≤ 0.05. Effect sizes and intervals are reported in all cases, including non-replication.

**Data exclusion.** Tumours excluded by dNdScv (more than 3,000 coding mutations) are left out of H2. Tumours with
missing stage are left out of H2; missing covariates, of the adjusted model (their number is reported).

**Missing data.** Complete-case analysis. If more than 20% of MSS tumours with stage miss a covariate, the model
without that covariate is reported as well (exploratory).

## 6. Other

**Deviations.** Any change to the code, data or definitions after registration is reported with its reason and its
effect on the results.

**Data governance.** Analyses run in the Genomics England Research Environment. Only aggregate results leave it:
no identifiers, and counts of tumours below 5 are suppressed.

**Reporting.** The results are reported whatever their direction.

**Discovery cohort through the same pipeline.** `results_tcga/` holds the output of this pipeline on TCGA (GRCh37,
MC3 exomes, 469 MSS tumours), as a check that the code reproduces the discovery estimates: H1 38 vs 5.1 expected
(ratio 7.45, 95% CI 4.78-11.75); H2 13 of 44 carriers stage IV against 59 of 407 others, adjusted odds ratio 2.77
(1.35-5.72). These are the estimates the replication is compared with; both are expected to shrink.
