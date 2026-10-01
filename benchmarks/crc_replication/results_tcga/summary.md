# Replication of the co-regulator findings

Cohort: 553 tumours ({'MSS': 469, 'MSI': 73, 'POLE': 11}), groups from cohort.

## Primary hypotheses (fixed sequence H1 then H2, one-sided, family-wise alpha 0.05)

| Hypothesis | Estimate | 95% CI | p (one-sided) | Threshold | Supported |
|---|---|---|---|---|---|
| H1 truncating excess in the nine co-regulators (MSS) | 38 vs 5.1 expected, ratio 7.45 | 4.78-11.75 | 1.68e-16 | 0.050 | yes |
| H2 carriers more often stage IV (MSS, adjusted) | OR 2.77 (n = 451) | 1.35-5.72 | 0.00285 | 0.050 | yes |

Stage IV: carriers 13 of 44 (30%), others 59 of 407 (14%); purity covariate: median allele fraction.

## Secondary

| Reactome set | Observed / expected | q | Without the nine genes | p |
|---|---|---|---|---|
| REACTOME_RUNX1_REGULATES_GENES_INVOLVED_IN_MEGAKARYOCYTE_DIFFERENTIATION_AND_PLATELET_FUNCTION | 30 / 9.9 | 9.44e-06 | 8 / 7.4 | 0.464 |
| REACTOME_TRANSCRIPTIONAL_REGULATION_BY_RUNX1 | 61 / 31.1 | 1.5e-05 | 36 / 27.6 | 0.0778 |
| REACTOME_STAT3_NUCLEAR_EVENTS_DOWNSTREAM_OF_ALK_SIGNALING | 12 / 2.0 | 1.62e-05 | 3 / 1.3 | 0.141 |
| REACTOME_ACTIVATION_OF_ANTERIOR_HOX_GENES_IN_HINDBRAIN_DEVELOPMENT_DURING_EARLY_EMBRYOGENESIS | 27 / 11.1 | 0.000165 | 14 / 8.8 | 0.0709 |
| REACTOME_FOXO_MEDIATED_TRANSCRIPTION | 23 / 8.9 | 0.000165 | 11 / 7.2 | 0.12 |
| REACTOME_SUMOYLATION_OF_TRANSCRIPTION_COFACTORS | 21 / 7.8 | 0.000165 | 9 / 6.1 | 0.165 |
| REACTOME_NEF_AND_SIGNAL_TRANSDUCTION | 10 / 2.0 | 0.000167 | 10 / 2.0 | 0.000133 |
| REACTOME_TRAF6_MEDIATED_IRF7_ACTIVATION | 14 / 4.0 | 0.000167 | 6 / 2.7 | 0.0611 |
| REACTOME_DEUBIQUITINATION | 65 / 39.8 | 0.000233 | 47 / 37.7 | 0.0856 |

Matched random genes: stage IV fraction among carriers 0.30, null median 0.14, 95th percentile 0.22, empirical p = 0.0012.

| Gene | Truncating obs / exp (substitutions) | Carriers | Stage IV | OR without it | p without it |
|---|---|---|---|---|---|
| EP300 | 5 / 0.39 | <5 | 0 | 3.20 | 0.00097 |
| CREBBP | 3 / 0.93 | <5 | <5 | 2.70 | 0.0045 |
| KMT2C | 5 / 1.02 | 8 | <5 | 2.78 | 0.00499 |
| KMT2B | 3 / 0.33 | <5 | <5 | 2.59 | 0.00595 |
| SIN3A | 4 / 0.38 | 5 | <5 | 2.83 | 0.00337 |
| ASXL1 | 5 / 0.52 | 5 | <5 | 2.24 | 0.0228 |
| ASXL2 | 3 / 0.20 | <5 | <5 | 2.82 | 0.00341 |
| USP9X | 5 / 0.94 | 7 | <5 | 2.30 | 0.0195 |
| TNRC6B | 5 / 0.40 | 7 | <5 | 3.41 | 0.000573 |

Sensitivity of H2: unadjusted OR 2.47 (p = 0.00584); without location OR 2.62 (p = 0.0041); substitution carriers only OR 3.37 (p = 0.000922).

Relative allele fraction (median, clonal fraction): core truncating 1.27, 83% (n = 48); apc tp53 truncating 1.24, 87% (n = 623); synonymous 1.00, 69% (n = 16059).

OS: unadjusted HR 1.59 (0.87-2.93), p = 0.135; stage stratified HR 1.33 (0.71-2.48), p = 0.37; stages i iii HR 1.12 (0.44-2.84), p = 0.806.

PFI: unadjusted HR 1.25 (0.70-2.23), p = 0.449; stage stratified HR 1.05 (0.58-1.88), p = 0.882; stages i iii HR 0.67 (0.24-1.83), p = 0.433.

MSI tumours (73): co-regulator truncating 18 vs 10.3 (p = 0.0238).

Commit 3cc0ca1a5ccccf15146b9550fb78f43b1c3cc970, spec.py sha256 86b25d40ff3c, gene_sets.tsv sha256 dc65b83ab1a6, dNdScv 0.0.1.0, run 2026-10-01T21:57:00+00:00.
