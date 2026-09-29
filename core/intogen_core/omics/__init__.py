"""
Support for non-mutational omics layers: DNA methylation (array beta values)
and transcriptomics (RNA-seq).

The omics layers are analysed per cohort and integrated with the
mutation-based driver discovery at three levels:

1. Cohort-matched expression is used to flag non-expressed genes,
   replacing the TCGA-based proxy used when no expression is available.
2. Every driver is annotated with omics features (epigenetic silencing,
   expression outliers, expression of mutant vs wild-type samples).
3. Optionally (``--integrate_omics``), the per-gene p-values of the
   epigenetic silencing and expression outlier tests are added as
   additional voters in the combination step.
"""
