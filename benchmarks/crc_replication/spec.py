"""
Frozen definitions of the pre-registered replication (PREREGISTRATION.md). Changing anything here after the
registration changes the registered analysis: record it as a deviation.
"""

# the nine transcriptional co-regulators found in TCGA MSS colorectal cancer (benchmarks/tcga_coadread)
CORE = ['EP300', 'CREBBP', 'KMT2C', 'KMT2B', 'SIN3A', 'ASXL1', 'ASXL2', 'USP9X', 'TNRC6B']

# the nine Reactome sets selected through truncating mutations in TCGA MSS tumours (members in gene_sets.tsv)
SETS = [
    'REACTOME_RUNX1_REGULATES_GENES_INVOLVED_IN_MEGAKARYOCYTE_DIFFERENTIATION_AND_PLATELET_FUNCTION',
    'REACTOME_TRANSCRIPTIONAL_REGULATION_BY_RUNX1',
    'REACTOME_STAT3_NUCLEAR_EVENTS_DOWNSTREAM_OF_ALK_SIGNALING',
    'REACTOME_ACTIVATION_OF_ANTERIOR_HOX_GENES_IN_HINDBRAIN_DEVELOPMENT_DURING_EARLY_EMBRYOGENESIS',
    'REACTOME_FOXO_MEDIATED_TRANSCRIPTION',
    'REACTOME_SUMOYLATION_OF_TRANSCRIPTION_COFACTORS',
    'REACTOME_NEF_AND_SIGNAL_TRANSDUCTION',
    'REACTOME_TRAF6_MEDIATED_IRF7_ACTIVATION',
    'REACTOME_DEUBIQUITINATION',
]

# dNdScv impacts counted as truncating. H1 counts substitutions only (as the TCGA test); a carrier (H2) also
# counts frameshift indels
TRUNCATING_SUBSTITUTIONS = ['Nonsense', 'Essential_Splice']

# drivers of a stratum: dNdScv qglobal_cv below this threshold, excluded from the background fit and the long tails
DRIVER_Q = 0.1

# groups derived from mutations when the cohort has no MSI / POLE calls
POLE_DOMAIN = (268, 471)               # exonuclease domain codons
HYPERMUTATED_CODING = 500              # coding mutations (TCGA: gap of the bimodal burden distribution)

# location from ICD-10 site codes (C18.1 appendix and unspecified sites are left out)
PROXIMAL = ['C18.0', 'C18.2', 'C18.3', 'C18.4']
DISTAL = ['C18.5', 'C18.6', 'C18.7']
RECTUM = ['C19', 'C19.9', 'C20', 'C20.9']

# H2 model: stage IV ~ carrier + covariates (logistic regression, complete cases)
H2_COVARIATES = ['age_10y', 'male', 'rectum', 'proximal', 'log_coding', 'purity']
# purity: the cohort's estimates when given for at least 90% of tumours, else the median allele fraction; when
# neither is available for at least half of the tumours (e.g. no read counts), purity is left out of the model
PURITY_MIN_COHORT = 0.9
PURITY_MIN_AVAILABLE = 0.5

# multiplicity: fixed sequence, H1 then H2, each one-sided at ALPHA; H2 is confirmatory only if H1 is supported
# (family-wise error rate ALPHA)
ALPHA = 0.05

# matched random-gene null of H2 (secondary)
NULL_DRAWS = 5000
SEED = 11

# aggregate outputs: counts of tumours below this are suppressed (check the current airlock rules)
MIN_CELL = 5
