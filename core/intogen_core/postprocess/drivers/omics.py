"""
Integration of the cohort omics features (omics-features output)
into the driver post-processing.
"""

import numpy as np
import pandas as pd


# Mutation-based driver identification methods
MUTATION_METHODS = {'oncodriveclustl', 'dndscv', 'oncodrivefml', 'hotmaps', 'smregions', 'cbase', 'mutpanning'}

# Columns added to the vetting (unfiltered drivers) output
OMICS_VET_COLUMNS = ['WARNING_EXPRESSION_SOURCE', 'QVALUE_METHYLATION', 'QVALUE_EXPRESSION']

# Columns added to the drivers output
OMICS_DRIVERS_COLUMNS = [
    'EXPRESSED', 'EXPRESSION_OUTLIER_DIRECTION', 'EXPRESSION_OUTLIER_SAMPLES', 'QVALUE_EXPRESSION',
    'EXPRESSION_LOG2FC_MUTANT', 'QVALUE_EXPRESSION_MUTANT',
    'EXPRESSION_LOG2FC_TRUNCATING', 'PVALUE_EXPRESSION_TRUNCATING',
    'METHYLATION_HYPER_SAMPLES', 'METHYLATION_FUNCTIONAL', 'QVALUE_METHYLATION',
    'SAMPLES_MUTATED_AND_HYPERMETHYLATED', 'SAMPLES_MUTATED_OR_HYPERMETHYLATED',
    'OMICS_ROLE_SUPPORT'
]


# Boolean columns (True/False/NA), read as text so that concatenating cohorts
# with and without them does not turn them into numbers
OMICS_BOOLEAN_COLUMNS = {'EXPRESSED': str, 'METHYLATION_FUNCTIONAL': str}


def load_omics(path):
    """Load the omics features of the cohort (None if not available)"""
    if path is None:
        return None
    df = pd.read_csv(path, sep='\t', dtype=OMICS_BOOLEAN_COLUMNS)
    for c in ['SYMBOL'] + OMICS_DRIVERS_COLUMNS:
        if c not in df.columns:
            df[c] = np.nan
    return df.drop_duplicates(subset=['SYMBOL'])


def _as_bool(value):
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, str) and value.lower() in ('true', 'false'):
        return value.lower() == 'true'
    return None


def apply_cohort_expression(df, omics):
    """
    Use the expression of the cohort to flag non-expressed genes.

    Genes measured in the cohort are flagged according to it; the rest keep
    the flag computed from TCGA.

    Args:
        df: dataframe with the GENE and Warning_Expression (TCGA-based) columns
        omics: omics features of the cohort

    Returns:
        df with the updated Warning_Expression and a Warning_Expression_Source column
    """
    expressed = {g: _as_bool(v) for g, v in zip(omics['SYMBOL'], omics['EXPRESSED'])}
    cohort = df['GENE'].map(lambda g: expressed.get(g))
    measured = cohort.notna()
    df = df.copy()
    df['Warning_Expression'] = np.where(measured, cohort.map(lambda v: v is False), df['Warning_Expression'])
    df['Warning_Expression'] = df['Warning_Expression'].astype(bool)
    df['Warning_Expression_Source'] = np.where(measured, 'cohort', 'TCGA')
    return df


def has_mutational_bidder(significant_bidders):
    """Whether any mutation-based method is among the significant bidders"""
    return any(m in MUTATION_METHODS for m in str(significant_bidders).split(','))
