"""
Gene-level omics features of a cohort (``omics-features``).

Combines the results of the omics analyses with the mutations of the cohort:

- expression status and expression outliers (``expression-analysis``)
- expression of mutated vs wild-type tumours: all non-synonymous mutations
  and truncating mutations only (the latter captures nonsense-mediated decay)
- promoter hypermethylation (``methylation-analysis``) and the number of
  tumours with a mutation and/or promoter hypermethylation of the same gene
  (alternative or complementary hits of a tumour suppressor)
- the mode of action supported by the omics data:

  - LoF: functional recurrent promoter hypermethylation or recurrent
    under-expression outliers
  - Act: recurrent over-expression outliers or higher expression in mutated
    tumours. Lower expression in mutated tumours is not used as LoF support
    because nonsense-mediated decay also affects truncating passenger mutations.
"""

import os

import click
import numpy as np
import pandas as pd

from intogen_core.omics import io
from intogen_core.omics.stats import fdr_bh, mannwhitney


SYNONYMOUS = {'synonymous_variant'}
TRUNCATING = {'stop_gained', 'frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant'}

OMICS_COLUMNS = [
    'SYMBOL',
    'EXPRESSED', 'FRACTION_EXPRESSED', 'MEDIAN_LOG2_EXPRESSION',
    'EXPRESSION_OUTLIER_DIRECTION', 'EXPRESSION_OUTLIER_SAMPLES', 'QVALUE_EXPRESSION',
    'EXPRESSION_MUTANT_SAMPLES', 'EXPRESSION_LOG2FC_MUTANT', 'PVALUE_EXPRESSION_MUTANT', 'QVALUE_EXPRESSION_MUTANT',
    'EXPRESSION_TRUNCATING_SAMPLES', 'EXPRESSION_LOG2FC_TRUNCATING', 'PVALUE_EXPRESSION_TRUNCATING',
    'METHYLATION_STATUS', 'METHYLATION_BASELINE_BETA', 'METHYLATION_HYPER_SAMPLES',
    'METHYLATION_HYPER_FREQUENCY', 'METHYLATION_FUNCTIONAL', 'QVALUE_METHYLATION',
    'SAMPLES_MUTATION_METHYLATION', 'SAMPLES_MUTATED_AND_HYPERMETHYLATED', 'SAMPLES_MUTATED_OR_HYPERMETHYLATED',
    'OMICS_ROLE_SUPPORT'
]


def read_mutations(path):
    """Read the processed VEP output of a cohort (parse-vep)"""
    df = pd.read_csv(path, sep='\t', usecols=['#Uploaded_variation', 'Consequence', 'SYMBOL'], dtype=str)
    df['SAMPLE'] = df['#Uploaded_variation'].str.split('__').str[1]
    return df.dropna(subset=['SYMBOL', 'SAMPLE'])


def mutated_samples(mutations):
    """
    Returns:
        (all samples, dict gene -> non-synonymously mutated samples,
        dict gene -> samples with truncating mutations)
    """
    nonsyn = mutations[~mutations['Consequence'].isin(SYNONYMOUS)]
    mutated = nonsyn.groupby('SYMBOL')['SAMPLE'].apply(set).to_dict()
    truncating = nonsyn[nonsyn['Consequence'].isin(TRUNCATING)].groupby('SYMBOL')['SAMPLE'].apply(set).to_dict()
    return set(mutations['SAMPLE']), mutated, truncating


def _compare(values, is_mutant, min_samples):
    mut, wt = values[is_mutant], values[~is_mutant]
    mut, wt = mut[np.isfinite(mut)], wt[np.isfinite(wt)]
    if len(mut) < min_samples or len(wt) < min_samples:
        return len(mut), np.nan, np.nan
    return len(mut), float(np.median(mut) - np.median(wt)), mannwhitney(mut, wt)


def expression_effects(expression, cohort_samples, mutated, truncating, min_samples=2):
    """Expression of mutated vs wild-type tumours for each mutated gene"""
    shared = [s for s in expression.columns if s in cohort_samples]
    columns = ['SYMBOL', 'EXPRESSION_MUTANT_SAMPLES', 'EXPRESSION_LOG2FC_MUTANT', 'PVALUE_EXPRESSION_MUTANT',
               'EXPRESSION_TRUNCATING_SAMPLES', 'EXPRESSION_LOG2FC_TRUNCATING', 'PVALUE_EXPRESSION_TRUNCATING']
    if len(shared) == 0:
        return pd.DataFrame(columns=columns + ['QVALUE_EXPRESSION_MUTANT']), 0

    E = expression[shared]
    genes = set(E.index)
    shared = np.array(shared)
    rows = []
    for gene, samples in mutated.items():
        if gene not in genes:
            continue
        values = E.loc[gene].values.astype(float)
        is_mutant = np.isin(shared, list(samples))
        n_mut, lfc, p = _compare(values, is_mutant, min_samples)

        trunc = truncating.get(gene, set())
        # wild-type: tumours without any non-synonymous mutation in the gene
        is_trunc = np.isin(shared, list(trunc))
        keep = is_trunc | ~is_mutant
        n_trunc, lfc_trunc, p_trunc = _compare(values[keep], is_trunc[keep], min_samples)
        rows.append([gene, n_mut, lfc, p, n_trunc, lfc_trunc, p_trunc])

    df = pd.DataFrame(rows, columns=columns)
    df['QVALUE_EXPRESSION_MUTANT'] = fdr_bh(df['PVALUE_EXPRESSION_MUTANT'].values)
    return df, len(shared)


def methylation_hits(events, methylated_samples, cohort_samples, mutated):
    """Tumours with a mutation and/or promoter hypermethylation of each gene"""
    profiled = set(methylated_samples) & cohort_samples
    hyper = events.groupby('SYMBOL')['SAMPLE'].apply(set).to_dict() if len(events) else {}
    rows = []
    for gene in set(hyper) | set(mutated):
        h = hyper.get(gene, set()) & profiled
        m = mutated.get(gene, set()) & profiled
        rows.append([gene, len(profiled), len(h & m), len(h | m)])
    columns = ['SYMBOL', 'SAMPLES_MUTATION_METHYLATION', 'SAMPLES_MUTATED_AND_HYPERMETHYLATED',
               'SAMPLES_MUTATED_OR_HYPERMETHYLATED']
    return pd.DataFrame(rows, columns=columns), len(profiled)


def as_bool(value):
    """Parse booleans read from TSV files (True/False/NA)"""
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, str) and value.lower() in ('true', 'false'):
        return value.lower() == 'true'
    return None


def role_support(row, threshold=0.1):
    def significant(column):
        value = row.get(column, np.nan)
        return pd.notna(value) and value < threshold

    lof = (significant('QVALUE_METHYLATION') and as_bool(row.get('METHYLATION_FUNCTIONAL')) is not False) or \
          (significant('QVALUE_EXPRESSION') and row.get('EXPRESSION_OUTLIER_DIRECTION') == 'under')
    act = (significant('QVALUE_EXPRESSION') and row.get('EXPRESSION_OUTLIER_DIRECTION') == 'over') or \
          (significant('QVALUE_EXPRESSION_MUTANT') and row.get('EXPRESSION_LOG2FC_MUTANT', 0) > 0)
    if lof and act:
        return 'conflicting'
    if lof:
        return 'LoF'
    if act:
        return 'Act'
    return None


def run(mutations_file, output, expression_matrix=None, expression_results=None,
        methylation_matrix=None, methylation_results=None, methylation_events=None,
        min_samples=2, threshold=0.1):

    stats = {}
    mutations = read_mutations(mutations_file)
    cohort_samples, mutated, truncating = mutated_samples(mutations)
    stats['samples_mutations'] = len(cohort_samples)

    tables = [pd.DataFrame({'SYMBOL': sorted(mutated)})]

    if expression_results is not None:
        expr = pd.read_csv(expression_results, sep='\t')
        expr = expr.rename(columns={
            'MEDIAN_LOG2': 'MEDIAN_LOG2_EXPRESSION',
            'DIRECTION': 'EXPRESSION_OUTLIER_DIRECTION',
            'Q_VALUE': 'QVALUE_EXPRESSION'})
        expr['EXPRESSION_OUTLIER_SAMPLES'] = np.where(
            expr['EXPRESSION_OUTLIER_DIRECTION'] == 'under', expr['UNDER_SAMPLES'], expr['OVER_SAMPLES'])
        tables.append(expr[['SYMBOL', 'EXPRESSED', 'FRACTION_EXPRESSED', 'MEDIAN_LOG2_EXPRESSION',
                            'EXPRESSION_OUTLIER_DIRECTION', 'EXPRESSION_OUTLIER_SAMPLES', 'QVALUE_EXPRESSION']])

    if expression_matrix is not None:
        effects, n = expression_effects(io.read_gene_matrix(expression_matrix), cohort_samples,
                                        mutated, truncating, min_samples=min_samples)
        stats['samples_mutations_expression'] = n
        if n == 0:
            stats['warning_expression_samples'] = 'No expression sample matches the samples with mutations'
        tables.append(effects)

    if methylation_results is not None:
        meth = pd.read_csv(methylation_results, sep='\t')
        meth = meth.rename(columns={
            'STATUS': 'METHYLATION_STATUS',
            'BASELINE_BETA': 'METHYLATION_BASELINE_BETA',
            'HYPER_SAMPLES': 'METHYLATION_HYPER_SAMPLES',
            'HYPER_FREQUENCY': 'METHYLATION_HYPER_FREQUENCY',
            'FUNCTIONAL': 'METHYLATION_FUNCTIONAL',
            'Q_VALUE': 'QVALUE_METHYLATION'})
        tables.append(meth[['SYMBOL', 'METHYLATION_STATUS', 'METHYLATION_BASELINE_BETA', 'METHYLATION_HYPER_SAMPLES',
                            'METHYLATION_HYPER_FREQUENCY', 'METHYLATION_FUNCTIONAL', 'QVALUE_METHYLATION']])

    if methylation_events is not None and methylation_matrix is not None:
        events = pd.read_csv(methylation_events, sep='\t', usecols=['SYMBOL', 'SAMPLE'], dtype=str)
        samples = pd.read_csv(methylation_matrix, sep='\t', nrows=0).columns[1:]
        hits, n = methylation_hits(events, samples, cohort_samples, mutated)
        stats['samples_mutations_methylation'] = n
        if n == 0:
            stats['warning_methylation_samples'] = 'No methylation sample matches the samples with mutations'
        tables.append(hits)

    df = tables[0]
    for t in tables[1:]:
        df = df.merge(t, on='SYMBOL', how='outer')
    for c in OMICS_COLUMNS:
        if c not in df.columns:
            df[c] = np.nan
    df['OMICS_ROLE_SUPPORT'] = [role_support(r, threshold) for r in df.to_dict('records')]
    df = df[OMICS_COLUMNS].sort_values('SYMBOL')

    stats['genes'] = int(len(df))
    stats['role_support'] = df['OMICS_ROLE_SUPPORT'].value_counts().to_dict()
    io.write_table(df, output)
    io.write_stats(output + '.stats.json', stats)


def summary(files, output, threshold=0.1):
    """
    Genes with significant epigenetic silencing or expression outliers
    in any cohort (candidate drivers altered by non-mutational mechanisms).
    """
    dfs = []
    for file in files:
        # booleans as text: concatenating them with missing values would turn them into numbers
        df = pd.read_csv(file, sep='\t', dtype={'EXPRESSED': str, 'METHYLATION_FUNCTIONAL': str})
        significant = (df['QVALUE_METHYLATION'] < threshold) | (df['QVALUE_EXPRESSION'] < threshold)
        df = df[significant].copy()
        df.insert(1, 'COHORT', os.path.basename(file).split('.')[0])
        dfs.append(df)
    columns = ['SYMBOL', 'COHORT'] + OMICS_COLUMNS[1:]
    df = pd.concat(dfs, sort=False) if len(dfs) > 0 else pd.DataFrame(columns=columns)
    io.write_table(df[columns].sort_values(['SYMBOL', 'COHORT']), output)


@click.command()
@click.option('-o', '--output', type=click.Path(), required=True)
@click.option('--threshold', type=float, default=0.1, show_default=True)
@click.argument('files', nargs=-1)
def summary_cli(output, threshold, files):
    summary(files, output, threshold=threshold)


@click.command()
@click.option('-m', '--mutations', 'mutations_file', type=click.Path(exists=True), required=True,
              help='Processed VEP output of the cohort (parse-vep)')
@click.option('-o', '--output', type=click.Path(), required=True)
@click.option('--expression-matrix', type=click.Path(exists=True), default=None)
@click.option('--expression', 'expression_results', type=click.Path(exists=True), default=None)
@click.option('--methylation-matrix', type=click.Path(exists=True), default=None)
@click.option('--methylation', 'methylation_results', type=click.Path(exists=True), default=None)
@click.option('--methylation-events', type=click.Path(exists=True), default=None)
@click.option('--min-samples', type=int, default=2, show_default=True,
              help='Minimum mutated (and wild-type) tumours to compare their expression')
@click.option('--threshold', type=float, default=0.1, show_default=True,
              help='Q-value threshold for the mode of action support')
def cli(**kwargs):
    run(**kwargs)
