"""
Transcriptomics layer.

``parse-expression``
    Normalises an RNA-seq matrix (raw counts, TPM/FPKM or log-transformed
    values) into log2(x + 1) values of the MANE protein coding genes
    (x = CPM for counts, the provided values for TPM/FPKM) of the tumour
    samples.

``expression-analysis``
    Per gene:

    - expression status: a gene is not expressed in the cohort when at least
      80% of the tumours have an expression below 1 (TPM/CPM/FPKM). This is
      the same criterion used with the TCGA data when no cohort expression
      is available.
    - expression outliers: tumours with a modified z-score above 3.5
      (Iglewicz and Hoaglin) and at least a 2-fold change with respect to
      the median of the cohort. The number of tumours with an over- (or
      under-) expression outlier is compared to the expectation given the
      genome-wide outlier rate of each tumour (Poisson-binomial recurrence
      test), in the spirit of the Cancer Outlier Profile Analysis.
"""

import warnings

import click
import numpy as np
import pandas as pd

from intogen_core.omics import io
from intogen_core.omics.stats import fdr_bh, recurrence_test, sample_background, robust_scale


RESULT_COLUMNS = [
    'SYMBOL', 'STATUS', 'N_SAMPLES', 'MEDIAN_LOG2', 'FRACTION_EXPRESSED', 'EXPRESSED',
    'OVER_SAMPLES', 'UNDER_SAMPLES', 'EXPECTED_OVER', 'EXPECTED_UNDER',
    'P_OVER', 'P_UNDER', 'DIRECTION', 'P_VALUE', 'Q_VALUE'
]

EVENT_COLUMNS = ['SYMBOL', 'SAMPLE', 'LOG2_EXPRESSION', 'ZSCORE', 'DIRECTION']

# Genes whose expression mostly reflects the sex of the patient
OUTLIER_EXCLUDED_CHROMOSOMES = {'Y', 'MT', 'M'}

COUNTS_SUMMARY_ROWS = r'^(__|N_(unmapped|multimapping|noFeature|ambiguous)$)'


def detect_units(values):
    """Guess the units of an expression matrix"""
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        raise io.OmicsError('The expression matrix has no values')
    if finite.min() < 0:
        return 'log2'
    if np.allclose(finite, np.round(finite)):
        return 'counts'
    if finite.max() <= 30:
        return 'log2'
    library_size = np.nanmedian(np.nansum(values, axis=0))
    if library_size > 5e6:
        return 'counts'  # e.g. RSEM expected counts
    return 'tpm'


def normalise(values, units):
    """Log2 transform (counts are first converted to counts per million)"""
    if units == 'counts':
        library = np.nansum(values, axis=0)
        if np.any(library <= 0):
            raise io.OmicsError('There are samples without counts')
        values = values / library * 1e6
    if units in ('counts', 'tpm', 'fpkm'):
        if np.nanmin(values) < 0:
            raise io.OmicsError(f'Negative values are not valid {units}')
        values = np.log2(values + 1)
    return values


def parse(input_file, output, samples=None, units='auto', gene_annotation=None, symbol_map=None):
    io.check_outputs(input_file, output)
    stats = {}
    sheet = io.read_samplesheet(samples)
    mapper = io.GeneMapper(io.load_gene_annotation(gene_annotation), io.load_symbol_map(symbol_map))

    df, dropped, total = io.read_matrix(input_file, dtype=np.float64)
    stats['genes_total'] = total
    stats['dropped_columns'] = dropped
    # summary rows of HTSeq (__no_feature...) and STAR (N_unmapped...) are not genes
    summary_rows = df.index.str.match(COUNTS_SUMMARY_ROWS)
    if summary_rows.any():
        stats['summary_rows_removed'] = list(df.index[summary_rows])
        df = df[~summary_rows]

    units = detect_units(df.values) if units == 'auto' else units
    stats['units'] = units
    # normalise before selecting genes so that library sizes include all genes
    df = pd.DataFrame(normalise(df.values, units), index=df.index, columns=df.columns)

    symbols = df.index.map(mapper.map)
    stats['genes_unmapped'] = int(pd.isna(symbols).sum())
    df = df[pd.notna(symbols)]
    df.index = symbols[pd.notna(symbols)]
    # for duplicated genes keep the most expressed entry
    duplicated = df.index.duplicated(keep=False)
    if duplicated.any():
        stats['genes_duplicated'] = int(df.index[duplicated].nunique())
        order = np.argsort(-df.mean(axis=1).values, kind='mergesort')
        df = df.iloc[order]
        df = df[~df.index.duplicated(keep='first')]

    tumors, normals, discarded = io.assign_samples(df.columns, sheet)
    if len(tumors) == 0:
        raise io.OmicsError('There are no tumour samples in the expression matrix')
    if len(discarded) > 0:
        stats['warning_duplicated_samples'] = discarded
    stats['normals_ignored'] = len(normals)

    out = df[list(tumors.values())]
    out.columns = list(tumors.keys())
    out = out[out.notna().any(axis=1)].sort_index()

    stats['genes'] = int(out.shape[0])
    stats['tumors'] = len(tumors)
    io.write_matrix(out, output)
    io.write_stats(output + '.stats.json', stats)


def analyse(matrix, output, events_output, gene_annotation=None, expressed_min=1.0,
            non_expressed_fraction=0.8, outlier_z=3.5, min_log2_change=1.0, min_samples=10):

    df = io.read_gene_matrix(matrix)
    genes, samples = df.index.values, df.columns.values
    X = df.values
    observed = np.isfinite(X)
    n = observed.sum(axis=1)
    threshold = np.log2(expressed_min + 1)

    n_expressing = ((X >= threshold) & observed).sum(axis=1)
    fraction_expressed = np.where(n > 0, n_expressing / np.maximum(n, 1), np.nan)
    # not expressed: at least `non_expressed_fraction` of the samples below the threshold
    not_expressed = (n - n_expressing) >= non_expressed_fraction * n - 1e-9
    expressed = pd.array(~not_expressed, dtype='boolean')
    expressed[n == 0] = pd.NA

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)
        median, scale = robust_scale(X)
        z = (X - median[:, None]) / np.where(scale > 0, scale, np.nan)[:, None]

    excluded = np.zeros(len(genes), dtype=bool)
    annotation = io.load_gene_annotation(gene_annotation)
    if annotation is not None:
        chromosome = dict(zip(annotation['SYMBOL'], annotation['CHROMOSOME']))
        excluded = np.array([chromosome.get(g, '') in OUTLIER_EXCLUDED_CHROMOSOMES for g in genes], dtype=bool)

    testable = (n >= min_samples) & (scale > 0) & ~excluded
    status = np.where(excluded, 'sex_chromosome',
                      np.where(n < min_samples, 'insufficient_samples',
                               np.where(scale > 0, 'tested', 'invariant')))

    Xf = np.nan_to_num(X, nan=0.0)
    zf = np.nan_to_num(z, nan=0.0)
    over = observed & (zf > outlier_z) & (Xf - median[:, None] >= min_log2_change) & (Xf >= threshold)
    under = observed & (zf < -outlier_z) & (median[:, None] - Xf >= min_log2_change) & (median[:, None] >= threshold)
    over &= testable[:, None]
    under &= testable[:, None]

    c_over, e_over, p_over = recurrence_test(over, observed, sample_background(over, observed, testable), genes=testable)
    c_under, e_under, p_under = recurrence_test(under, observed, sample_background(under, observed, testable), genes=testable)

    p_value = np.minimum(1.0, 2 * np.fmin(p_over, p_under))
    direction = np.where(~testable, None, np.where(p_over <= p_under, 'over', 'under'))
    q_value = fdr_bh(p_value)

    result = pd.DataFrame({
        'SYMBOL': genes,
        'STATUS': status,
        'N_SAMPLES': n,
        'MEDIAN_LOG2': median,
        'FRACTION_EXPRESSED': fraction_expressed,
        'EXPRESSED': expressed,
        'OVER_SAMPLES': c_over,
        'UNDER_SAMPLES': c_under,
        'EXPECTED_OVER': e_over,
        'EXPECTED_UNDER': e_under,
        'P_OVER': p_over,
        'P_UNDER': p_under,
        'DIRECTION': direction,
        'P_VALUE': p_value,
        'Q_VALUE': q_value,
    }, columns=RESULT_COLUMNS)
    result.sort_values(['P_VALUE', 'SYMBOL'], inplace=True, na_position='last')
    io.write_table(result, output)

    events = over | under
    gi, si = np.nonzero(events)
    ev = pd.DataFrame({
        'SYMBOL': genes[gi],
        'SAMPLE': samples[si],
        'LOG2_EXPRESSION': X[gi, si],
        'ZSCORE': z[gi, si],
        'DIRECTION': np.where(over[gi, si], 'over', 'under'),
    }, columns=EVENT_COLUMNS)
    io.write_table(ev.sort_values(['SYMBOL', 'SAMPLE']), events_output)

    stats = {
        'genes': int(len(genes)),
        'genes_tested': int(testable.sum()),
        'genes_not_expressed': int((~expressed.fillna(True)).sum()),
        'tumors': int(len(samples)),
        'over_events': int(over.sum()),
        'under_events': int(under.sum()),
        'genes_q_lt_0.1': int(np.nansum(q_value < 0.1)),
    }
    io.write_stats(output + '.stats.json', stats)


@click.command()
@click.option('-i', '--input', 'input_file', type=click.Path(exists=True), required=True,
              help='Expression matrix (genes x samples)')
@click.option('-o', '--output', type=click.Path(), required=True, help='Normalised tumour expression matrix')
@click.option('--samples', type=click.Path(exists=True), default=None, help='Sample sheet (ID, SAMPLE, TYPE)')
@click.option('--units', type=click.Choice(['auto', 'counts', 'tpm', 'fpkm', 'log2']), default='auto',
              help='Units of the expression values. log2 values are used as provided')
@click.option('--genes', 'gene_annotation', type=click.Path(exists=True), default=None,
              help='Gene annotation. Default: $INTOGEN_DATASETS/regions/cds_biomart.tsv')
@click.option('--symbols', 'symbol_map', type=click.Path(exists=True), default=None,
              help='Outdated to current HUGO symbols (JSON)')
def parse_cli(input_file, output, samples, units, gene_annotation, symbol_map):
    parse(input_file, output, samples=samples, units=units,
          gene_annotation=gene_annotation or io.default_dataset('regions', 'cds_biomart.tsv'),
          symbol_map=symbol_map or io.default_dataset('others', 'mapping_new_hugo_symbols.json'))


@click.command()
@click.option('-i', '--input', 'matrix', type=click.Path(exists=True), required=True,
              help='Expression matrix (parse-expression output)')
@click.option('-o', '--output', type=click.Path(), required=True)
@click.option('--events', 'events_output', type=click.Path(), required=True)
@click.option('--genes', 'gene_annotation', type=click.Path(exists=True), default=None,
              help='Gene annotation. Default: $INTOGEN_DATASETS/regions/cds_biomart.tsv')
@click.option('--expressed-min', type=float, default=1.0, show_default=True,
              help='Minimum expression (TPM/CPM/FPKM) of an expressing sample')
@click.option('--non-expressed-fraction', type=float, default=0.8, show_default=True,
              help='Fraction of non-expressing samples of a non-expressed gene')
@click.option('--outlier-z', type=float, default=3.5, show_default=True, help='Modified z-score of an outlier')
@click.option('--min-log2-change', type=float, default=1.0, show_default=True,
              help='Minimum log2 difference between an outlier and the median')
@click.option('--min-samples', type=int, default=10, show_default=True)
def analysis_cli(gene_annotation, **kwargs):
    analyse(gene_annotation=gene_annotation or io.default_dataset('regions', 'cds_biomart.tsv'), **kwargs)
