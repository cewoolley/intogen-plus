"""
DNA methylation layer.

``parse-methylation``
    Summarises array probes (Illumina 450K, EPIC or EPICv2; beta values,
    M-values or percentages) into the promoter methylation (beta value) of
    each gene, using the probes around the TSS of its MANE transcript.
    Gene-level matrices are also accepted. Tumour and normal samples are
    written into separate matrices.

``methylation-analysis``
    Detects genes recurrently silenced by promoter hypermethylation:

    1. The reference (unmethylated) state of each promoter is the mean beta
       value in normal samples, or, when not enough normal samples are
       available, a low quantile of the tumour beta values. Only genes with
       an unmethylated reference promoter are tested.
    2. A tumour is hypermethylated at a gene when its promoter beta value is
       above both an absolute threshold and the reference plus a minimum
       difference.
    3. The number of hypermethylated tumours is compared with the expectation
       given the genome-wide hypermethylation rate of each tumour
       (Poisson-binomial recurrence test).
    4. When expression data are available, a gene is considered
       functionally silenced only if hypermethylated tumours show
       lower expression than the rest (one-sided Mann-Whitney test).
       Genes with evidence against silencing get a p-value of 1.
"""

import re
import warnings

import click
import numpy as np
import pandas as pd
from scipy import sparse

from intogen_core.omics import io
from intogen_core.omics.stats import fdr_bh, recurrence_test, sample_background, \
    mannwhitney, spearman


PROBE_ID = re.compile(r'^(cg\d+|ch\.\w+|rs\d+|nv_\w+)')
EPICV2_SUFFIX = re.compile(r'_[TB][CO]\d+$')

RESULT_COLUMNS = [
    'SYMBOL', 'STATUS', 'N_SAMPLES', 'N_NORMALS', 'BASELINE_BETA', 'BASELINE_SOURCE',
    'MEDIAN_BETA', 'THRESHOLD_BETA', 'HYPER_SAMPLES', 'HYPER_FREQUENCY', 'EXPECTED_HYPER',
    'MEAN_DELTA_BETA', 'P_RECURRENCE', 'N_SAMPLES_EXPRESSION', 'SPEARMAN_RHO',
    'LOG2FC_SILENCING', 'P_SILENCING', 'FUNCTIONAL', 'P_VALUE', 'Q_VALUE'
]

EVENT_COLUMNS = ['SYMBOL', 'SAMPLE', 'BETA', 'BASELINE_BETA', 'DELTA_BETA', 'LOG2_EXPRESSION']


def strip_suffix(probe):
    """Remove the EPICv2 replicate suffix (e.g. cg00000029_TC21 -> cg00000029)"""
    return EPICV2_SUFFIX.sub('', probe)


def is_probe_level(identifiers):
    identifiers = [str(i) for i in identifiers]
    if len(identifiers) == 0:
        return False
    matches = sum(1 for i in identifiers if PROBE_ID.match(i))
    return matches >= 0.5 * len(identifiers)


class ProbeIndex:
    """Map array probes to the genes whose promoter they cover"""

    def __init__(self, path):
        df = pd.read_csv(path, sep='\t', dtype=str)
        df = df.dropna(subset=['PROBE', 'SYMBOL'])
        self.chromosome = {}
        if 'CHROMOSOME' in df.columns:
            self.chromosome = dict(zip(df['SYMBOL'], df['CHROMOSOME']))
        df = df[['PROBE', 'SYMBOL']].drop_duplicates()
        self.exact = df.groupby('PROBE')['SYMBOL'].apply(list).to_dict()
        stripped = df.assign(PROBE=df['PROBE'].map(strip_suffix)).drop_duplicates()
        self.stripped = stripped.groupby('PROBE')['SYMBOL'].apply(list).to_dict()

    def genes(self, probe):
        genes = self.exact.get(probe)
        if genes is None:
            genes = self.stripped.get(strip_suffix(probe))
        return genes

    def contains(self, probes):
        return np.array([self.genes(p) is not None for p in probes], dtype=bool)


def detect_values(values):
    """Guess whether the values are beta values, percentages or M-values"""
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        raise io.OmicsError('The methylation matrix has no values')
    lo, hi = finite.min(), finite.max()
    if lo >= 0 and hi <= 1:
        return 'beta'
    if lo >= 0 and hi <= 100:
        return 'percent'
    return 'm'


def to_beta(values, kind):
    if kind == 'm':
        values = 1.0 / (1.0 + np.power(2.0, -values))
    elif kind == 'percent':
        values = values / 100.0
    finite = values[np.isfinite(values)]
    if finite.size and (finite.min() < -1e-6 or finite.max() > 1 + 1e-6):
        raise io.OmicsError(f'Methylation values are not in the expected range for "{kind}"')
    return np.clip(values, 0, 1)


def aggregate_probes(df, index, min_probes=1):
    """
    Average the beta values of the promoter probes of each gene.

    Returns:
        (gene x sample DataFrame, Series with the number of probes per gene)
    """
    rows, cols, genes = [], [], {}
    for j, probe in enumerate(df.index):
        for gene in index.genes(probe) or []:
            rows.append(genes.setdefault(gene, len(genes)))
            cols.append(j)
    names = list(genes)
    if len(names) == 0:
        raise io.OmicsError('None of the probes maps to a gene promoter')

    A = sparse.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)),
                          shape=(len(names), df.shape[0]))
    values = df.values.astype(np.float32, copy=False)
    finite = np.isfinite(values)
    sums = np.asarray(A @ np.where(finite, values, np.float32(0)), dtype=float)
    counts = np.asarray(A @ finite.astype(np.float32), dtype=float)
    with np.errstate(invalid='ignore', divide='ignore'):
        mean = sums / counts
    mean[counts < min_probes] = np.nan

    n_probes = pd.Series(np.asarray(A.sum(axis=1)).ravel().astype(int), index=names)
    return pd.DataFrame(mean, index=names, columns=df.columns), n_probes


def parse(input_file, output, output_normal, probes=None, samples=None, values='auto',
          gene_annotation=None, symbol_map=None, min_probes=1, keep_sex_chromosomes=False):

    io.check_outputs(input_file, output, output_normal)
    stats = {}
    sheet = io.read_samplesheet(samples)
    mapper = io.GeneMapper(io.load_gene_annotation(gene_annotation), io.load_symbol_map(symbol_map))

    first = pd.read_csv(input_file, sep=io._separator(input_file), usecols=[0], nrows=5000, dtype=str)
    probe_level = is_probe_level(first.iloc[:, 0].dropna())
    stats['input_level'] = 'probe' if probe_level else 'gene'

    chromosome = mapper.chromosome
    if probe_level:
        if probes is None:
            raise io.OmicsError('Probe-level methylation requires the promoter probes annotation (--probes)')
        index = ProbeIndex(probes)
        df, dropped, total = io.read_matrix(input_file, keep=index.contains)
        stats['probes_total'] = total
        stats['probes_promoter'] = int(df.shape[0])
        if df.shape[0] == 0:
            raise io.OmicsError('None of the probes maps to a gene promoter')
        chromosome = index.chromosome or chromosome
    else:
        df, dropped, total = io.read_matrix(input_file)
        symbols = df.index.map(mapper.map)
        stats['genes_total'] = total
        stats['genes_unmapped'] = int(pd.isna(symbols).sum())
        df = df[pd.notna(symbols)]
        df.index = symbols[pd.notna(symbols)]
        df = df.groupby(level=0).mean()

    stats['dropped_columns'] = dropped

    kind = detect_values(df.values) if values == 'auto' else values
    stats['values'] = kind
    # values are kept as float32 (large arrays)
    df = pd.DataFrame(to_beta(df.values.astype(np.float32, copy=False), kind), index=df.index, columns=df.columns)

    if probe_level:
        genes, n_probes = aggregate_probes(df, index, min_probes=min_probes)
        stats['probes_per_gene_median'] = float(n_probes.median())
    else:
        genes = df

    if not keep_sex_chromosomes:
        if chromosome:
            sex = [g for g in genes.index if chromosome.get(g, '') in io.SEX_CHROMOSOMES]
            genes = genes.drop(index=sex)
            stats['genes_sex_chromosomes_removed'] = len(sex)
        else:
            stats['warning_sex_chromosomes'] = 'Chromosome annotation unavailable: genes in sex chromosomes were not removed'

    tumors, normals, discarded = io.assign_samples(genes.columns, sheet)
    if len(tumors) == 0:
        raise io.OmicsError('There are no tumour samples in the methylation matrix')
    if len(discarded) > 0:
        stats['warning_duplicated_samples'] = discarded

    tumor_df = genes[list(tumors.values())]
    tumor_df.columns = list(tumors.keys())
    tumor_df = tumor_df[tumor_df.notna().any(axis=1)].sort_index()

    normal_df = genes[list(normals.values())].reindex(tumor_df.index)
    normal_df.columns = list(normals.keys())

    stats['genes'] = int(tumor_df.shape[0])
    stats['tumors'] = len(tumors)
    stats['normals'] = len(normals)

    io.write_matrix(tumor_df, output)
    io.write_matrix(normal_df, output_normal)
    io.write_stats(output + '.stats.json', stats)


def analyse(tumor_file, output, events_output, normal_file=None, expression_file=None,
            unmethylated_max=0.2, hypermethylated_min=0.3, min_delta=0.2, baseline_quantile=0.25,
            min_normals=3, min_samples=10, min_expression_samples=3, silencing_alpha=0.05):

    stats = {}
    tumors = io.read_gene_matrix(tumor_file)
    genes, samples = tumors.index.values, tumors.columns.values
    X = tumors.values
    observed = np.isfinite(X)
    n_samples = observed.sum(axis=1)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)
        baseline_tumor = np.nanquantile(X, baseline_quantile, axis=1)
        median_beta = np.nanmedian(X, axis=1)

        n_normals = np.zeros(len(genes), dtype=int)
        baseline_normal = np.full(len(genes), np.nan)
        if normal_file is not None:
            normals = io.read_gene_matrix(normal_file).reindex(tumors.index)
            if normals.shape[1] > 0:
                n_normals = np.isfinite(normals.values).sum(axis=1)
                baseline_normal = np.nanmean(normals.values, axis=1)

    use_normals = n_normals >= min_normals
    baseline = np.where(use_normals, baseline_normal, baseline_tumor)
    source = np.where(use_normals, 'normal', 'tumor_quantile')
    threshold = np.maximum(hypermethylated_min, baseline + min_delta)

    enough = n_samples >= min_samples
    eligible = enough & (baseline < unmethylated_max)
    status = np.where(~enough, 'insufficient_samples',
                      np.where(eligible, 'tested', 'baseline_methylated'))

    events = observed & (np.nan_to_num(X, nan=-1.0) >= threshold[:, None]) & eligible[:, None]
    background = sample_background(events, observed, eligible)
    counts, expected, p_recurrence = recurrence_test(events, observed, background, genes=eligible)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)
        delta = np.where(events, X - baseline[:, None], np.nan)
        mean_delta = np.nanmean(delta, axis=1)

    # Functional silencing: expression of hypermethylated vs other tumours
    n_expression = np.zeros(len(genes), dtype=int)
    rho = np.full(len(genes), np.nan)
    lfc = np.full(len(genes), np.nan)
    p_silencing = np.full(len(genes), np.nan)
    functional = np.full(len(genes), None, dtype=object)
    expression_values = None
    if expression_file is not None:
        expression = io.read_gene_matrix(expression_file)
        shared = [s for s in samples if s in set(expression.columns)]
        stats['samples_with_expression'] = len(shared)
        if len(shared) > 0:
            expression = expression.reindex(index=tumors.index, columns=shared)
            cols = np.array([tumors.columns.get_loc(s) for s in shared])
            E = expression.values
            expression_values = (expression, cols)
            for i in np.flatnonzero(eligible):
                x, e = X[i, cols], E[i]
                ok = np.isfinite(x) & np.isfinite(e)
                n_expression[i] = ok.sum()
                if ok.sum() < min_samples:
                    continue
                rho[i] = spearman(x[ok], e[ok])
                hyper = events[i, cols] & ok
                other = ~events[i, cols] & ok
                if hyper.sum() >= min_expression_samples and other.sum() >= min_expression_samples:
                    p_silencing[i] = mannwhitney(e[hyper], e[other], alternative='less')
                    lfc[i] = np.median(e[hyper]) - np.median(e[other])
                    functional[i] = bool(p_silencing[i] < silencing_alpha)
    else:
        stats['samples_with_expression'] = 0

    # Genes with evidence against silencing do not contribute as candidates
    p_value = p_recurrence.copy()
    not_functional = np.array([f is False for f in functional])
    p_value[not_functional] = 1.0
    q_value = fdr_bh(np.where(eligible, p_value, np.nan))

    with np.errstate(invalid='ignore', divide='ignore'):
        frequency = np.where(n_samples > 0, counts / np.maximum(n_samples, 1), np.nan)

    result = pd.DataFrame({
        'SYMBOL': genes,
        'STATUS': status,
        'N_SAMPLES': n_samples,
        'N_NORMALS': n_normals,
        'BASELINE_BETA': baseline,
        'BASELINE_SOURCE': source,
        'MEDIAN_BETA': median_beta,
        'THRESHOLD_BETA': threshold,
        'HYPER_SAMPLES': np.where(eligible, counts, 0),
        'HYPER_FREQUENCY': np.where(eligible, frequency, 0.0),
        'EXPECTED_HYPER': expected,
        'MEAN_DELTA_BETA': mean_delta,
        'P_RECURRENCE': p_recurrence,
        'N_SAMPLES_EXPRESSION': n_expression,
        'SPEARMAN_RHO': rho,
        'LOG2FC_SILENCING': lfc,
        'P_SILENCING': p_silencing,
        'FUNCTIONAL': pd.array(functional, dtype='boolean'),
        'P_VALUE': p_value,
        'Q_VALUE': q_value,
    }, columns=RESULT_COLUMNS)
    result.sort_values(['P_VALUE', 'SYMBOL'], inplace=True, na_position='last')
    io.write_table(result, output)

    # Events: one row per hypermethylated gene and tumour
    gi, si = np.nonzero(events)
    ev = pd.DataFrame({
        'SYMBOL': genes[gi],
        'SAMPLE': samples[si],
        'BETA': X[gi, si],
        'BASELINE_BETA': baseline[gi],
        'DELTA_BETA': X[gi, si] - baseline[gi],
        'LOG2_EXPRESSION': np.nan,
    }, columns=EVENT_COLUMNS)
    if expression_values is not None and len(ev) > 0:
        expression, cols = expression_values
        position = np.full(len(samples), -1)
        position[cols] = np.arange(len(cols))
        p = position[si]
        E = expression.values
        ev['LOG2_EXPRESSION'] = np.where(p >= 0, E[gi, np.maximum(p, 0)], np.nan)
    io.write_table(ev.sort_values(['SYMBOL', 'SAMPLE']), events_output)

    q75, q25 = np.percentile(background, [75, 25])
    high = samples[background > q75 + 1.5 * (q75 - q25)]
    stats.update({
        'genes': int(len(genes)),
        'genes_tested': int(eligible.sum()),
        'genes_baseline_methylated': int((status == 'baseline_methylated').sum()),
        'genes_insufficient_samples': int((status == 'insufficient_samples').sum()),
        'genes_baseline_from_normals': int((use_normals & eligible).sum()),
        'tumors': int(len(samples)),
        'events': int(events.sum()),
        'hypermethylation_rate_median': float(np.median(background)) if len(background) else None,
        'samples_high_hypermethylation_rate': list(high),
        'genes_not_functional': int(not_functional.sum()),
        'genes_q_lt_0.1': int(np.nansum(q_value < 0.1)),
    })
    io.write_stats(output + '.stats.json', stats)


@click.command()
@click.option('-i', '--input', 'input_file', type=click.Path(exists=True), required=True,
              help='Methylation matrix (probes or genes x samples)')
@click.option('-o', '--output', type=click.Path(), required=True, help='Tumour promoter methylation matrix')
@click.option('--output-normal', type=click.Path(), required=True, help='Normal promoter methylation matrix')
@click.option('--probes', type=click.Path(exists=True), default=None,
              help='Promoter probes annotation. Default: $INTOGEN_DATASETS/methylation/promoter_probes.tsv.gz')
@click.option('--samples', type=click.Path(exists=True), default=None, help='Sample sheet (ID, SAMPLE, TYPE)')
@click.option('--values', type=click.Choice(['auto', 'beta', 'm', 'percent']), default='auto',
              help='Type of the methylation values')
@click.option('--genes', 'gene_annotation', type=click.Path(exists=True), default=None,
              help='Gene annotation. Default: $INTOGEN_DATASETS/regions/cds_biomart.tsv')
@click.option('--symbols', 'symbol_map', type=click.Path(exists=True), default=None,
              help='Outdated to current HUGO symbols (JSON)')
@click.option('--min-probes', type=int, default=1, help='Minimum promoter probes with data per gene and sample')
@click.option('--keep-sex-chromosomes', is_flag=True, help='Do not remove genes in sex chromosomes')
def parse_cli(input_file, output, output_normal, probes, samples, values, gene_annotation, symbol_map,
              min_probes, keep_sex_chromosomes):
    parse(input_file, output, output_normal,
          probes=probes or io.default_dataset('methylation', 'promoter_probes.tsv.gz'),
          samples=samples, values=values,
          gene_annotation=gene_annotation or io.default_dataset('regions', 'cds_biomart.tsv'),
          symbol_map=symbol_map or io.default_dataset('others', 'mapping_new_hugo_symbols.json'),
          min_probes=min_probes, keep_sex_chromosomes=keep_sex_chromosomes)


@click.command()
@click.option('-i', '--input', 'tumor_file', type=click.Path(exists=True), required=True,
              help='Tumour promoter methylation matrix (parse-methylation output)')
@click.option('--normal', 'normal_file', type=click.Path(exists=True), default=None,
              help='Normal promoter methylation matrix')
@click.option('--expression', 'expression_file', type=click.Path(exists=True), default=None,
              help='Expression matrix (parse-expression output)')
@click.option('-o', '--output', type=click.Path(), required=True)
@click.option('--events', 'events_output', type=click.Path(), required=True)
@click.option('--unmethylated-max', type=float, default=0.2, show_default=True,
              help='Maximum reference beta value of a testable (unmethylated) promoter')
@click.option('--hypermethylated-min', type=float, default=0.3, show_default=True,
              help='Minimum beta value of a hypermethylated promoter')
@click.option('--min-delta', type=float, default=0.2, show_default=True,
              help='Minimum beta value increase over the reference')
@click.option('--baseline-quantile', type=float, default=0.25, show_default=True,
              help='Quantile of the tumour beta values used as reference without normal samples')
@click.option('--min-normals', type=int, default=3, show_default=True)
@click.option('--min-samples', type=int, default=10, show_default=True)
@click.option('--min-expression-samples', type=int, default=3, show_default=True,
              help='Minimum hypermethylated (and other) tumours with expression to test silencing')
@click.option('--silencing-alpha', type=float, default=0.05, show_default=True)
def analysis_cli(**kwargs):
    analyse(**kwargs)
