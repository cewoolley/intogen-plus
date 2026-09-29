"""
Pathway analysis of a cohort (``pathway-analysis``) and summary of all the
cohorts (``pathway-summary``).

For each gene set and layer (mutations, epigenetic silencing, expression
outliers) two scopes are tested:

- ``all``: all the genes of the set.
- ``long_tail``: without the genes that are individually significant in the
  layer (for mutations, genes with tier 1-3 in the combination of driver
  identification methods; for omics layers, genes with q-value < threshold).
  A significant long tail means that the set is altered beyond its known
  drivers, through genes that are rarely altered individually.

The layers of a set are combined with Fisher's method (``combined``): the
mutations, the epigenetic silencing and the expression outliers (the most
significant direction, Bonferroni corrected).

Then, the co-occurrence in the same tumours of the drivers, the
individually significant omics genes and the significant long tails of the
gene sets is tested (see :mod:`intogen_core.pathways.cooccurrence`).
"""

import os

import click
import numpy as np
import pandas as pd

from intogen_core.omics import io
from intogen_core.omics.features import SYNONYMOUS, TRUNCATING, read_mutations
from intogen_core.omics.stats import fdr_bh, fisher_combine
from intogen_core.pathways import cooccurrence
from intogen_core.pathways.genesets import GeneSets, read_gene_sets
from intogen_core.pathways.selection import MUTATION_LAYERS, MutationLayer, expression_layers, silencing_layer


SCOPES = ['all', 'long_tail']

RESULT_COLUMNS = ['SET', 'SOURCE', 'NAME', 'SCOPE', 'LAYER', 'N_GENES', 'N_GENES_EXCLUDED', 'N_GENES_ALTERED',
                  'OBSERVED', 'EXPECTED', 'RATIO', 'TUMOURS', 'P_VALUE', 'Q_VALUE', 'TOP_GENES']

# layer whose tumour alterations are used for co-occurrence, per tested layer
EVENT_LAYER = {'mutation': 'mutation', 'mutation_missense': 'mutation', 'mutation_truncating': 'mutation',
               'silencing': 'silencing', 'expression_over': 'expression_over',
               'expression_under': 'expression_under'}


def read_vet(path):
    """Genes individually significant (tier 1-3) and drivers of the cohort"""
    df = pd.read_csv(path, sep='\t')
    if df.empty or not {'SYMBOL', 'TIER', 'FILTER'} <= set(df.columns):
        return set(), {}
    tier = pd.to_numeric(df['TIER'], errors='coerce')
    significant = set(df.loc[tier <= 3, 'SYMBOL'])
    drivers = df[df['FILTER'] == 'PASS']
    qvalues = pd.to_numeric(drivers['QVALUE_COMBINATION'], errors='coerce').fillna(0) \
        if 'QVALUE_COMBINATION' in drivers.columns else pd.Series(0.0, index=drivers.index)
    return significant, dict(zip(drivers['SYMBOL'], qvalues))


def mutation_hits(mutations_file):
    """Tumours and altered genes per mutation layer"""
    muts = read_mutations(mutations_file)
    tumours = sorted(muts['SAMPLE'].unique())
    nonsyn = muts[~muts['Consequence'].isin(SYNONYMOUS)]
    return tumours, {
        'mutation': nonsyn,
        'mutation_missense': nonsyn[nonsyn['Consequence'] == 'missense_variant'],
        'mutation_truncating': nonsyn[nonsyn['Consequence'].isin(TRUNCATING)],
    }


def events_table(path, layer_genes, direction=None):
    df = pd.read_csv(path, sep='\t', dtype={'SYMBOL': str, 'SAMPLE': str})
    if direction is not None:
        df = df[df['DIRECTION'] == direction]
    return df[df['SYMBOL'].isin(layer_genes)][['SYMBOL', 'SAMPLE']]


def matrix_samples(path):
    return list(pd.read_csv(path, sep='\t', nrows=0).columns[1:])


def top_genes(layer, genes, n=5):
    counts = sorted(((layer.observed(g), g) for g in genes if layer.observed(g) > 0), key=lambda x: (-x[0], x[1]))
    return ','.join(f'{g}:{c}' for c, g in counts[:n])


def load_layers(genemuts, mutations, vet, methylation=None, methylation_events=None, methylation_matrix=None,
                expression=None, expression_events=None, expression_matrix=None, threshold=0.1):
    """
    Returns:
        tests: dict layer -> (set test layer, individually significant genes)
        hits: dict layer -> cooccurrence.Layer (tumour alterations)
        drivers: dict driver -> q-value
        stats: dict
    """
    stats = {}
    tests, hits = {}, {}
    significant, drivers = read_vet(vet)

    genemuts = pd.read_csv(genemuts, sep='\t')
    tumours, mutated = mutation_hits(mutations)
    if len(genemuts) > 0:
        for name in MUTATION_LAYERS:
            layer = MutationLayer(genemuts, name)
            tests[name] = (layer, significant)
            hits[name] = cooccurrence.Layer(name, mutated[name][['SYMBOL', 'SAMPLE']], tumours)
        stats['mutation_method'] = layer.method
        stats['mutation_theta'] = layer.theta
    else:
        stats['warning_mutations'] = 'No dNdScv results: mutations not tested'
    stats['tumours_mutations'] = len(tumours)

    if methylation is not None:
        layer, sig = silencing_layer(pd.read_csv(methylation, sep='\t'), threshold)
        tests['silencing'] = (layer, sig)
        samples = matrix_samples(methylation_matrix)
        hits['silencing'] = cooccurrence.Layer('silencing', events_table(methylation_events, layer.universe), samples)
        stats['tumours_methylation'] = len(samples)

    if expression is not None:
        samples = matrix_samples(expression_matrix)
        for name, (layer, sig) in expression_layers(pd.read_csv(expression, sep='\t'), threshold).items():
            tests[name] = (layer, sig)
            direction = name.split('_')[1]
            hits[name] = cooccurrence.Layer(name, events_table(expression_events, layer.universe, direction), samples)
        stats['tumours_expression'] = len(samples)

    stats['layers'] = list(tests)
    return tests, hits, drivers, stats


def test_sets(gene_sets, tests, hits):
    rows = []
    for set_id in gene_sets.ids:
        genes = gene_sets.genes[set_id]
        info = [set_id, gene_sets.source[set_id], gene_sets.name[set_id]]
        for scope in SCOPES:
            pvalues = {}
            for name, (layer, significant) in tests.items():
                excluded = genes & significant & layer.universe if scope == 'long_tail' else set()
                selected = genes - excluded
                r = layer.test(selected)
                pvalues[name] = r['P_VALUE']
                rows.append(info + [scope, name, r['N_GENES'], len(excluded), r['N_GENES_ALTERED'],
                                    r['OBSERVED'], r['EXPECTED'], r['RATIO'], hits[name].altered_tumours(selected),
                                    r['P_VALUE'], np.nan, top_genes(layer, selected)])
            expression = [pvalues.get(f'expression_{d}', np.nan) for d in ['over', 'under']]
            expression = min(1.0, 2 * np.nanmin(expression)) if np.isfinite(expression).any() else np.nan
            combined = fisher_combine([pvalues.get('mutation', np.nan), pvalues.get('silencing', np.nan), expression])
            rows.append(info + [scope, 'combined', len(genes), np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                                combined, np.nan, ''])

    df = pd.DataFrame(rows, columns=RESULT_COLUMNS)
    for (scope, layer), idx in df.groupby(['SCOPE', 'LAYER']).groups.items():
        df.loc[idx, 'Q_VALUE'] = fdr_bh(df.loc[idx, 'P_VALUE'].values)
    return df


def candidate_events(results, gene_sets, tests, drivers, threshold):
    """Drivers, individually significant omics genes and significant long tails of gene sets"""
    candidates = []
    for gene, q in drivers.items():
        candidates.append(cooccurrence.Event('mutation', 'gene', gene, gene, [gene], q))
    for name in ['silencing', 'expression_over', 'expression_under']:
        if name in tests:
            layer, significant = tests[name]
            for gene in significant:
                candidates.append(cooccurrence.Event(name, 'gene', gene, gene, [gene], 0.0))

    long_tail = results[(results['SCOPE'] == 'long_tail') & (results['LAYER'] != 'combined') &
                        (results['Q_VALUE'] < threshold)]
    best = {}
    for _, row in long_tail.iterrows():
        key = (row['SET'], EVENT_LAYER[row['LAYER']])
        if key not in best or row['Q_VALUE'] < best[key]['Q_VALUE']:
            best[key] = row
    for (set_id, layer_name), row in best.items():
        layer, significant = tests[row['LAYER']]
        genes = (gene_sets.genes[set_id] - significant) & layer.universe
        candidates.append(cooccurrence.Event(layer_name, 'pathway', set_id, row['NAME'], genes, row['Q_VALUE']))
    return candidates


def pathways_by_gene(results, gene_sets, pairs, events, threshold):
    """
    For each gene, the gene sets whose long tail is significant (combined) and
    the co-occurrence modules with an event of the gene
    """
    significant = results[(results['SCOPE'] == 'long_tail') & (results['LAYER'] == 'combined') &
                          (results['Q_VALUE'] < threshold)].sort_values(['Q_VALUE', 'SET'])
    by_gene = {}
    for set_id, name in zip(significant['SET'], significant['NAME']):
        for gene in gene_sets.genes[set_id]:
            by_gene.setdefault(gene, []).append(name)
    pathways = pd.DataFrame({'SYMBOL': list(by_gene), 'PATHWAYS': [';'.join(v) for v in by_gene.values()]},
                            columns=['SYMBOL', 'PATHWAYS'])

    genes = {e.identifier: next(iter(e.genes)) for e in events if e.kind == 'gene'}
    in_module = pairs[pairs['MODULE'].notnull()]
    by_gene = {}
    for column in ['EVENT_1', 'EVENT_2']:
        for event, module in zip(in_module[column], in_module['MODULE']):
            if event in genes:
                by_gene.setdefault(genes[event], set()).add(module)
    modules = pd.DataFrame({'SYMBOL': list(by_gene),
                            'MODULES': [';'.join(sorted(v, key=lambda m: int(m[1:]))) for v in by_gene.values()]},
                           columns=['SYMBOL', 'MODULES'])
    return pathways.merge(modules, on='SYMBOL', how='outer').sort_values('SYMBOL')


def run(genemuts, mutations, vet, gene_sets, output, cooccurrence_output, genes_output,
        min_size=10, max_size=500, threshold=0.1, min_tumours=3, max_events=100, symbol_map=None, **omics):
    tests, hits, drivers, stats = load_layers(genemuts, mutations, vet, threshold=threshold, **omics)

    universe = set().union(*[layer.universe for layer, _ in tests.values()]) if tests else set()
    sets = GeneSets(read_gene_sets(gene_sets), universe, min_size=min_size, max_size=max_size,
                    symbol_map=io.load_symbol_map(symbol_map))
    stats['gene_sets'] = len(sets)
    stats['gene_sets_discarded_by_size'] = sets.discarded

    results = test_sets(sets, tests, hits) if len(sets) and tests else pd.DataFrame(columns=RESULT_COLUMNS)
    results = results.sort_values(['SCOPE', 'LAYER', 'P_VALUE', 'SET'], na_position='last')
    io.write_table(results, output)

    candidates = candidate_events(results, sets, tests, drivers, threshold) if len(results) else []
    event_layers = {name: hits[name] for name in set(EVENT_LAYER.values()) if name in hits}
    pairs, events = cooccurrence.run(candidates, event_layers, max_events=max_events, min_tumours=min_tumours,
                                     threshold=threshold)
    io.write_table(pairs, cooccurrence_output)
    io.write_table(pathways_by_gene(results, sets, pairs, events, threshold), genes_output)

    significant = results[results['Q_VALUE'] < threshold]
    stats['significant'] = {f'{s}:{l}': int(n) for (s, l), n in significant.groupby(['SCOPE', 'LAYER']).size().items()}
    stats['cooccurrence_events'] = len(events)
    stats['cooccurrence_pairs'] = int(len(pairs))
    stats['cooccurrence_testable_pairs'] = int(pairs['TESTABLE'].sum()) if len(pairs) else 0
    stats['cooccurrence_significant_pairs'] = int((pairs['Q_VALUE'] < threshold).sum()) if len(pairs) else 0
    io.write_stats(output + '.stats.json', stats)


def summary(files, pathways_output, cooccurrence_output, modules_output, threshold=0.1):
    """Significant gene sets, co-occurring events and modules of all cohorts"""
    pathways, pairs = [], []
    for file in files:
        cohort = os.path.basename(file).split('.')[0]
        if file.endswith('.pathways.tsv.gz'):
            pathways.append(pd.read_csv(file, sep='\t').assign(COHORT=cohort))
        elif file.endswith('.pathway_cooccurrence.tsv.gz'):
            pairs.append(pd.read_csv(file, sep='\t').assign(COHORT=cohort))

    # gene sets: one row per cohort and set, with the ratio and q-value of each scope and layer
    columns = ['COHORT', 'SET', 'SOURCE', 'NAME', 'N_GENES']
    df = pd.concat(pathways) if pathways else pd.DataFrame(columns=RESULT_COLUMNS + ['COHORT'])
    if len(df):
        df['KEY'] = df['LAYER'].str.upper() + np.where(df['SCOPE'] == 'long_tail', '_LONG_TAIL', '')
        index = ['COHORT', 'SET', 'SOURCE', 'NAME']
        wide = df.pivot_table(index=index, columns='KEY', values=['RATIO', 'Q_VALUE'], aggfunc='first')
        wide.columns = [f'{key}_{value}' if value == 'RATIO' else f'{key}_Q' for value, key in wide.columns]
        wide = wide.reset_index()
        n_genes = df[(df['SCOPE'] == 'all') & (df['LAYER'] == 'combined')][index + ['N_GENES']]
        long_tail_genes = df[(df['SCOPE'] == 'long_tail') & (df['LAYER'] == 'mutation')][index + ['TOP_GENES']]
        wide = wide.merge(n_genes, on=index, how='left').merge(
            long_tail_genes.rename(columns={'TOP_GENES': 'MUTATION_LONG_TAIL_TOP_GENES'}), on=index, how='left')
        q_columns = [c for c in wide.columns if c.endswith('_Q')]
        wide = wide[(wide[q_columns] < threshold).any(axis=1)]
        ordered = [c for c in wide.columns if c not in columns]
        ordered = sorted(ordered, key=lambda c: (not c.startswith('COMBINED'), c))
        wide = wide[columns + ordered].sort_values(['COHORT', 'COMBINED_LONG_TAIL_Q' if 'COMBINED_LONG_TAIL_Q' in wide
                                                    else 'SET'])
    else:
        wide = pd.DataFrame(columns=columns)
    io.write_table(wide, pathways_output)

    pairs = pd.concat(pairs) if pairs else pd.DataFrame(columns=cooccurrence.COLUMNS + ['COHORT'])
    significant = pairs[pairs['Q_VALUE'] < threshold]
    io.write_table(significant[['COHORT'] + cooccurrence.COLUMNS], cooccurrence_output)

    rows = []
    for (cohort, module), group in significant.groupby(['COHORT', 'MODULE']):
        events = sorted(set(group['EVENT_1']) | set(group['EVENT_2']))
        rows.append([cohort, module, len(events), len(group), ';'.join(events)])
    io.write_table(pd.DataFrame(rows, columns=['COHORT', 'MODULE', 'N_EVENTS', 'N_PAIRS', 'EVENTS']), modules_output)


@click.command()
@click.option('--genemuts', type=click.Path(exists=True), required=True, help='dNdScv genemuts of the cohort')
@click.option('--mutations', type=click.Path(exists=True), required=True, help='Processed VEP output (parse-vep)')
@click.option('--vet', type=click.Path(exists=True), required=True, help='Vetting of the drivers (drivers-discovery)')
@click.option('--gene-sets', type=click.Path(exists=True), default=None,
              help='Gene sets (TSV or GMT). Default: $INTOGEN_DATASETS/pathways/gene_sets.tsv.gz')
@click.option('-o', '--output', type=click.Path(), required=True, help='Gene set tests')
@click.option('--cooccurrence', 'cooccurrence_output', type=click.Path(), required=True)
@click.option('--genes', 'genes_output', type=click.Path(), required=True,
              help='Significant gene sets and co-occurrence modules of each gene')
@click.option('--methylation', type=click.Path(exists=True), default=None)
@click.option('--methylation-events', type=click.Path(exists=True), default=None)
@click.option('--methylation-matrix', type=click.Path(exists=True), default=None)
@click.option('--expression', type=click.Path(exists=True), default=None)
@click.option('--expression-events', type=click.Path(exists=True), default=None)
@click.option('--expression-matrix', type=click.Path(exists=True), default=None)
@click.option('--min-size', type=int, default=10, show_default=True, help='Minimum genes of a gene set')
@click.option('--max-size', type=int, default=500, show_default=True, help='Maximum genes of a gene set')
@click.option('--threshold', type=float, default=0.1, show_default=True, help='Q-value threshold')
@click.option('--min-tumours', type=int, default=3, show_default=True,
              help='Minimum tumours with an event to test its co-occurrence')
@click.option('--max-events', type=int, default=100, show_default=True,
              help='Maximum events tested for co-occurrence (the most significant)')
def cli(gene_sets, **kwargs):
    methylation = [kwargs[k] for k in ['methylation', 'methylation_events', 'methylation_matrix']]
    expression = [kwargs[k] for k in ['expression', 'expression_events', 'expression_matrix']]
    for layer in [methylation, expression]:
        if any(layer) and not all(layer):
            raise click.UsageError('Omics layers need the results, events and matrix files')
    run(gene_sets=gene_sets or io.default_dataset('pathways', 'gene_sets.tsv.gz'),
        symbol_map=io.default_dataset('others', 'mapping_new_hugo_symbols.json'), **kwargs)


@click.command()
@click.option('--pathways', 'pathways_output', type=click.Path(), required=True)
@click.option('--cooccurrence', 'cooccurrence_output', type=click.Path(), required=True)
@click.option('--modules', 'modules_output', type=click.Path(), required=True)
@click.option('--threshold', type=float, default=0.1, show_default=True)
@click.argument('files', nargs=-1)
def summary_cli(files, **kwargs):
    summary(files, **kwargs)
