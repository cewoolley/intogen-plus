import json

import numpy as np
import pandas as pd
import pytest
from scipy.special import expit

import synthetic_pathways
from intogen_core.omics.stats import (binomial_pmf, fisher_combine, poisson_binomial_tail, simes_combine, sum_pmfs,
                                     upper_tail)
from intogen_core.pathways import analysis, cooccurrence
from intogen_core.pathways.genesets import GeneSets, read_gene_sets
from intogen_core.pathways.selection import MUTATION_LAYERS, EventLayer, MutationLayer, estimate_theta


def neutral_genemuts(rng, n=600, theta=8.0, covariates=True):
    length = rng.lognormal(7.2, 0.5, n)
    pred = rng.lognormal(0, 0.3, n)
    rate = rng.gamma(theta, pred / theta)
    df = pd.DataFrame({'gene_name': [f'G{i}' for i in range(n)]})
    for k, f in {'syn': 0.25, 'mis': 0.65, 'non': 0.06, 'spl': 0.04}.items():
        df[f'n_{k}'] = rng.poisson(length * 2.5e-3 * f * rate)
        df[f'exp_{k}'] = length * 2.5e-3 * f
    if covariates:
        df['exp_syn_cv'] = df['exp_syn'] * pred
    return df


def test_sum_of_binomials_is_exact():
    rng = np.random.default_rng(0)
    n, p = rng.integers(0, 6, 30), rng.uniform(0.3, 0.9, 30)
    pmf = sum_pmfs([binomial_pmf(k, q) for k, q in zip(n, p)])
    tail = poisson_binomial_tail(np.repeat(p, n))
    assert all(abs(upper_tail(pmf, x) - tail[x]) < 1e-12 for x in range(len(tail)))
    assert upper_tail(pmf, 0) == 1.0 and upper_tail(pmf, len(pmf)) == 0.0
    assert fisher_combine([np.nan]) != fisher_combine([np.nan])   # NaN
    assert fisher_combine([0.01, 0.01]) < 0.01
    assert simes_combine([0.01, 0.5, np.nan]) == pytest.approx(0.02)
    assert simes_combine([0.04, 0.05, 0.9]) == pytest.approx(0.075)
    assert np.isnan(simes_combine([np.nan]))


def test_theta_estimate():
    rng = np.random.default_rng(1)
    mean = rng.lognormal(0.5, 0.5, 20000)
    counts = rng.poisson(rng.gamma(4.0, mean / 4.0))
    assert estimate_theta(counts, mean) == pytest.approx(4.0, rel=0.25)
    assert np.isinf(estimate_theta(rng.poisson(mean), mean)) or estimate_theta(rng.poisson(mean), mean) > 50


def test_theta_lower_bound():
    """The lower bound used by the tests covers the true overdispersion with few genes"""
    covered, ratio = [], []
    for rep in range(60):
        rng = np.random.default_rng(100 + rep)
        mean = rng.lognormal(0.3, 0.5, 800)
        counts = rng.poisson(rng.gamma(8.0, mean / 8.0))
        bound, mle = estimate_theta(counts, mean, confidence=0.95), estimate_theta(counts, mean)
        covered.append(bound <= 8.0)
        ratio.append(bound / mle)
    assert np.mean(covered) >= 0.9
    assert max(ratio) < 1
    # the bound is close to the estimate with a whole exome
    rng = np.random.default_rng(7)
    mean = rng.lognormal(0.3, 0.5, 19000)
    counts = rng.poisson(rng.gamma(8.0, mean / 8.0))
    assert estimate_theta(counts, mean, confidence=0.95) / estimate_theta(counts, mean) > 0.75


@pytest.mark.parametrize('method', ['nb', 'conditional'])
def test_mutation_set_test_is_calibrated(method):
    pvalues = []
    for rep in range(250):
        rng = np.random.default_rng(500 + rep)
        df = neutral_genemuts(rng)
        genes = rng.choice(df['gene_name'], rng.integers(15, 40), replace=False)
        pvalues.append(MutationLayer(df, 'mutation', method=method).test(genes)['P_VALUE'])
    pvalues = np.array(pvalues)
    assert (pvalues < 0.05).mean() <= 0.08
    assert (pvalues < 0.01).mean() <= 0.025


def test_mutation_layer_methods():
    rng = np.random.default_rng(2)
    assert MutationLayer(neutral_genemuts(rng), 'mutation').method == 'nb'
    layer = MutationLayer(neutral_genemuts(rng, covariates=False), 'mutation')
    assert layer.method == 'conditional' and layer.theta is None
    r = layer.test(['G1', 'G2', 'NOT_A_GENE'])
    assert r['N_GENES'] == 2


def test_event_layer():
    layer = EventLayer('silencing', pd.Series([5, 0, 3], index=['A', 'B', 'C']), pd.Series([1.0, 1.0, 1.0],
                                                                                       index=['A', 'B', 'C']))
    r = layer.test(['A', 'B', 'C', 'X'])
    assert r['OBSERVED'] == 8 and r['EXPECTED'] == 3 and r['N_GENES'] == 3 and r['N_GENES_ALTERED'] == 2
    assert r['P_VALUE'] == pytest.approx(1 - 0.988095, abs=1e-5)   # Poisson(3) P(X >= 8)
    assert np.isnan(layer.test(['X'])['P_VALUE'])


def test_gene_sets(tmp_path):
    gmt = tmp_path / 'sets.gmt'
    gmt.write_text('SMALL\tdesc\tA\tB\n'
                   'BIG\tdesc\t' + '\t'.join(f'G{i}' for i in range(20)) + '\tOLD\tNOPE\n')
    table = read_gene_sets(str(gmt))
    assert set(table['SOURCE']) == {'custom'}
    universe = {f'G{i}' for i in range(20)} | {'A', 'B', 'NEW'}
    sets = GeneSets(table, universe, min_size=3, max_size=25, symbol_map={'OLD': 'NEW'})
    assert sets.ids == ['BIG'] and sets.discarded == 1
    assert 'NEW' in sets.genes['BIG'] and 'NOPE' not in sets.genes['BIG']


def test_background_model_keeps_margins():
    rng = np.random.default_rng(3)
    hits = rng.uniform(size=(200, 80)) < expit(rng.normal(-3, 1, 200)[:, None] + rng.normal(0, 1, 80)[None, :])
    hits[0] = True           # always altered
    hits[1] = False          # never altered
    p = cooccurrence.fit_background(hits)
    assert np.allclose(p.sum(axis=1), hits.sum(axis=1), atol=1e-4)
    assert np.allclose(p.sum(axis=0), hits.sum(axis=0), atol=1e-4)
    assert (p[0] == 1).all() and (p[1] == 0).all()


def test_tarone_bh():
    p = np.array([0.001, 0.01, 0.5, 0.2])
    p_min = np.array([1e-6, 1e-4, 0.3, 0.2])    # the last two can never be significant
    q, testable = cooccurrence.tarone_bh(p, p_min, 0.1)
    assert list(testable) == [True, True, False, False]
    assert q[0] == pytest.approx(0.002) and q[1] == pytest.approx(0.01) and q[2] == 1


def test_overlapping_events():
    tumours = [f'T{i}' for i in range(20)]
    # A is silenced, under-expressed and mutated in the same 10 tumours
    pairs = pd.DataFrame({'SYMBOL': ['A'] * 10 + ['B'] * 3 + ['C'] * 3, 'SAMPLE': tumours[:10] + tumours[:3] + tumours[10:13]})
    layers = {name: cooccurrence.Layer(name, pairs, tumours) for name in ['silencing', 'expression_under', 'mutation']}
    events = [cooccurrence.Event(layer, 'gene', 'A', 'A', ['A'], 0.0) for layer in layers]
    tested = cooccurrence.test_pairs(events, layers)
    tested_pairs = set(zip(tested['EVENT_1'], tested['EVENT_2']))
    # expression changes of a gene are not co-occurrence with its own silencing or mutation
    assert tested_pairs == {('silencing:A', 'mutation:A')}
    assert cooccurrence.overlap_removed('mutation', 'mutation')
    # truncating mutations are also mutations
    assert cooccurrence.overlap_removed('mutation', 'mutation_truncating')
    assert not cooccurrence.overlap_removed('mutation', 'silencing')


@pytest.fixture(scope='module')
def pathway_results(tmp_path_factory):
    folder = str(tmp_path_factory.mktemp('pathways'))
    files = synthetic_pathways.PathwayCohort(seed=2).write(folder)
    out = {k: f'{folder}/C.{k}.tsv.gz' for k in ['pathways', 'pathway_cooccurrence', 'pathway_genes']}
    analysis.run(output=out['pathways'], cooccurrence_output=out['pathway_cooccurrence'],
                 genes_output=out['pathway_genes'], **files)
    out['folder'] = folder
    return out


def test_pathway_selection(pathway_results):
    raw = pd.read_csv(pathway_results['pathways'], sep='\t')
    assert list(raw.columns) == analysis.RESULT_COLUMNS
    df = raw.set_index(['SET', 'SCOPE', 'LAYER'])
    # selected through genes that are not individually significant
    assert df.loc[('LONGTAIL', 'long_tail', 'mutation'), 'Q_VALUE'] < 1e-3
    assert df.loc[('LONGTAIL', 'long_tail', 'mutation'), 'RATIO'] > 1.5
    # significant only because of its driver
    assert df.loc[('DRIVER_SET', 'all', 'mutation'), 'Q_VALUE'] < 0.05
    assert df.loc[('DRIVER_SET', 'long_tail', 'mutation'), 'Q_VALUE'] > 0.1
    assert df.loc[('DRIVER_SET', 'long_tail', 'mutation'), 'N_GENES_EXCLUDED'] == 1
    assert df.loc[('SILENCED', 'long_tail', 'silencing'), 'Q_VALUE'] < 1e-3
    assert df.loc[('SILENCED', 'long_tail', 'combined'), 'Q_VALUE'] < 1e-3
    neutral = df.reset_index()
    neutral = neutral[neutral['SET'].str.startswith('NEUTRAL') & (neutral['LAYER'] == 'combined')]
    assert (neutral['Q_VALUE'] < 0.01).sum() == 0
    stats = json.load(open(pathway_results['pathways'] + '.stats.json'))
    assert stats['mutation_method'] == 'nb' and stats['gene_sets'] == 23


def test_pathway_cooccurrence(pathway_results):
    pairs = pd.read_csv(pathway_results['pathway_cooccurrence'], sep='\t')
    assert list(pairs.columns) == cooccurrence.COLUMNS
    hit = pairs[pairs['EVENT_1'].isin(['mutation:DRV1', 'silencing:SILENCED']) &
                pairs['EVENT_2'].isin(['mutation:DRV1', 'silencing:SILENCED'])].iloc[0]
    assert hit['Q_VALUE'] < 0.1 and hit['MODULE'] == 'M1'
    assert hit['TUMOURS_BOTH'] > hit['EXPECTED_BOTH']


def test_pathway_genes_and_summary(pathway_results, tmp_path):
    genes = pd.read_csv(pathway_results['pathway_genes'], sep='\t')
    assert list(genes.columns) == ['SYMBOL', 'PATHWAYS', 'MODULES']
    genes = genes.set_index('SYMBOL')
    assert genes.loc['G5', 'PATHWAYS'] == 'LONGTAIL'
    # DRIVER_SET is not significant beyond DRV1, but DRV1 co-occurs with the silencing of SILENCED
    assert pd.isna(genes.loc['DRV1', 'PATHWAYS']) and genes.loc['DRV1', 'MODULES'] == 'M1'
    assert genes['MODULES'].notnull().sum() == 1     # pathway events are not assigned to their genes

    files = [pathway_results['pathways'], pathway_results['pathway_cooccurrence']]
    out = {k: str(tmp_path / f'{k}.tsv') for k in ['pathways', 'cooccurrence', 'modules']}
    analysis.summary(files, out['pathways'], out['cooccurrence'], out['modules'])
    summary = pd.read_csv(out['pathways'], sep='\t')
    assert list(summary.columns[:5]) == ['COHORT', 'SET', 'SOURCE', 'NAME', 'N_GENES']
    assert {'COMBINED_Q', 'COMBINED_LONG_TAIL_Q', 'MUTATION_LONG_TAIL_RATIO', 'SILENCING_LONG_TAIL_Q',
            'MUTATION_LONG_TAIL_TOP_GENES'} <= set(summary.columns)
    assert {'LONGTAIL', 'SILENCED'} <= set(summary['SET'])
    modules = pd.read_csv(out['modules'], sep='\t')
    assert modules.iloc[0]['EVENTS'] == 'mutation:DRV1;silencing:SILENCED'
    assert set(pd.read_csv(out['cooccurrence'], sep='\t')['COHORT']) == {'C'}


def test_empty_inputs(tmp_path):
    """Cohorts without drivers or dNdScv results"""
    files = synthetic_pathways.PathwayCohort(seed=4).write(str(tmp_path))
    pd.DataFrame(columns=['SYMBOL', 'TIER', 'FILTER']).to_csv(files['vet'], sep='\t', index=False)
    pd.DataFrame(columns=['gene_name', 'n_syn', 'n_mis', 'n_non', 'n_spl', 'exp_syn', 'exp_mis', 'exp_non',
                          'exp_spl', 'exp_syn_cv']).to_csv(files['genemuts'], sep='\t', index=False)
    out = str(tmp_path / 'x.pathways.tsv.gz')
    analysis.run(output=out, cooccurrence_output=str(tmp_path / 'x.co.tsv.gz'),
                 genes_output=str(tmp_path / 'x.genes.tsv.gz'), **files)
    df = pd.read_csv(out, sep='\t')
    assert set(df['LAYER']) == {'silencing', 'combined'}
    assert 'warning_mutations' in json.load(open(out + '.stats.json'))


@pytest.fixture(scope='module')
def sequencing_results(tmp_path_factory):
    """Mutations only: the selection network is DRV2 and the long tails of LOF_TAIL and MIS_TAIL"""
    folder = str(tmp_path_factory.mktemp('sequencing'))
    files = synthetic_pathways.SequencingCohort(seed=1).write(folder)
    out = {k: f'{folder}/S.{k}.tsv.gz' for k in ['pathways', 'pathway_cooccurrence', 'pathway_genes']}
    analysis.run(output=out['pathways'], cooccurrence_output=out['pathway_cooccurrence'],
                 genes_output=out['pathway_genes'], **files)
    return out


def test_sequencing_only_selection(sequencing_results):
    df = pd.read_csv(sequencing_results['pathways'], sep='\t')
    assert set(df['LAYER']) == set(MUTATION_LAYERS) | {'combined'}
    df = df.set_index(['SET', 'SCOPE', 'LAYER'])
    # selected through truncating mutations only: diluted among all the mutations
    lof = df.loc['LOF_TAIL'].xs('long_tail')
    assert lof.loc['mutation_truncating', 'Q_VALUE'] < 1e-6
    assert lof.loc['combined', 'Q_VALUE'] < 1e-6
    assert lof.loc['combined', 'Q_VALUE'] < lof.loc['mutation', 'Q_VALUE']
    for name in ['LONGTAIL', 'MIS_TAIL']:
        assert df.loc[(name, 'long_tail', 'combined'), 'Q_VALUE'] < 0.01
    assert df.loc[('DRIVER_SET', 'long_tail', 'combined'), 'Q_VALUE'] > 0.1
    stats = json.load(open(sequencing_results['pathways'] + '.stats.json'))
    assert stats['layers'] == list(MUTATION_LAYERS)
    assert stats['mutation_theta'] < stats['mutation_theta_mle']
    genes = pd.read_csv(sequencing_results['pathway_genes'], sep='\t').set_index('SYMBOL')
    assert genes.loc['G50', 'PATHWAYS'] == 'LOF_TAIL'


def test_sequencing_only_network(sequencing_results):
    pairs = pd.read_csv(sequencing_results['pathway_cooccurrence'], sep='\t')
    significant = pairs[pairs['Q_VALUE'] < 0.1]
    modules = dict(zip(significant['EVENT_1'], significant['MODULE']))
    modules.update(zip(significant['EVENT_2'], significant['MODULE']))
    # gene sets are events with the mutations they are selected through
    assert modules['mutation:DRV2'] == modules['mutation_truncating:LOF_TAIL'] == 'M1'
    assert modules['mutation_missense:MIS_TAIL'] == 'M1'
    assert set(modules) == {'mutation:DRV2', 'mutation_truncating:LOF_TAIL', 'mutation_missense:MIS_TAIL'}
    # selected, but not together with the network
    events = set(pairs['EVENT_1']) | set(pairs['EVENT_2'])
    assert {'mutation:DRV1', 'mutation_missense:LONGTAIL'} <= events
    hit = significant[significant['EVENT_1'].isin(['mutation:DRV2', 'mutation_truncating:LOF_TAIL']) &
                      significant['EVENT_2'].isin(['mutation:DRV2', 'mutation_truncating:LOF_TAIL'])].iloc[0]
    assert hit['Q_VALUE'] < 1e-4 and hit['TUMOURS_BOTH'] > 2 * hit['EXPECTED_BOTH']
    genes = pd.read_csv(sequencing_results['pathway_genes'], sep='\t').set_index('SYMBOL')
    assert genes.loc['DRV2', 'MODULES'] == 'M1'
    assert 'DRV1' not in genes.index or pd.isna(genes.loc['DRV1', 'MODULES'])
