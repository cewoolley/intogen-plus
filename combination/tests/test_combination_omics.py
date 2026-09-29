import itertools

import numpy as np
import pandas as pd
import pytest

import synthetic_methods as synthetic
from intogen_combination import config, grid_optimizer, parser


def test_methods_configuration():
    assert config.METHODS == ['oncodriveclustl', 'dndscv', 'oncodrivefml', 'hotmaps', 'smregions', 'cbase',
                              'mutpanning']
    assert config.OPTIONAL_METHODS == ['methylation', 'expression']
    assert config.RESTRICTED_METHODS == ['methylation', 'expression']
    files = {m: 'x' for m in config.METHODS}
    assert config.active_methods(files) == config.METHODS
    files['expression'] = 'y'
    assert config.active_methods(files) == config.METHODS + ['expression']


def legacy_grid_optimize(func, methods, low_quality=None):
    """Grid search as implemented before omics were added (full grid)"""
    optimum = {k: None for k in methods}
    optimum.update({'Objective_Function': 0})
    low_quality_index = set()
    if low_quality is not None:
        low_quality_index = [i for i, m in enumerate(methods) if m in low_quality]
    dim = len(methods) - len(low_quality_index)
    for w in itertools.product(np.linspace(0, 1, 21), repeat=dim - 1):
        if sum(w) <= 1 - grid_optimizer.LOWER_BOUND:
            w_dim = list(np.append(w, [1 - sum(w)]))
            if grid_optimizer.all_constraints(w_dim):
                w_all = grid_optimizer.fill_with_zeros(w_dim, low_quality_index, methods)
                f = func(w_all)
                if optimum['Objective_Function'] > f:
                    optimum['Objective_Function'] = f
                    for i, v in enumerate(methods):
                        optimum[v] = w_all[i]
    return optimum


@pytest.mark.parametrize('low_quality', [None, {'m2'}])
def test_grid_visits_the_same_points(low_quality):
    methods = ['m0', 'm1', 'm2', 'm3', 'm4', 'm5']
    visited_new, visited_old = [], []
    rng = np.random.default_rng(0)
    coefficients = rng.normal(size=len(methods))

    def objective(visited):
        def f(w):
            visited.append(tuple(w))
            # many ties, to check that the first optimum is kept in both cases
            return -round(float(np.dot(coefficients, w)), 1)
        return f

    new = grid_optimizer.grid_optimize(objective(visited_new), low_quality=low_quality, methods=methods)
    old = legacy_grid_optimize(objective(visited_old), methods, low_quality=low_quality)
    assert visited_new == visited_old
    assert len(visited_new) > 0
    assert new == old


def test_candidate_restriction(tmp_path):
    files = synthetic.make_method_outputs(str(tmp_path), seed=1)
    candidates = parser.candidate_genes(files)
    fml = pd.read_csv(files['oncodrivefml'], sep='\t')
    assert candidates == set(fml.loc[fml['Q_VALUE'].notna(), 'SYMBOL'])

    ranking, pvalues = parser.parse(**files)
    for method in ['methylation', 'expression']:
        assert set(ranking[method]) <= candidates
        assert len(ranking[method]) > 0
    # non candidate genes do not get omics p-values (the multiple testing universe is unchanged)
    assert not any(g.startswith('NOMUT') for g in pvalues)
    genes_with_omics = {g for g, v in pvalues.items() if 'methylation' in v}
    assert genes_with_omics <= candidates

    # mutation-based methods are not restricted
    all_dndscv = {g for g, v in pvalues.items() if 'dndscv' in v}
    assert all_dndscv - candidates


def test_qc_filter_uses_the_same_restriction(tmp_path):
    files = synthetic.make_method_outputs(str(tmp_path), seed=2)
    qc = grid_optimizer.Filter(**files)
    assert set(qc.data.index) == set(files)
    from intogen_combination.qc.parser import Parser
    df = Parser('methylation', config.REGIONS, candidates=qc.candidates).read(files['methylation'])
    assert set(df['GENE_ID']) <= qc.candidates
