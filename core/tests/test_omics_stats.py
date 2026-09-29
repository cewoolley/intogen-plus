import itertools

import numpy as np
import pytest
from scipy.stats import binom

from intogen_core.omics.stats import fdr_bh, poisson_binomial_tail, recurrence_test, robust_scale, \
    sample_background


def brute_force_tail(p):
    pmf = np.zeros(len(p) + 1)
    for bits in itertools.product([0, 1], repeat=len(p)):
        b = np.array(bits)
        pmf[b.sum()] += np.prod(np.where(b, p, 1 - p))
    return np.cumsum(pmf[::-1])[::-1]


def test_poisson_binomial_matches_brute_force():
    p = np.random.default_rng(0).uniform(0, 0.5, 11)
    assert np.allclose(poisson_binomial_tail(p), brute_force_tail(p), atol=1e-14)


def test_poisson_binomial_reduces_to_binomial_including_far_tail():
    tail = poisson_binomial_tail([0.05] * 300)
    for k in [0, 1, 15, 40, 80]:
        assert tail[k] == pytest.approx(binom.sf(k - 1, 300, 0.05), rel=1e-9)
    assert tail[0] == 1.0


def test_poisson_binomial_degenerate_probabilities():
    tail = poisson_binomial_tail([0, 1, 1])
    assert list(tail) == [1, 1, 1, 0]


def test_fdr_bh():
    p = np.array([0.01, 0.04, np.nan, 0.03, 0.5])
    q = fdr_bh(p)
    # manual BH on the 4 finite values: sorted 0.01, 0.03, 0.04, 0.5
    expected = {0: 0.04, 1: 0.04 * 4 / 3, 3: 0.04 * 4 / 3, 4: 0.5}
    for i, v in expected.items():
        assert q[i] == pytest.approx(v)
    assert np.isnan(q[2])
    assert np.all(np.isnan(fdr_bh([np.nan, np.nan])))


def test_recurrence_test_uses_sample_background():
    # 20 genes x 10 samples; sample 0 has events in every gene (e.g. methylator phenotype)
    events = np.zeros((20, 10), dtype=bool)
    events[:, 0] = True
    events[0, 1:6] = True      # gene 0 recurrent in 5 other samples
    observed = np.ones_like(events)
    observed[1, 9] = False     # a gene with a missing value
    background = sample_background(events, observed)
    assert background[0] == 1.0
    counts, expected, pvalues = recurrence_test(events, observed, background)
    assert counts[0] == 6 and counts[5] == 1
    # the event in the methylator sample is expected: no evidence for genes 1..19
    assert pvalues[5] == pytest.approx(1.0)
    assert pvalues[0] < 1e-4
    assert expected[1] == pytest.approx(background[:9].sum())


def test_recurrence_test_only_selected_genes():
    events = np.eye(4, dtype=bool)
    observed = np.ones_like(events)
    genes = np.array([True, False, True, True])
    counts, expected, pvalues = recurrence_test(events, observed, np.full(4, 0.25), genes=genes)
    assert np.isnan(pvalues[1]) and np.isfinite(pvalues[0])


def test_robust_scale_falls_back_to_mean_absolute_deviation():
    x = np.array([[0, 0, 0, 0, 0, 10.0],     # MAD = 0
                  [1, 2, 3, 4, 5, 6.0]])
    median, scale = robust_scale(x)
    assert median[0] == 0
    assert scale[0] == pytest.approx(1.2533 * 10 / 6)
    assert scale[1] == pytest.approx(1.4826 * 1.5)
