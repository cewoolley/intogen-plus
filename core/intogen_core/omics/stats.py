"""
Statistics shared by the omics analyses.

The recurrence test used for the omics layers mirrors the way mutation-based
methods account for the mutation burden of each sample: an omics event
(e.g. promoter hypermethylation, expression outlier) is expected to occur in
sample ``s`` with probability ``pi_s``, the fraction of the analysed genes
with an event in that sample. Under the null hypothesis the number of samples
with an event in a gene follows a Poisson-binomial distribution with
parameters ``pi_1 ... pi_n``. Samples with a genome-wide excess of events
(e.g. CpG island methylator phenotype) therefore contribute little evidence,
analogously to hypermutated samples in mutation analyses.
"""

import numpy as np
from scipy import stats as sps


def fdr_bh(pvalues):
    """
    Benjamini-Hochberg q-values. Non-finite p-values are ignored
    and get a NaN q-value.
    """
    p = np.asarray(pvalues, dtype=float)
    q = np.full(p.shape, np.nan)
    mask = np.isfinite(p)
    n = int(mask.sum())
    if n == 0:
        return q
    pv = p[mask]
    order = np.argsort(pv, kind='mergesort')
    ranked = pv[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    q[mask] = out
    return q


def poisson_binomial_tail(probs):
    """
    Upper tail of the Poisson-binomial distribution.

    Args:
        probs: success probability of each Bernoulli trial

    Returns:
        array ``tail`` of length n + 1 with ``tail[k] = P(X >= k)``
    """
    probs = np.clip(np.asarray(probs, dtype=float), 0.0, 1.0)
    n = len(probs)
    pmf = np.zeros(n + 1)
    pmf[0] = 1.0
    for i, p in enumerate(probs):
        pmf[1:i + 2] = pmf[1:i + 2] * (1.0 - p) + pmf[0:i + 1] * p
        pmf[0] *= (1.0 - p)
    tail = np.clip(np.cumsum(pmf[::-1])[::-1], 0.0, 1.0)
    tail[0] = 1.0  # exact (avoids rounding errors: p-values of 1 are ignored when combining methods)
    return tail


def sample_background(events, observed, genes=None):
    """
    Per-sample event probability: fraction of the (selected) genes with an
    event among the genes with data in the sample.

    Args:
        events: bool array genes x samples
        observed: bool array genes x samples (non-missing values)
        genes: optional bool mask of the genes to use

    Returns:
        array of probabilities (one per sample)
    """
    if genes is not None:
        events, observed = events[genes], observed[genes]
    n_obs = observed.sum(axis=0)
    n_events = (events & observed).sum(axis=0)
    return np.where(n_obs > 0, n_events / np.maximum(n_obs, 1), 0.0)


def recurrence_test(events, observed, background, genes=None, max_cache=1024):
    """
    Sample-aware recurrence test for each gene.

    Args:
        events: bool array genes x samples
        observed: bool array genes x samples
        background: per-sample event probabilities
        genes: optional bool mask of the genes to test (others get NaN)
        max_cache: maximum number of distributions kept in memory

    Returns:
        (observed events, expected events, p-values) arrays
    """
    n_genes = events.shape[0]
    counts = (events & observed).sum(axis=1)
    expected = np.full(n_genes, np.nan)
    pvalues = np.full(n_genes, np.nan)
    if genes is None:
        genes = np.ones(n_genes, dtype=bool)

    cache = {}
    for i in np.flatnonzero(genes):
        obs = observed[i]
        key = np.packbits(obs).tobytes()
        tail = cache.get(key)
        if tail is None:
            if len(cache) >= max_cache:
                cache.clear()
            tail = poisson_binomial_tail(background[obs])
            cache[key] = tail
        expected[i] = background[obs].sum()
        pvalues[i] = tail[counts[i]]
    return counts, expected, pvalues


def robust_scale(values):
    """
    Robust location and scale of each row.

    The scale is the median absolute deviation (MAD) scaled to be consistent
    with the standard deviation (x 1.4826). When the MAD is 0 (more than half of
    the values are identical, e.g. genes not expressed in most samples), the mean
    absolute deviation scaled by 1.2533 is used instead.

    Returns:
        (median, scale) arrays
    """
    median = np.nanmedian(values, axis=1)
    deviation = np.abs(values - median[:, None])
    mad = np.nanmedian(deviation, axis=1) * 1.4826
    meanad = np.nanmean(deviation, axis=1) * 1.2533
    scale = np.where(mad > 0, mad, meanad)
    return median, scale


def mannwhitney(x, y, alternative='two-sided'):
    """Mann-Whitney U test p-value (NaN when any group is empty)"""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    if len(x) == 0 or len(y) == 0:
        return np.nan
    if np.all(x == x[0]) and np.all(y == x[0]):
        return 1.0
    return float(sps.mannwhitneyu(x, y, alternative=alternative).pvalue)


def spearman(x, y):
    """Spearman correlation of paired values (NaN if not computable)"""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3 or np.all(x[mask] == x[mask][0]) or np.all(y[mask] == y[mask][0]):
        return np.nan
    return float(sps.spearmanr(x[mask], y[mask])[0])
