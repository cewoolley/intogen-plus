"""
Co-occurrence of dysregulation events in the same tumours.

Events are alterations of a gene or of a gene set in a layer:

- ``mutation``: non-synonymous mutations
- ``silencing``: promoter hypermethylation with loss of expression
- ``expression_over`` / ``expression_under``: expression outliers

Tumours with many alterations (e.g. hypermutated or methylator tumours) have
many events, which makes any pair of events co-occur in them. As in DISCOVER
(Canisius et al., Genome Biology 2016), the probability that gene ``g`` is
altered in tumour ``s`` is estimated with a model that preserves both the
alteration frequency of each gene and the alteration burden of each tumour::

    p_gs = 1 / (1 + exp(-(a_g + b_s)))

fitted independently for each layer. The probability that an event (a set of
genes) happens in a tumour is ``1 - prod(1 - p_gs)``. Under independence, the
number of tumours with both events follows a Poisson-binomial distribution
with parameters ``q1_s * q2_s``, which gives the co-occurrence p-value.

Within a layer, the genes shared by two events are removed from both, so that
overlapping gene sets do not co-occur trivially. The same applies to expression
events and events of other layers, since expression changes can be a direct
consequence of the silencing or mutation of the same gene. Mutations and
silencing of the same gene (two hits) are tested.

Pairs of events that co-occur significantly are grouped into modules
(connected components), i.e. clusters of dysregulation that happen together.
"""

import itertools
import warnings

import numpy as np
import pandas as pd
from scipy.special import expit, logit

from intogen_core.omics.stats import fdr_bh, poisson_binomial_tail


COLUMNS = ['EVENT_1', 'NAME_1', 'EVENT_2', 'NAME_2', 'TUMOURS', 'TUMOURS_1', 'TUMOURS_2', 'TUMOURS_BOTH',
           'EXPECTED_BOTH', 'RATIO', 'P_VALUE', 'Q_VALUE', 'SHARED_GENES_REMOVED', 'MODULE', 'TESTABLE']


def fit_background(hits, max_iter=500, tol=1e-6):
    """
    Probability of alteration of each gene in each tumour, preserving the
    number of altered tumours of each gene and of altered genes of each tumour.

    Args:
        hits: bool array genes x tumours

    Returns:
        array of probabilities (genes x tumours)
    """
    hits = np.asarray(hits, dtype=bool)
    n_genes, n_tumours = hits.shape
    rows, cols = hits.sum(axis=1), hits.sum(axis=0)
    p = np.zeros(hits.shape)
    p[rows == n_tumours, :] = 1.0
    p[:, cols == n_genes] = 1.0
    fit_rows = (rows > 0) & (rows < n_tumours)
    fit_cols = (cols > 0) & (cols < n_genes)
    if not fit_rows.any() or not fit_cols.any():
        return p

    sub = hits[np.ix_(fit_rows, fit_cols)].astype(float)
    r, c = sub.sum(axis=1), sub.sum(axis=0)
    a = logit(np.clip(r / sub.shape[1], 1e-6, 1 - 1e-6))
    b = np.zeros(sub.shape[1])
    for _ in range(max_iter):
        q = expit(a[:, None] + b[None, :])
        w = q * (1 - q)
        a += np.clip((r - q.sum(axis=1)) / np.maximum(w.sum(axis=1), 1e-12), -5, 5)
        q = expit(a[:, None] + b[None, :])
        w = q * (1 - q)
        b += np.clip((c - q.sum(axis=0)) / np.maximum(w.sum(axis=0), 1e-12), -5, 5)
        q = expit(a[:, None] + b[None, :])
        if max(np.abs(q.sum(axis=1) - r).max(), np.abs(q.sum(axis=0) - c).max()) < tol:
            break
    p[np.ix_(fit_rows, fit_cols)] = q
    return p


class Layer:
    """Alterations of the genes of a layer in its tumours"""

    def __init__(self, name, pairs, tumours):
        """
        Args:
            name: layer name
            pairs: DataFrame with the altered SYMBOL and SAMPLE
            tumours: tumours analysed in the layer
        """
        self.name = name
        self.tumours = list(tumours)
        tumour_index = {t: i for i, t in enumerate(self.tumours)}
        pairs = pairs[pairs['SAMPLE'].isin(tumour_index)].drop_duplicates()
        self.genes = sorted(pairs['SYMBOL'].unique())
        self.gene_index = {g: i for i, g in enumerate(self.genes)}
        self.hits = np.zeros((len(self.genes), len(self.tumours)), dtype=bool)
        self.hits[pairs['SYMBOL'].map(self.gene_index).values, pairs['SAMPLE'].map(tumour_index).values] = True
        self._log_absent = None

    @property
    def log_absent(self):
        """log(1 - p) of the background model (fitted when first needed)"""
        if self._log_absent is None:
            self.p = fit_background(self.hits)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', category=RuntimeWarning)
                self._log_absent = np.log1p(-self.p)
        return self._log_absent

    def altered_tumours(self, genes):
        """Number of tumours with an alteration in any of the genes"""
        idx = [self.gene_index[g] for g in genes if g in self.gene_index]
        return int(self.hits[idx].any(axis=0).sum()) if idx else 0

    def event(self, genes):
        """(hits, probabilities) of an event made of some genes, per tumour"""
        idx = [self.gene_index[g] for g in genes if g in self.gene_index]
        if len(idx) == 0:
            return np.zeros(len(self.tumours), dtype=bool), np.zeros(len(self.tumours))
        hits = self.hits[idx].any(axis=0)
        prob = -np.expm1(self.log_absent[idx].sum(axis=0))
        return hits, prob


class Event:

    def __init__(self, layer, kind, identifier, name, genes, significance):
        self.layer = layer
        self.kind = kind                   # gene or pathway
        self.identifier = f'{layer}:{identifier}'
        self.name = name
        self.genes = frozenset(genes)
        self.significance = significance


def overlap_removed(layer1, layer2):
    """
    Whether the genes shared by two events must be removed before testing them.
    Within a layer, shared genes make events co-occur trivially. The expression of
    a gene is a consequence of its other alterations (e.g. silencing, truncating
    mutations), so expression events do not count the genes of the other event.
    A mutation and the silencing of the same gene (two hits) are kept.
    """
    return layer1 == layer2 or layer1.startswith('expression') or layer2.startswith('expression')


def test_pairs(events, layers, min_tumours=3):
    """
    Co-occurrence test of every pair of events.

    Args:
        events: list of Event
        layers: dict layer name -> Layer

    Returns:
        DataFrame (COLUMNS)
    """
    cache = {e.identifier: layers[e.layer].event(e.genes) for e in events}
    common = {}
    rows = []
    for e1, e2 in itertools.combinations(events, 2):
        shared = 0
        v1, v2 = cache[e1.identifier], cache[e2.identifier]
        if overlap_removed(e1.layer, e2.layer):
            shared = len(e1.genes & e2.genes)
            if shared > 0:
                g1, g2 = e1.genes - e2.genes, e2.genes - e1.genes
                if len(g1) == 0 or len(g2) == 0:
                    continue
                v1, v2 = layers[e1.layer].event(g1), layers[e2.layer].event(g2)

        key = (e1.layer, e2.layer)
        if key not in common:
            t2 = {t: i for i, t in enumerate(layers[e2.layer].tumours)}
            pairs = [(i, t2[t]) for i, t in enumerate(layers[e1.layer].tumours) if t in t2]
            common[key] = (np.array([p[0] for p in pairs], dtype=int), np.array([p[1] for p in pairs], dtype=int))
        i1, i2 = common[key]
        if len(i1) == 0:
            continue
        h1, q1 = v1[0][i1], v1[1][i1]
        h2, q2 = v2[0][i2], v2[1][i2]
        if h1.sum() < min_tumours or h2.sum() < min_tumours:
            continue
        both = int((h1 & h2).sum())
        probs = q1 * q2
        expected = float(probs.sum())
        tail = poisson_binomial_tail(probs)
        p_value = float(tail[both])
        # smallest p-value achievable given how often each event happens
        p_min = float(tail[min(int(h1.sum()), int(h2.sum()))])
        rows.append([e1.identifier, e1.name, e2.identifier, e2.name, len(i1), int(h1.sum()), int(h2.sum()),
                     both, expected, both / expected if expected > 0 else np.nan, p_value, np.nan, shared, None,
                     None, p_min])

    return pd.DataFrame(rows, columns=COLUMNS + ['P_MIN'])


def tarone_bh(p_values, p_min, alpha):
    """
    Benjamini-Hochberg correction restricted to the testable hypotheses
    (Tarone 1990; Gilbert 2005): discrete tests whose smallest achievable p-value
    cannot reach significance do not count in the correction, which gives power
    to detect co-occurrence of rare events. Untestable hypotheses get q = 1.
    """
    p_values, p_min = np.asarray(p_values, dtype=float), np.asarray(p_min, dtype=float)
    q = np.ones(len(p_values))
    if len(p_values) == 0:
        return q, np.zeros(0, dtype=bool)
    # smallest K such that the number of hypotheses with p_min <= alpha / K is at most K
    k = 1
    while (p_min <= alpha / k).sum() > k:
        k += 1
    testable = p_min <= alpha / k
    if testable.any():
        q[testable] = fdr_bh(p_values[testable])
    return q, testable


def modules(pairs, threshold=0.1):
    """Connected components of the significantly co-occurring events"""
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    significant = pairs[pairs['Q_VALUE'] < threshold]
    for a, b in zip(significant['EVENT_1'], significant['EVENT_2']):
        parent[find(a)] = find(b)
    groups = {}
    for event in parent:
        groups.setdefault(find(event), set()).add(event)
    ordered = sorted(groups.values(), key=lambda g: (-len(g), sorted(g)))
    return {event: f'M{i + 1}' for i, group in enumerate(ordered) for event in group}


def select_events(candidates, layers, max_events=100, min_tumours=3):
    """Most significant events altered in at least min_tumours tumours"""
    selected = []
    for event in sorted(candidates, key=lambda e: (e.significance, e.identifier)):
        if layers[event.layer].altered_tumours(event.genes) >= min_tumours:
            selected.append(event)
        if len(selected) >= max_events:
            break
    return selected


def run(candidates, layers, max_events=100, min_tumours=3, threshold=0.1):
    events = select_events(candidates, layers, max_events=max_events, min_tumours=min_tumours)
    pairs = test_pairs(events, layers, min_tumours=min_tumours)
    pairs['Q_VALUE'], pairs['TESTABLE'] = tarone_bh(pairs['P_VALUE'].values, pairs['P_MIN'].values, threshold)
    pairs = pairs[COLUMNS]
    membership = modules(pairs, threshold=threshold)
    significant = pairs['Q_VALUE'] < threshold
    pairs.loc[significant, 'MODULE'] = pairs.loc[significant, 'EVENT_1'].map(membership)
    return pairs.sort_values(['P_VALUE', 'EVENT_1', 'EVENT_2']), events
