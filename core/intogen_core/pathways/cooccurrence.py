"""
Co-occurrence and mutual exclusivity of dysregulation events in the same tumours.

Events are alterations of a gene or of a gene set in a layer:

- ``mutation``: mutations that alter the protein (including indels)
- ``mutation_missense`` / ``mutation_truncating``: missense or truncating mutations
- ``silencing``: promoter hypermethylation with loss of expression
- ``expression_over`` / ``expression_under``: expression outliers

Tumours with many alterations (e.g. hypermutated or methylator tumours) have
many events, which makes any pair of events co-occur in them. How much the
probability of an event depends on the alteration burden of the tumour differs
between events: passenger-rich events are proportional to the burden, while
the selected alterations of drivers barely depend on it, or even decrease with it
(e.g. EGFR mutations in lung adenocarcinomas of non-smokers). Models that
impose the same dependence on all the genes, such as the additive model of
DISCOVER (Canisius et al., Genome Biology 2016), make driver co-occurrence
undetectable and passengers co-occur in tumours with heterogeneous burdens.
Each event is therefore given its own burden elasticity: the probability that
the event happens in tumour ``s`` is estimated with a logistic regression on
the burden of the tumour in the layer of the event (number of altered genes,
without the genes of the event)::

    q_s = 1 / (1 + exp(-(a + b * log(1 + burden_s))))

Under independence given the burden, the number of tumours with both events
follows a Poisson-binomial distribution with probabilities ``q1_s * q2_s``. Its
upper tail gives the co-occurrence p-value and its lower tail the mutual
exclusivity p-value.

Within a layer (the mutation layers count as one), the genes shared by two
events are removed from both, so that overlapping gene sets do not co-occur
trivially. The same applies to expression
events and events of other layers, since expression changes can be a direct
consequence of the silencing or mutation of the same gene. Mutations and
silencing of the same gene (two hits) are tested.

Pairs of events that co-occur significantly are grouped into modules
(connected components), i.e. clusters of dysregulation that happen together.
"""

import itertools

import numpy as np
import pandas as pd
from scipy.special import expit, logit

from intogen_core.omics.stats import fdr_bh, poisson_binomial_tail


COLUMNS = ['EVENT_1', 'NAME_1', 'EVENT_2', 'NAME_2', 'TUMOURS', 'TUMOURS_1', 'TUMOURS_2', 'TUMOURS_BOTH',
           'EXPECTED_BOTH', 'RATIO', 'P_VALUE', 'Q_VALUE', 'P_VALUE_EXCLUSIVITY', 'Q_VALUE_EXCLUSIVITY',
           'SHARED_GENES_REMOVED', 'MODULE', 'TESTABLE', 'TESTABLE_EXCLUSIVITY']


def event_probabilities(hits, log_burden, ridge=1e-4, max_iter=100):
    """
    Probability of an event in each tumour, from a logistic regression of its
    occurrence on the (log) alteration burden of the tumours.
    """
    h = np.asarray(hits, dtype=float)
    if h.sum() == 0 or h.sum() == len(h):
        return h
    X = np.column_stack([np.ones(len(h)), log_burden])
    beta = np.array([logit(h.mean()), 0.0])
    penalty = np.diag([0.0, ridge])
    for _ in range(max_iter):
        p = expit(X @ beta)
        step = np.linalg.solve((X * (p * (1 - p))[:, None]).T @ X + penalty, X.T @ (h - p) - penalty @ beta)
        beta += step
        if np.abs(step).max() < 1e-9:
            break
    return expit(X @ beta)


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
        self.burden = self.hits.sum(axis=0)

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
        own = self.hits[idx].sum(axis=0)
        return hits, event_probabilities(hits, np.log1p(self.burden - own))


class Event:

    def __init__(self, layer, kind, identifier, name, genes, significance):
        self.layer = layer
        self.kind = kind                   # gene or pathway
        self.identifier = f'{layer}:{identifier}'
        self.name = name
        self.genes = frozenset(genes)
        self.significance = significance


def family(layer):
    """Layers with the same alterations (all the mutation layers are mutations)"""
    return 'mutation' if layer.startswith('mutation') else layer


def overlap_removed(layer1, layer2):
    """
    Whether the genes shared by two events must be removed before testing them.
    Within a layer, shared genes make events co-occur trivially. The expression of
    a gene is a consequence of its other alterations (e.g. silencing, truncating
    mutations), so expression events do not count the genes of the other event.
    A mutation and the silencing of the same gene (two hits) are kept.
    """
    return family(layer1) == family(layer2) or layer1.startswith('expression') or layer2.startswith('expression')


def test_pairs(events, layers, min_tumours=3):
    """
    Co-occurrence and mutual exclusivity tests of every pair of events.

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
        tail = np.append(poisson_binomial_tail(probs), 0.0)     # tail[k] = P(X >= k)
        n1, n2 = int(h1.sum()), int(h2.sum())
        # smallest p-values achievable given how often each event happens
        p_min = float(tail[min(n1, n2)])
        p_min_exclusivity = float(1 - tail[max(0, n1 + n2 - len(i1)) + 1])
        rows.append([e1.identifier, e1.name, e2.identifier, e2.name, len(i1), n1, n2, both, expected,
                     both / expected if expected > 0 else np.nan, float(tail[both]), np.nan,
                     float(1 - tail[both + 1]), np.nan, shared, None, None, None, p_min, p_min_exclusivity])

    return pd.DataFrame(rows, columns=COLUMNS + ['P_MIN', 'P_MIN_EXCLUSIVITY'])


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
    pairs['Q_VALUE_EXCLUSIVITY'], pairs['TESTABLE_EXCLUSIVITY'] = tarone_bh(
        pairs['P_VALUE_EXCLUSIVITY'].values, pairs['P_MIN_EXCLUSIVITY'].values, threshold)
    pairs['MIN_P'] = pairs[['P_VALUE', 'P_VALUE_EXCLUSIVITY']].min(axis=1)
    pairs = pairs.sort_values(['MIN_P', 'EVENT_1', 'EVENT_2'])
    pairs = pairs[COLUMNS]
    membership = modules(pairs, threshold=threshold)
    significant = pairs['Q_VALUE'] < threshold
    pairs.loc[significant, 'MODULE'] = pairs.loc[significant, 'EVENT_1'].map(membership)
    return pairs, events
