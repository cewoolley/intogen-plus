"""
Set-level tests of the layers of alteration.

Mutations (``mutation``, ``mutation_missense``, ``mutation_truncating``)
    Uses the observed and expected (neutral) numbers of synonymous and
    non-synonymous substitutions per gene computed by dNdScv (genemuts).

    Default (``nb``), following the model of dNdScv: the local mutation rate of
    each gene is predicted from genomic covariates (``exp_syn_cv``) with
    gamma-distributed variation (shape ``theta``, estimated from the synonymous
    mutations of all genes; the lower bound of its 95% confidence interval is
    used, to account for its uncertainty). Updated with the synonymous mutations observed in
    the gene, the number of non-synonymous mutations expected under neutrality
    follows a negative binomial distribution. The number of non-synonymous
    mutations of a set is compared with the sum of these distributions
    (exact convolution; one-sided test).

    Fallback (``conditional``, without covariates): conditional on the number
    of mutations of a gene, its number of non-synonymous mutations is binomial
    with probability ``e_ns / (e_ns + e_syn)`` under neutrality, so that the
    local mutation rate cancels out. More robust but less powerful.

    The ratio reported is the observed/expected number of non-synonymous
    mutations (dN/dS of the set).

Omics layers (``silencing``, ``expression_over``, ``expression_under``)
    The number of events (e.g. hypermethylated tumours) of the genes of a set
    is compared with the number expected given the event rate of each tumour
    (the sum of the per-gene expectations of the omics analyses). Since the
    events are many independent rare events, their sum is approximated by a
    Poisson distribution, which is conservative.
"""

import numpy as np
from scipy import optimize
from scipy import stats as sps

from intogen_core.omics.stats import binomial_pmf, sum_pmfs, upper_tail


MUTATION_LAYERS = {
    'mutation': (['n_mis', 'n_non', 'n_spl'], ['exp_mis', 'exp_non', 'exp_spl']),
    'mutation_missense': (['n_mis'], ['exp_mis']),
    'mutation_truncating': (['n_non', 'n_spl'], ['exp_non', 'exp_spl']),
}

OMICS_LAYERS = ['silencing', 'expression_over', 'expression_under']


def estimate_theta(n_syn, mean, confidence=None):
    """
    Maximum likelihood shape (overdispersion) of a negative binomial model
    of the synonymous mutations per gene given their covariate-based expectation.
    Returns inf when there is no overdispersion.

    With a confidence level, returns instead the lower bound of its one-sided
    profile likelihood confidence interval. Smaller shapes mean more
    variation of the mutation rate between genes, so that set-level tests are
    conservative with respect to the uncertainty of the estimate, which is large
    when few genes are analysed.
    """
    ok = np.isfinite(mean) & (mean > 0)
    n_syn, mean = n_syn[ok], mean[ok]
    if len(n_syn) < 10 or n_syn.sum() == 0:
        return np.inf

    def nll(log_theta):
        theta = np.exp(log_theta)
        return -sps.nbinom.logpmf(n_syn, theta, theta / (theta + mean)).sum()

    lower, upper = -6, 12
    res = optimize.minimize_scalar(nll, bounds=(lower, upper), method='bounded')
    if confidence is None:
        return np.inf if res.x > upper - 0.5 else float(np.exp(res.x))

    cutoff = res.fun + sps.chi2.ppf(2 * confidence - 1, 1) / 2
    if nll(lower) <= cutoff:
        return float(np.exp(lower))
    return float(np.exp(optimize.brentq(lambda x: nll(x) - cutoff, lower, res.x)))


class MutationLayer:
    """
    Selection on the substitutions of the genes of a set (from dNdScv genemuts).

    Args:
        genemuts: dNdScv genemuts table
        layer: mutation, mutation_missense or mutation_truncating
        method: nb (covariates, default when exp_syn_cv is available) or conditional
        theta: overdispersion of the nb method (default: lower bound of its 95% confidence interval)
    """

    def __init__(self, genemuts, layer, method=None, theta=None):
        self.layer = layer
        observed, expected = MUTATION_LAYERS[layer]
        if method is None:
            method = 'nb' if 'exp_syn_cv' in genemuts.columns and genemuts['exp_syn_cv'].notna().any() else 'conditional'
        self.method = method
        x = genemuts[observed].sum(axis=1).values.astype(int)
        s = genemuts['n_syn'].values.astype(int)
        e_ns = genemuts[expected].sum(axis=1).values.astype(float)
        e_syn = genemuts['exp_syn'].values.astype(float)
        mu = genemuts['exp_syn_cv'].values.astype(float) if method == 'nb' else e_syn
        valid = np.isfinite(e_ns) & np.isfinite(e_syn) & (e_ns + e_syn > 0)
        if method == 'nb':
            valid &= np.isfinite(mu) & (mu > 0) & (e_syn > 0)
            self.theta = estimate_theta(s[valid], mu[valid], confidence=0.95) if theta is None else theta
        else:
            self.theta = None
        self.data = {}
        for g, xi, si, ei, esi, mui, ok in zip(genemuts['gene_name'], x, s, e_ns, e_syn, mu, valid):
            if ok:
                self.data[g] = (xi, si, ei, esi, mui)
        self.universe = set(self.data)
        self._null = {}

    def observed(self, gene):
        return self.data[gene][0] if gene in self.data else 0

    def _gene_null(self, gene):
        """Null distribution (pmf) and mean of the non-synonymous mutations of a gene"""
        null = self._null.get(gene)
        if null is None:
            x, s, e, es, mu = self.data[gene]
            if self.method == 'nb':
                # posterior local rate: Gamma(theta + s, rate = theta / mu_rate + es), mu_rate = mu / es
                if np.isinf(self.theta):
                    mean = e * mu / es
                    dist = sps.poisson(mean)
                else:
                    shape = self.theta + s
                    rate = self.theta * es / mu + es
                    dist = sps.nbinom(shape, rate / (rate + e))
                    mean = e * shape / rate
                top = int(max(x, dist.ppf(1 - 1e-12))) + 1
                pmf = dist.pmf(np.arange(top + 1))
            else:
                m = x + s
                pmf = binomial_pmf(m, e / (e + es))
                mean = m * e / (e + es)
            null = (pmf, float(mean))
            self._null[gene] = null
        return null

    def test(self, genes):
        genes = [g for g in genes if g in self.data]
        if self.method == 'conditional':
            genes_tested = [g for g in genes if self.data[g][0] + self.data[g][1] > 0]
        else:
            genes_tested = genes
        observed = int(sum(self.data[g][0] for g in genes_tested))
        if len(genes_tested) == 0:
            return dict(N_GENES=len(genes), N_GENES_ALTERED=0, OBSERVED=0, EXPECTED=0.0, RATIO=np.nan, P_VALUE=np.nan)
        nulls = [self._gene_null(g) for g in genes_tested]
        pmf = sum_pmfs([n[0] for n in nulls])
        expected = sum(n[1] for n in nulls)
        return dict(
            N_GENES=len(genes),
            N_GENES_ALTERED=sum(1 for g in genes_tested if self.data[g][0] > 0),
            OBSERVED=observed,
            EXPECTED=float(expected),
            RATIO=observed / expected if expected > 0 else np.nan,
            P_VALUE=upper_tail(pmf, observed),
        )


class EventLayer:
    """Recurrence of omics events of the genes of a set"""

    def __init__(self, layer, observed, expected):
        self.layer = layer
        self.data = {g: (int(o), float(e)) for g, o, e in zip(observed.index, observed.values, expected.values)
                     if np.isfinite(e)}
        self.universe = set(self.data)

    def observed(self, gene):
        return self.data[gene][0] if gene in self.data else 0

    def test(self, genes):
        genes = [g for g in genes if g in self.data]
        observed = sum(self.data[g][0] for g in genes)
        expected = sum(self.data[g][1] for g in genes)
        if len(genes) == 0:
            return dict(N_GENES=0, N_GENES_ALTERED=0, OBSERVED=0, EXPECTED=0.0, RATIO=np.nan, P_VALUE=np.nan)
        if expected > 0:
            p = float(sps.poisson.sf(observed - 1, expected))
        else:
            p = 1.0 if observed == 0 else 0.0
        return dict(
            N_GENES=len(genes),
            N_GENES_ALTERED=sum(1 for g in genes if self.data[g][0] > 0),
            OBSERVED=int(observed),
            EXPECTED=float(expected),
            RATIO=observed / expected if expected > 0 else np.nan,
            P_VALUE=min(1.0, p),
        )


def silencing_layer(results, threshold=0.1):
    """Epigenetic silencing layer from the methylation-analysis results"""
    df = results[results['STATUS'] == 'tested']
    functional = df['FUNCTIONAL'].astype(str).str.lower()
    # hypermethylation without loss of expression is not silencing
    df = df[functional != 'false'].set_index('SYMBOL')
    significant = set(df.index[df['Q_VALUE'] < threshold])
    return EventLayer('silencing', df['HYPER_SAMPLES'], df['EXPECTED_HYPER']), significant


def expression_layers(results, threshold=0.1):
    """Over- and under-expression layers from the expression-analysis results"""
    df = results[results['STATUS'] == 'tested'].set_index('SYMBOL')
    layers = {}
    for direction in ['over', 'under']:
        significant = set(df.index[(df['Q_VALUE'] < threshold) & (df['DIRECTION'] == direction)])
        layer = EventLayer(f'expression_{direction}', df[f'{direction.upper()}_SAMPLES'],
                           df[f'EXPECTED_{direction.upper()}'])
        layers[f'expression_{direction}'] = (layer, significant)
    return layers
