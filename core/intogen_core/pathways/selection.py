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

    Background dN/dS (competitive null, ``nb``): in real exomes, the dN/dS of
    genes that are not under positive selection departs from 1 in a way that
    depends on their genomic context (e.g. genes not expressed in the tissue,
    with no purifying selection nor transcription-coupled repair, have dN/dS
    of 1.1-1.2 with the global mutational model). Summed over hundreds of genes,
    these small deviations make large gene sets significant (neuronal, muscle or
    extracellular matrix genes). The expected non-synonymous mutations of each
    gene are therefore multiplied by the background dN/dS of genes with a
    similar context, predicted with a Poisson regression on genomic covariates
    (by default, the epigenomic covariates used by dNdScv) fitted on the genes
    that are not individually significant. Gene sets are thus compared with
    comparable genes rather than with the neutral model.

Omics layers (``silencing``, ``expression_over``, ``expression_under``)
    The number of events (e.g. hypermethylated tumours) of the genes of a set
    is compared with the number expected given the event rate of each tumour
    (the sum of the per-gene expectations of the omics analyses). Since the
    events are many independent rare events, their sum is approximated by a
    Poisson distribution, which is conservative.
"""

import numpy as np
import pandas as pd
from scipy import optimize
from scipy import stats as sps

from intogen_core.omics.stats import binomial_pmf, sum_pmfs, upper_tail


MUTATION_LAYERS = {
    'mutation': (['n_mis', 'n_non', 'n_spl'], ['exp_mis', 'exp_non', 'exp_spl']),
    'mutation_missense': (['n_mis'], ['exp_mis']),
    'mutation_truncating': (['n_non', 'n_spl'], ['exp_non', 'exp_spl']),
}

# classes of mutations with their own background dN/dS
BACKGROUND_CLASSES = {'mutation_missense': ['exp_mis'], 'mutation_truncating': ['exp_non', 'exp_spl']}

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

    def neutral_mean(self, gene):
        """Expected non-synonymous mutations of a gene under neutrality given its local rate (nb)"""
        x, s, e, es, mu = self.data[gene]
        if np.isinf(self.theta):
            return e * mu / es
        return e * (self.theta + s) / (self.theta * es / mu + es)

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
        genes = sorted(g for g in genes if g in self.data)
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
        genes = sorted(g for g in genes if g in self.data)
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


def poisson_regression(y, offset, X, ridge=1.0, max_iter=100):
    """Poisson regression with an offset, an intercept and ridge-penalised coefficients (Newton-Raphson)"""
    X = np.column_stack([np.ones(len(y)), X])
    penalty = ridge * np.diag(np.r_[0.0, np.ones(X.shape[1] - 1)])
    beta = np.zeros(X.shape[1])
    beta[0] = np.log(max(y.sum(), 0.5) / np.exp(offset).sum())
    for _ in range(max_iter):
        mu = np.exp(np.clip(offset + X @ beta, -50, 50))
        step = np.linalg.solve((X * mu[:, None]).T @ X + penalty, X.T @ (y - mu) - penalty @ beta)
        beta += step
        if np.abs(step).max() < 1e-8:
            break
    return beta


def genemuts_features(genemuts):
    """Covariates available in genemuts: covariate-predicted relative mutation rate and gene size"""
    rel = np.log(genemuts['exp_syn_cv'] / genemuts['exp_syn'])
    size = np.log(genemuts['exp_syn'])
    return pd.DataFrame({'rate': rel, 'size': size, 'rate2': rel ** 2, 'size2': size ** 2, 'rate_size': rel * size})


def background_omega(genemuts, exclude=(), covariates=None, ridge=1.0, bounds=(0.5, 2.0)):
    """
    Background dN/dS of each gene given its genomic context.

    Args:
        genemuts: dNdScv genemuts (with exp_syn_cv)
        exclude: genes left out of the fit (individually significant)
        covariates: DataFrame of numeric covariates indexed by gene (default: genemuts_features)

    Returns:
        genemuts with the expected non-synonymous mutations multiplied by the background dN/dS,
        and a dict describing the fit
    """
    genes = genemuts['gene_name'].values
    if covariates is None:
        features, source = genemuts_features(genemuts), 'genemuts'
        features.index = genes
    else:
        features, source = covariates.reindex(genes), 'covariates'
    features = features.replace([np.inf, -np.inf], np.nan)
    missing = features.isna().any(axis=1).values
    features = ((features - features.mean()) / features.std().replace(0, 1)).fillna(0.0)
    info = {'source': source, 'features': int(features.shape[1]), 'genes_without_covariates': int(missing.sum())}

    adjusted = genemuts.copy()
    excluded = set(exclude)
    for name, columns in BACKGROUND_CLASSES.items():
        layer = MutationLayer(genemuts, name, method='nb')
        offset = np.full(len(genes), np.nan)
        observed = np.zeros(len(genes))
        for i, g in enumerate(genes):
            if g in layer.data:
                offset[i] = np.log(layer.neutral_mean(g)) if layer.neutral_mean(g) > 0 else np.nan
                observed[i] = layer.data[g][0]
        fit = np.isfinite(offset) & ~missing & ~np.isin(genes, list(excluded))
        beta = poisson_regression(observed[fit], offset[fit], features.values[fit], ridge=ridge)
        omega = np.clip(np.exp(np.column_stack([np.ones(len(genes)), features.values]) @ beta), *bounds)
        for c in columns:
            adjusted[c] = adjusted[c] * omega
        info[name] = {'genes_fitted': int(fit.sum()),
                      'omega_quantiles': [float(q) for q in np.quantile(omega[fit], [0.05, 0.5, 0.95])]}
    return adjusted, info
