"""
Synthetic outputs of the driver identification methods (and of the omics
analyses) used to test the combination step.
"""

import os

import numpy as np
import pandas as pd


def bh(p):
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    mask = np.isfinite(p)
    pv = p[mask]
    n = len(pv)
    order = np.argsort(pv)
    ranked = np.minimum.accumulate((pv[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1)
    q[mask] = out
    return q


def make_datasets(folder, n_genes=300, n_cgc=60):
    """Minimal IntOGen datasets needed by the combination"""
    genes = [f'GENE{i}' for i in range(n_genes)]
    os.makedirs(os.path.join(folder, 'cgc'), exist_ok=True)
    os.makedirs(os.path.join(folder, 'others'), exist_ok=True)
    os.makedirs(os.path.join(folder, 'regions'), exist_ok=True)
    pd.DataFrame({
        'Gene Symbol': genes[:n_cgc],
        'cancer_type': 'PRAD',
        'Tier': 1,
        'Role in Cancer': 'TSG',
    }).to_csv(os.path.join(folder, 'cgc', 'cancer_gene_census_parsed.tsv'), sep='\t', index=False)
    with open(os.path.join(folder, 'others', 'negative_gene_set.tsv'), 'w') as fd:
        fd.write('PANCANCER\t{}\n'.format(','.join(genes[-10:])))
    return genes


def _pvalues(rng, n, signal, fraction=0.8):
    p = rng.uniform(size=n)
    hits = signal & (rng.uniform(size=n) < fraction)
    p[hits] = 10 ** -rng.uniform(3, 10, hits.sum())
    return p


def make_method_outputs(folder, seed=0, n_genes=300, omics=True):
    """
    Write the results of all methods for a synthetic cohort.

    Returns:
        dict method -> file path
    """
    rng = np.random.default_rng(seed)
    genes = np.array([f'GENE{i}' for i in range(n_genes)])
    drivers = np.zeros(n_genes, dtype=bool)
    drivers[:20] = True       # CGC drivers
    drivers[100:106] = True   # non CGC drivers
    samples = rng.integers(0, 12, n_genes)
    samples[drivers] = rng.integers(5, 40, drivers.sum())
    candidates = samples >= 2
    os.makedirs(folder, exist_ok=True)
    files = {}

    def write(method, df, name, compression=None):
        path = os.path.join(folder, name)
        df.to_csv(path, sep='\t', index=False, compression=compression)
        files[method] = path

    p = _pvalues(rng, n_genes, drivers)
    q = np.where(candidates, bh(np.where(candidates, p, np.nan)), np.nan)
    write('oncodrivefml', pd.DataFrame({
        'GENE_ID': [f'ENSG{i:011d}' for i in range(n_genes)], 'SYMBOL': genes,
        'MUTS': samples + 1, 'MUTS_RECURRENCE': samples, 'SAMPLES': samples,
        'P_VALUE': p, 'Q_VALUE': q}), 'fml.tsv.gz', compression='gzip')

    p = _pvalues(rng, n_genes, drivers)
    write('dndscv', pd.DataFrame({
        'gene_name': genes, 'n_syn': 1, 'n_mis': samples, 'n_non': samples // 3, 'n_spl': 0, 'n_ind': 1,
        'wmis_cv': rng.uniform(0.5, 5, n_genes), 'wnon_cv': rng.uniform(0.5, 5, n_genes),
        'wspl_cv': 1.0, 'wind_cv': 1.0, 'pallsubs_cv': p, 'qallsubs_cv': bh(p)}),
        'dndscv.tsv.gz', compression='gzip')

    p = _pvalues(rng, n_genes, drivers)
    write('oncodriveclustl', pd.DataFrame({'ENSID': genes, 'SYMBOL': genes, 'P_ANALYTICAL': p,
                                           'Q_ANALYTICAL': bh(p)}), 'clustl.txt')
    p = _pvalues(rng, n_genes, drivers, 0.4)
    write('hotmaps', pd.DataFrame({'GENE': genes, 'Min p-value': p, 'q-value': bh(p)}), 'hotmaps.tsv')
    p = _pvalues(rng, n_genes, drivers, 0.5)
    write('smregions', pd.DataFrame({'HUGO_SYMBOL': genes, 'P_VALUE': p, 'Q_VALUE': bh(p),
                                     'REGION': [f'ENST{i}:PF1:1:10' for i in range(n_genes)]}), 'smregions.tsv')
    p = _pvalues(rng, n_genes, drivers)
    write('cbase', pd.DataFrame({'gene': genes, 'p_pos': p, 'q_pos': bh(p)}), 'cbase.tsv')
    p = _pvalues(rng, n_genes, drivers)
    write('mutpanning', pd.DataFrame({'Name': genes, 'Significance': p, 'FDR': bh(p)}), 'mutpanning.txt')

    if omics:
        # omics results include genes without mutations (not candidates) and untested genes
        extra = np.array([f'NOMUT{i}' for i in range(50)])
        all_genes = np.concatenate([genes, extra])
        signal = np.concatenate([drivers, np.ones(50, dtype=bool)])
        for method in ['methylation', 'expression']:
            p = _pvalues(rng, len(all_genes), signal, 0.5)
            p[rng.uniform(size=len(p)) < 0.1] = np.nan
            write(method, pd.DataFrame({'SYMBOL': all_genes, 'STATUS': 'tested', 'P_VALUE': p, 'Q_VALUE': bh(p)}),
                  f'{method}.tsv.gz', compression='gzip')
    return files
