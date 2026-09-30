"""Calibration of gene-set selection tests on real cohorts, with random gene sets of non-driver genes"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats as sps

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'core'))
from intogen_core.omics.stats import simes_combine            # noqa: E402
from intogen_core.pathways.selection import MUTATION_LAYERS, MutationLayer   # noqa: E402

COHORTS = ['KIRC', 'LUAD', 'BRCA', 'COAD']


def naive_poisson(gm, genes):
    """Covariate-predicted neutral expectation without overdispersion nor synonymous update"""
    sub = gm.loc[[g for g in genes if g in gm.index]]
    exp = ((sub['exp_mis'] + sub['exp_non'] + sub['exp_spl']) * sub['exp_syn_cv'] / sub['exp_syn']).sum()
    obs = (sub['n_mis'] + sub['n_non'] + sub['n_spl']).sum()
    return float(sps.poisson.sf(obs - 1, exp))


def tests(cohort, random_sets):
    gm = pd.read_csv(f'data/{cohort}/genemuts.tsv', sep='\t')
    adj = pd.read_csv(f'results/{cohort}/genemuts_final.tsv', sep='\t')
    layers = {
        'shipped': {name: MutationLayer(gm, name) for name in MUTATION_LAYERS},
        'final': {name: MutationLayer(adj, name) for name in MUTATION_LAYERS},
        'conditional': {'mutation': MutationLayer(gm, 'mutation', method='conditional')},
    }
    gmi = gm.set_index('gene_name')
    rows = []
    for i, genes in enumerate(random_sets):
        row = {'SET': i, 'N_GENES': len(genes)}
        for name in ['shipped', 'final']:
            ps = {layer: L.test(genes)['P_VALUE'] for layer, L in layers[name].items()}
            row[f'{name}_mutation'] = ps['mutation']
            row[f'{name}_combined'] = simes_combine(list(ps.values()))
            row[f'{name}_ratio'] = layers[name]['mutation'].test(genes)['RATIO']
        row['conditional'] = layers['conditional']['mutation'].test(genes)['P_VALUE']
        row['naive_poisson'] = naive_poisson(gmi, genes)
        rows.append(row)
    df = pd.DataFrame(rows)
    gsd_file = f'results/{cohort}/gsd.tsv'
    if os.path.exists(gsd_file):
        gsd = pd.read_csv(gsd_file, sep='\t')
        gsd = gsd[gsd['UNIT'].str.startswith('RND|')].copy()
        gsd['SET'] = gsd['UNIT'].str.split('|').str[1].astype(int)
        # one-sided (excess of non-synonymous mutations) from the two-sided LRT
        gsd['genesetdnds'] = np.where(gsd['wall'] > 1, gsd['p_wall'] / 2, 1 - gsd['p_wall'] / 2)
        df = df.merge(gsd[['SET', 'genesetdnds']], on='SET', how='left')
    return df


if __name__ == '__main__':
    summary = {}
    for cohort in COHORTS:
        random_sets = json.load(open(f'results/{cohort}/random_sets.json'))
        df = tests(cohort, random_sets)
        df.to_csv(f'results/{cohort}/calibration_final.tsv', sep='\t', index=False)
        methods = [c for c in ['naive_poisson', 'conditional', 'genesetdnds', 'shipped_mutation', 'shipped_combined',
                               'final_mutation', 'final_combined'] if c in df.columns]
        big = df['N_GENES'] >= df['N_GENES'].quantile(2 / 3)
        summary[cohort] = {m: {'p05': float((df[m] < 0.05).mean()), 'p01': float((df[m] < 0.01).mean()),
                               'p05_large_sets': float((df.loc[big, m] < 0.05).mean()),
                               'n': int(df[m].notna().sum())} for m in methods}
        print(cohort, ' '.join(f"{m}: {v['p05']:.3f}/{v['p01']:.3f} (large {v['p05_large_sets']:.2f})"
                               for m, v in summary[cohort].items()), flush=True)
        print('   median random-set ratio shipped %.3f, bgomega %.3f' % (df['shipped_ratio'].median(),
                                                                         df['final_ratio'].median()))
    json.dump(summary, open('results/calibration_final.json', 'w'), indent=2)
