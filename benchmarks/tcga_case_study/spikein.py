"""
Semi-synthetic power: selection spread over the long tail of a gene set, spiked into real cohorts.

A random set of 40 non-driver genes receives extra missense (or truncating) mutations,
omega-fold their neutral expectation. Detection (q < 0.1 among all the Reactome/hallmark sets):

- pre-fork: gene-level significant genes (q < 0.1 over all genes), then over-representation (ORA)
- fork as shipped: long-tail gene-set test (self-contained null)
- fork with covariate-matched background
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats as sps

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'core'))
from intogen_core.omics.stats import fdr_bh, simes_combine          # noqa: E402
from intogen_core.pathways.genesets import GeneSets, read_gene_sets   # noqa: E402
from intogen_core.pathways.selection import MUTATION_LAYERS, MutationLayer   # noqa: E402

COHORTS = ['KIRC', 'LUAD', 'BRCA', 'COAD']
SCENARIOS = [('missense', 1.3), ('missense', 1.6), ('missense', 2.0), ('truncating', 2.0), ('truncating', 3.0)]
REPS = 25
SIZE = 40


def gene_pvalues(layers, genes):
    return {g: simes_combine([L.test([g])['P_VALUE'] for L in layers.values()]) for g in genes}


def ora(hits, universe, sets):
    n, N = len(hits), len(universe)
    out = {}
    for s, genes in sets.items():
        K, k = len(genes), len(genes & hits)
        out[s] = float(sps.hypergeom.sf(k - 1, N, K, n)) if k > 0 else 1.0
    return out


def run(cohort):
    rng = np.random.default_rng(7)
    gm = pd.read_csv(f'data/{cohort}/genemuts.tsv', sep='\t')
    adj = pd.read_csv(f'results/{cohort}/genemuts_final.tsv', sep='\t')
    drivers = set(pd.read_csv(f'results/{cohort}/vet.tsv', sep='\t')['SYMBOL'])
    base = {name: MutationLayer(gm, name) for name in MUTATION_LAYERS}
    theta = base['mutation'].theta
    universe = sorted(base['mutation'].universe)
    sets = GeneSets(read_gene_sets('data/gene_sets.tsv.gz'), set(universe), 10, 500)
    collection = {s: sets.genes[s] for s in sets.ids}
    # p-values of the real gene sets (long tail, combined), as shipped and with the matched background
    real = {}
    for tag in ['fork', 'final']:
        df = pd.read_csv(f'results/{cohort}/{tag}.pathways.tsv.gz', sep='\t')
        df = df[(df['SCOPE'] == 'long_tail') & (df['LAYER'] == 'combined')]
        real[tag] = dict(zip(df['SET'], df['P_VALUE']))
    # gene-level p-values of all genes (pre-fork proxy), computed once
    print(cohort, 'gene-level tests...', flush=True)
    gene_p = gene_pvalues(base, universe)
    pool = np.array([g for g in universe if g not in drivers])
    gmi, adji = gm.set_index('gene_name'), adj.set_index('gene_name')
    rows = []
    for kind, omega in SCENARIOS:
        cols = ['n_mis'] if kind == 'missense' else ['n_non', 'n_spl']
        for rep in range(REPS):
            genes = list(rng.choice(pool, SIZE, replace=False))
            spiked = gmi.copy()
            for g in genes:
                for c in cols:
                    # neutral expectation of the gene given its local rate (posterior mean of the rate)
                    e = spiked.at[g, c.replace('n_', 'exp_')]
                    shape = theta + spiked.at[g, 'n_syn']
                    rate = theta * spiked.at[g, 'exp_syn'] / spiked.at[g, 'exp_syn_cv'] + spiked.at[g, 'exp_syn']
                    spiked.at[g, c] += rng.poisson(e * shape / rate * (omega - 1))
            spiked_adj = adji.copy()
            spiked_adj.loc[genes, ['n_mis', 'n_non', 'n_spl']] = spiked.loc[genes, ['n_mis', 'n_non', 'n_spl']].values
            res = {'COHORT': cohort, 'KIND': kind, 'OMEGA': omega, 'REP': rep}
            # pre-fork: gene-level significance, then ORA of the significant genes
            layers = {name: MutationLayer(spiked.reset_index(), name, theta=theta) for name in MUTATION_LAYERS}
            gp = dict(gene_p)
            gp.update(gene_pvalues(layers, genes))
            names = list(gp)
            q = dict(zip(names, fdr_bh(np.array([gp[g] for g in names]))))
            hits = {g for g in names if q[g] < 0.1} | drivers
            res['GENES_SIGNIFICANT'] = sum(q[g] < 0.1 for g in genes)
            coll = dict(collection)
            coll['SPIKED'] = set(genes)
            ora_p = ora(hits, set(universe), coll)
            ora_q = dict(zip(ora_p, fdr_bh(np.array(list(ora_p.values())))))
            res['PREFORK_ORA'] = ora_q['SPIKED'] < 0.1
            # fork: long tail test of the spiked set among the real sets
            for tag, table in [('fork', spiked), ('final', spiked_adj)]:
                lay = {name: MutationLayer(table.reset_index(), name, theta=theta) for name in MUTATION_LAYERS}
                p = simes_combine([L.test(genes)['P_VALUE'] for L in lay.values()])
                allp = dict(real[tag])
                allp['SPIKED'] = p
                qq = dict(zip(allp, fdr_bh(np.array(list(allp.values())))))
                res[f'{tag.upper()}_P'] = p
                res[f'{tag.upper()}'] = qq['SPIKED'] < 0.1
            rows.append(res)
        sub = pd.DataFrame([r for r in rows if r['KIND'] == kind and r['OMEGA'] == omega])
        print(f'  {kind} omega={omega}: genes significant (mean) {sub.GENES_SIGNIFICANT.mean():.2f} | '
              f'pre-fork ORA {sub.PREFORK_ORA.mean():.2f} | fork shipped {sub.FORK.mean():.2f} | '
              f'fork matched {sub.FINAL.mean():.2f}', flush=True)
    return rows


if __name__ == '__main__':
    cohorts = sys.argv[1:] or COHORTS
    for cohort in cohorts:
        pd.DataFrame(run(cohort)).to_csv(f'results/{cohort}/spikein_final.tsv', sep='\t', index=False)
