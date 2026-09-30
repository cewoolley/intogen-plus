"""
Co-occurrence: common practice (pairwise Fisher's exact test on the most mutated genes, as in
maftools somaticInteractions / cBioPortal) vs the fork (burden-aware DISCOVER-like null with
Tarone-BH), on real TCGA cohorts and with co-occurrence planted in real tumours.
"""
import itertools
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats as sps
from scipy.optimize import brentq
from scipy.special import expit, logit

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'core'))
from intogen_core.omics.stats import fdr_bh, poisson_binomial_tail   # noqa: E402
from intogen_core.pathways import analysis, cooccurrence             # noqa: E402

COHORTS = ['KIRC', 'LUAD', 'BRCA', 'COAD']
CONTROLS = {
    'LUAD': [('STK11', 'KEAP1', 'co-occurrence (KL subtype)'), ('KRAS', 'STK11', 'co-occurrence'),
             ('KRAS', 'EGFR', 'mutual exclusivity')],
    'BRCA': [('CDH1', 'PIK3CA', 'co-occurrence (lobular)'), ('TP53', 'PIK3CA', 'mutual exclusivity'),
             ('TP53', 'CDH1', 'mutual exclusivity')],
    'COAD': [('APC', 'KRAS', 'co-occurrence'), ('BRAF', 'KRAS', 'mutual exclusivity'),
             ('BRAF', 'APC', 'mutual exclusivity')],
    'KIRC': [('VHL', 'PBRM1', 'co-occurrence (3p loss)'), ('PBRM1', 'BAP1', 'mutual exclusivity')],
}


def load(cohort):
    tumours, hits = analysis.mutation_hits(f'results/{cohort}/mutations.tsv.gz')
    drivers = set(pd.read_csv(f'results/{cohort}/vet.tsv', sep='\t')['SYMBOL'])
    return tumours, hits['mutation'][['SYMBOL', 'SAMPLE']], drivers


def fisher_pairs(layer, genes):
    rows = []
    n = len(layer.tumours)
    for a, b in itertools.combinations(genes, 2):
        h1 = layer.hits[layer.gene_index[a]]
        h2 = layer.hits[layer.gene_index[b]]
        both = int((h1 & h2).sum())
        table = [[both, int(h1.sum()) - both], [int(h2.sum()) - both, n - int((h1 | h2).sum())]]
        odds, p = sps.fisher_exact(table)                 # two-sided, as maftools
        rows.append((a, b, both, h1.sum() * h2.sum() / n, odds, p))
    df = pd.DataFrame(rows, columns=['A', 'B', 'BOTH', 'EXPECTED_INDEPENDENT', 'ODDS', 'P'])
    df['Q'] = fdr_bh(df['P'].values)
    return df


def discover_pairs(layer, genes, iters=300):
    """DISCOVER (Canisius et al. 2016): logit p_gs = a_g + b_s, keeping the margins of genes and tumours"""
    hits = layer.hits.astype(float)
    r, c = hits.sum(axis=1), hits.sum(axis=0)
    ok = (r > 0) & (r < hits.shape[1])
    a = np.full(len(r), -np.inf)
    a[ok] = logit(r[ok] / hits.shape[1])
    b = np.zeros(hits.shape[1])
    for _ in range(iters):
        q = expit(a[ok, None] + b[None, :])
        a[ok] += np.clip((r[ok] - q.sum(axis=1)) / np.maximum((q * (1 - q)).sum(axis=1), 1e-12), -5, 5)
        q = expit(a[ok, None] + b[None, :])
        b += np.clip((c - q.sum(axis=0)) / np.maximum((q * (1 - q)).sum(axis=0), 1e-12), -5, 5)
    rows = []
    for g1, g2 in itertools.combinations(genes, 2):
        i, j = layer.gene_index[g1], layer.gene_index[g2]
        probs = expit(a[i] + b) * expit(a[j] + b)
        tail = np.append(poisson_binomial_tail(probs), 0.0)
        both = int((layer.hits[i] & layer.hits[j]).sum())
        rows.append((g1, g2, float(probs.sum()), float(tail[both]),
                     float(tail[min(int(layer.hits[i].sum()), int(layer.hits[j].sum()))])))
    df = pd.DataFrame(rows, columns=['A', 'B', 'EXPECTED_BOTH', 'P_VALUE', 'P_MIN'])
    df['Q_VALUE'], df['TESTABLE'] = cooccurrence.tarone_bh(df['P_VALUE'].values, df['P_MIN'].values, 0.1)
    return df


def real_data():
    summary = {}
    for cohort in COHORTS:
        tumours, pairs, drivers = load(cohort)
        layer = cooccurrence.Layer('mutation', pairs, tumours)
        counts = pairs.drop_duplicates().groupby('SYMBOL').size().sort_values(ascending=False)
        top = list(counts.index[:25])
        genes = sorted(set(top) | (drivers & set(counts.index[counts >= 3])))
        fisher = fisher_pairs(layer, genes)
        disc = discover_pairs(layer, genes)
        elas = elastic_pairs(layer, genes)
        def kind(a, b):
            return {0: 'passenger-passenger', 1: 'driver-passenger', 2: 'driver-driver'}[(a in drivers) + (b in drivers)]
        fisher['KIND'] = [kind(a, b) for a, b in zip(fisher['A'], fisher['B'])]
        disc['KIND'] = [kind(a, b) for a, b in zip(disc['A'], disc['B'])]
        elas['KIND'] = [kind(a, b) for a, b in zip(elas['A'], elas['B'])]
        e_sig = elas[elas['Q_VALUE'] < 0.1]
        f_sig = fisher[(fisher['P'] < 0.05) & (fisher['ODDS'] > 1)]
        f_sigq = fisher[(fisher['Q'] < 0.1) & (fisher['ODDS'] > 1)]
        d_sig = disc[disc['Q_VALUE'] < 0.1]
        res = {'genes': len(genes), 'drivers_tested': len(set(genes) & drivers), 'pairs': len(fisher),
               'fisher_p05_cooccurring': len(f_sig), 'fisher_q10_cooccurring': len(f_sigq),
               'fisher_q10_by_kind': f_sigq['KIND'].value_counts().to_dict(),
               'discover_q10': len(d_sig), 'discover_q10_by_kind': d_sig['KIND'].value_counts().to_dict(),
               'discover_significant_pairs': [f'{a}-{b}' for a, b in zip(d_sig['A'], d_sig['B'])],
               'elastic_q10': len(e_sig), 'elastic_q10_by_kind': e_sig['KIND'].value_counts().to_dict(),
               'elastic_significant_pairs': [f'{a}-{b}' for a, b in zip(e_sig['A'], e_sig['B'])],
               'controls': []}
        for a, b, label in CONTROLS[cohort]:
            f = fisher[((fisher.A == a) & (fisher.B == b)) | ((fisher.A == b) & (fisher.B == a))]
            d = disc[((disc.A == a) & (disc.B == b)) | ((disc.A == b) & (disc.B == a))]
            e = elas[((elas.A == a) & (elas.B == b)) | ((elas.A == b) & (elas.B == a))]
            if len(f) and len(d):
                f, d, e = f.iloc[0], d.iloc[0], e.iloc[0]
                res['controls'].append({'pair': f'{a}-{b}', 'expected': label, 'both': int(f['BOTH']),
                                        'expected_independent': round(float(f['EXPECTED_INDEPENDENT']), 1),
                                        'expected_burden_aware': round(float(d['EXPECTED_BOTH']), 1),
                                        'fisher_p': float(f['P']), 'fisher_odds': float(f['ODDS']),
                                        'discover_p': float(d['P_VALUE']), 'discover_q': float(d['Q_VALUE']),
                                        'expected_elastic': round(float(e['EXPECTED_BOTH']), 1),
                                        'elastic_p': float(e['P_VALUE']), 'elastic_q': float(e['Q_VALUE'])})
        summary[cohort] = res
        print(cohort, json.dumps({k: v for k, v in res.items() if k != 'controls'}), flush=True)
        for c in res['controls']:
            print('   ', c, flush=True)
    return summary


def tumour_offsets(layer, iters=300):
    """b_s of the additive background model (logit p_gs = a_g + b_s) fitted on the genes altered in 1-99% tumours"""
    hits = layer.hits.astype(float)
    r = hits.sum(axis=1)
    hits = hits[(r > 0) & (r < hits.shape[1])]
    r, c = hits.sum(axis=1), hits.sum(axis=0)
    a = logit(r / hits.shape[1])
    b = np.zeros(hits.shape[1])
    for _ in range(iters):
        q = expit(a[:, None] + b[None, :])
        a += np.clip((r - q.sum(axis=1)) / np.maximum((q * (1 - q)).sum(axis=1), 1e-12), -5, 5)
        q = expit(a[:, None] + b[None, :])
        b += np.clip((c - q.sum(axis=0)) / np.maximum((q * (1 - q)).sum(axis=0), 1e-12), -5, 5)
    return b


def elastic_probs(h, log_burden):
    """P(event in each tumour) from a logistic regression on the tumour burden (own elasticity)"""
    X = np.column_stack([np.ones(len(h)), log_burden])
    beta = np.array([logit(np.clip(h.mean(), 1e-3, 1 - 1e-3)), 0.0])
    for _ in range(100):
        p = expit(X @ beta)
        step = np.linalg.solve((X * (p * (1 - p))[:, None]).T @ X + 1e-6 * np.eye(2), X.T @ (h - p))
        beta += step
        if np.abs(step).max() < 1e-10:
            break
    return expit(X @ beta)


def elastic_test(h1, h2, log_burden):
    q = elastic_probs(h1.astype(float), log_burden) * elastic_probs(h2.astype(float), log_burden)
    tail = poisson_binomial_tail(q)
    both = int((h1 & h2).sum())
    return float(tail[both]), float(tail[min(int(h1.sum()), int(h2.sum()))]), float(q.sum())


def elastic_pairs(layer, genes):
    log_burden = np.log1p(layer.hits.sum(axis=0))
    rows = []
    for a, b in itertools.combinations(genes, 2):
        h1, h2 = layer.hits[layer.gene_index[a]], layer.hits[layer.gene_index[b]]
        # burden without the two tested genes
        lb = np.log1p(layer.hits.sum(axis=0) - h1 - h2)
        p, pmin, exp = elastic_test(h1, h2, lb)
        rows.append((a, b, p, pmin, exp))
    df = pd.DataFrame(rows, columns=['A', 'B', 'P_VALUE', 'P_MIN', 'EXPECTED_BOTH'])
    df['Q_VALUE'], df['TESTABLE'] = cooccurrence.tarone_bh(df['P_VALUE'].values, df['P_MIN'].values, 0.1)
    return df


def synthetic_gene(b, n):
    a = brentq(lambda x: expit(x + b).sum() - n, -50, 50)
    return expit(a + b)


def planted(cohort, alphas=(0.0, 0.5, 1.0), psis=(1.0, 2.0, 3.0), freq=0.15, reps=300):
    """Two synthetic genes mutated in ~15% of the real tumours, with a tumour preference burden^alpha"""
    rng = np.random.default_rng(11)
    tumours, pairs, drivers = load(cohort)
    layer = cooccurrence.Layer('mutation', pairs, tumours)
    b = tumour_offsets(layer)
    burden = layer.hits.sum(axis=0).astype(float)
    log_burden = np.log1p(burden)
    n = int(freq * len(tumours))
    rows = []
    for alpha in alphas:
        w = burden ** alpha
        for psi in psis:
            fis, dis, ela = [], [], []
            for _ in range(reps):
                x = np.zeros(len(tumours), bool)
                x[rng.choice(len(tumours), n, replace=False, p=w / w.sum())] = True
                wy = w * np.where(x, psi, 1.0)
                y = np.zeros(len(tumours), bool)
                y[rng.choice(len(tumours), n, replace=False, p=wy / wy.sum())] = True
                both = int((x & y).sum())
                fis.append(sps.fisher_exact([[both, n - both], [n - both, len(tumours) - 2 * n + both]],
                                            alternative='greater')[1])
                q = synthetic_gene(b, n) ** 2              # both synthetic genes have the same margins
                dis.append(poisson_binomial_tail(q)[both])
                ela.append(elastic_test(x, y, log_burden)[0])
            fis, dis, ela = np.array(fis), np.array(dis), np.array(ela)
            rows.append({'COHORT': cohort, 'ALPHA': alpha, 'PSI': psi, 'FISHER_P05': float((fis < 0.05).mean()),
                         'DISCOVER_P05': float((dis < 0.05).mean()), 'ELASTIC_P05': float((ela < 0.05).mean())})
            print('   ', rows[-1], flush=True)
    return rows


if __name__ == '__main__':
    out = {'real': real_data(), 'planted': []}
    for cohort in ['LUAD', 'COAD']:
        out['planted'] += planted(cohort)
    json.dump(out, open('results/cooccurrence.json', 'w'), indent=2, default=float)
