"""
Detectable effect sizes for co-occurrence / exclusivity between two driver-like events in MSS CRC, with the
fork's burden-elastic test, at the TCGA size (469 MSS tumours) and at the size of Cornish et al. (2,023 tumours,
resampling TCGA MSS tumours). Threshold p < 5e-5 ~ the first discovery of a BH correction over ~2,000 pairs.
"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'core'))
from intogen_core.omics.stats import poisson_binomial_tail     # noqa: E402
from intogen_core.pathways import analysis, cooccurrence        # noqa: E402

rng = np.random.default_rng(3)
tumours, hits = analysis.mutation_hits('results/MSS/mutations.tsv.gz')
layer = cooccurrence.Layer('mutation', hits['mutation'][['SYMBOL', 'SAMPLE']], tumours)
burden_all = layer.burden.astype(float)
THRESHOLD = 5e-5
REPS = 100


def detect(n, f1, f2, psi):
    burden = rng.choice(burden_all, n, replace=True) if n != len(burden_all) else burden_all
    lb = np.log1p(burden)
    co, ex = 0, 0
    for _ in range(REPS):
        x = rng.uniform(size=n) < f1
        base = f2 / (1 + f1 * (psi - 1))           # keeps the marginal frequency of y close to f2
        py = np.clip(np.where(x, base * psi, base), 0, 1)
        y = rng.uniform(size=n) < py
        if x.sum() < 3 or y.sum() < 3:
            continue
        q = cooccurrence.event_probabilities(x, lb) * cooccurrence.event_probabilities(y, lb)
        tail = np.append(poisson_binomial_tail(q), 0.0)
        both = int((x & y).sum())
        co += tail[both] < THRESHOLD
        ex += (1 - tail[both + 1]) < THRESHOLD
    return co / REPS, ex / REPS


rows = []
for n in [469, 2023]:
    for f1, f2 in [(0.05, 0.05), (0.1, 0.1), (0.2, 0.2), (0.1, 0.4)]:
        for psi in [0.2, 0.4, 0.6, 1.5, 2, 3, 4, 6]:
            co, ex = detect(n, f1, f2, psi)
            rows.append({'n': n, 'f1': f1, 'f2': f2, 'odds': psi, 'power': ex if psi < 1 else co})
        print(n, f1, f2, [(r['odds'], r['power']) for r in rows if r['n'] == n and r['f1'] == f1 and r['f2'] == f2], flush=True)
df = pd.DataFrame(rows)
df.to_csv('results/power.tsv', sep='\t', index=False)
