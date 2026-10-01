"""
Power of the two primary hypotheses in a replication cohort (PREREGISTRATION.md, section Power).

H1: truncating substitutions in the nine co-regulators against the neutral expectation, which grows with the number
    of MSS tumours (TCGA: 5.1 expected in 469 exomes). Poisson approximation of the fork's null.
H2: stage IV in carriers versus other MSS tumours, normal approximation of the log odds ratio from the expected
    2 x 2 table, checked by simulation for the planning scenario.
Both one-sided at spec.ALPHA (fixed sequence: H1, then H2 if H1 is supported).
"""
import os
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spec  # noqa: E402

ALPHA = spec.ALPHA
EXPECTED_PER_TUMOUR = 5.1 / 469          # TCGA MSS, neutral truncating substitutions in the nine genes


def power_h1(n, ratio):
    e0 = EXPECTED_PER_TUMOUR * n
    crit = int(stats.poisson.ppf(1 - ALPHA, e0)) + 1          # smallest count with P(X >= crit) <= alpha
    while stats.poisson.sf(crit - 1, e0) > ALPHA:
        crit += 1
    return float(stats.poisson.sf(crit - 1, ratio * e0))


def power_h2(n, carrier_rate, p0, odds):
    p1 = odds * p0 / (1 - p0 + odds * p0)
    n1, n0 = n * carrier_rate, n * (1 - carrier_rate)
    se = np.sqrt(1 / (n1 * p1) + 1 / (n1 * (1 - p1)) + 1 / (n0 * p0) + 1 / (n0 * (1 - p0)))
    return float(stats.norm.cdf(np.log(odds) / se - stats.norm.ppf(1 - ALPHA)))


def simulate_h2(n, carrier_rate, p0, odds, reps=2000, seed=spec.SEED):
    rng = np.random.default_rng(seed)
    p1 = odds * p0 / (1 - p0 + odds * p0)
    hits = 0
    for _ in range(reps):
        x = (rng.uniform(size=n) < carrier_rate).astype(float)
        y = (rng.uniform(size=n) < np.where(x == 1, p1, p0)).astype(float)
        fit = sm.Logit(y, sm.add_constant(x)).fit(disp=0)
        hits += (fit.params[1] > 0) and (fit.pvalues[1] / 2 < ALPHA)
    return hits / reps


def main():
    rows = []
    for n in [1000, 1400, 1700]:
        for ratio in [1.5, 2, 3]:
            rows.append({'hypothesis': 'H1', 'mss_tumours': n, 'ratio': ratio, 'power': power_h1(n, ratio)})
        for c in [0.06, 0.09]:
            for p0 in [0.05, 0.08, 0.12]:
                for odds in [1.5, 2.0, 2.8]:
                    rows.append({'hypothesis': 'H2', 'mss_tumours': n, 'carrier_rate': c, 'stage_iv_others': p0,
                                 'odds_ratio': odds, 'power': power_h2(n, c, p0, odds)})
    df = pd.DataFrame(rows)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
    os.makedirs(out, exist_ok=True)
    df.to_csv(os.path.join(out, 'power.tsv'), sep='\t', index=False, float_format='%.3f')
    print(df[df.hypothesis == 'H1'].pivot(index='mss_tumours', columns='ratio', values='power').round(3))
    h2 = df[df.hypothesis == 'H2']
    print(h2.pivot_table(index=['mss_tumours', 'carrier_rate', 'stage_iv_others'], columns='odds_ratio', values='power').round(2))
    check = simulate_h2(1400, 0.09, 0.08, 2.0)
    print(f'simulation check, 1,400 MSS tumours, 9% carriers, 8% stage IV, OR 2.0: {check:.2f} '
          f'(approximation {power_h2(1400, 0.09, 0.08, 2.0):.2f})')


if __name__ == '__main__':
    main()
