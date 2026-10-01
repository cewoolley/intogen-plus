"""
Pre-registered replication of the co-regulator findings of TCGA MSS colorectal cancer (PREREGISTRATION.md).

    python replicate.py --work <work dir of run.sh> --out <results dir> [--cornish <driver table>]

Writes replication.json and summary.md with aggregate results only: no sample identifiers, and counts of tumours
below spec.MIN_CELL are suppressed.
"""
import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import scipy
import statsmodels
import statsmodels.api as sm
from scipy import stats
from statsmodels.duration.hazard_regression import PHReg

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'core'))
sys.path.insert(0, HERE)
import spec                                                                   # noqa: E402
from intogen_core.omics.stats import fdr_bh, sum_pmfs                         # noqa: E402
from intogen_core.pathways import analysis                                    # noqa: E402
from intogen_core.pathways.selection import MutationLayer, background_omega   # noqa: E402

TRUNC_LAYER = 'mutation_truncating'


# ------------------------------------------------------------------------------------------------ helpers
def cell(n):
    """Count of tumours for an aggregate output: suppressed below spec.MIN_CELL (but zero is reported)"""
    n = int(n)
    return n if n == 0 or n >= spec.MIN_CELL else f'<{spec.MIN_CELL}'


def one_sided(coef, p_two):
    return p_two / 2 if coef > 0 else 1 - p_two / 2


def sha256(path):
    with open(path, 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()


def clean(x):
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, (np.integer, np.bool_)):
        return x.item()
    if isinstance(x, (float, np.floating)):
        return None if not np.isfinite(x) else float(x)
    return x


def stage_group(stage):
    """Stage group I-IV from registry or AJCC strings (Stage IIIB, IVA, III...); anything else is missing"""
    s = stage.astype(str).str.upper()
    return s.str.extract(r'(?<![A-Z])(IV|III|II|I)(?=[ABC]?(?![A-Z]))')[0]


def truncating(annot):
    """Truncating mutations of a dNdScv annotation: nonsense, essential splice and frameshift indels"""
    a = annot.copy()
    diff = (a.ref.fillna('').str.replace('-', '', regex=False).str.len()
            - a.mut.fillna('').str.replace('-', '', regex=False).str.len())
    a['FRAMESHIFT'] = (a.impact == 'no-SNV') & (diff % 3 != 0)
    a['SUBSTITUTION'] = a.impact.isin(spec.TRUNCATING_SUBSTITUTIONS)
    return a[a.FRAMESHIFT | a.SUBSTITUTION]


class Stratum:
    """dNdScv results of a stratum with the fork's null of truncating substitutions"""

    def __init__(self, work, name):
        d = os.path.join(work, 'dnds', name)
        self.ok = os.path.exists(os.path.join(d, 'genemuts.tsv'))
        if not self.ok:
            return
        self.genemuts = pd.read_csv(os.path.join(d, 'genemuts.tsv'), sep='\t')
        self.sel = pd.read_csv(os.path.join(d, 'sel_cv.tsv'), sep='\t').set_index('gene_name')
        self.annot = pd.read_csv(os.path.join(d, 'annotmuts.tsv.gz'), sep='\t', dtype={'sampleID': str, 'chr': str})
        self.excluded = set(pd.read_csv(os.path.join(d, 'excluded.tsv'), sep='\t', dtype={'SAMPLE': str}).SAMPLE)
        self.drivers = set(self.sel.index[self.sel.qglobal_cv < spec.DRIVER_Q])
        covs = analysis.read_covariates(os.path.join(work, 'dnds', 'covariates_pcs.tsv'))
        self.adjusted, self.fit = background_omega(self.genemuts, exclude=self.drivers, covariates=covs)
        self.layer = MutationLayer(self.adjusted, TRUNC_LAYER)

    def test(self, genes):
        return self.layer.test([g for g in genes if g in self.layer.data])

    def ratio_ci(self, genes, level=0.95):
        """Exact interval of the observed / expected ratio: expected scaled by omega until a tail reaches (1-level)/2"""
        genes = [g for g in genes if g in self.layer.data]
        sub = self.adjusted[self.adjusted.gene_name.isin(genes)]
        obs = int(sum(self.layer.observed(g) for g in genes))
        a = (1 - level) / 2

        def tails(omega):
            scaled = sub.copy()
            scaled[['exp_non', 'exp_spl']] *= omega
            lay = MutationLayer(scaled, TRUNC_LAYER, theta=self.layer.theta)
            pmf = sum_pmfs([lay._gene_null(g)[0] for g in genes if g in lay.data])
            upper = pmf[obs:].sum() if obs < len(pmf) else 0.0          # P(X >= obs)
            lower = pmf[:obs + 1].sum()                                  # P(X <= obs)
            return upper, lower

        def solve(f, lo, hi):
            for _ in range(60):
                mid = np.sqrt(lo * hi)
                if f(mid):
                    hi = mid
                else:
                    lo = mid
            return np.sqrt(lo * hi)

        # lower bound: smallest omega with P(X >= obs) >= a ; upper bound: largest omega with P(X <= obs) >= a
        low = solve(lambda w: tails(w)[0] >= a, 1e-3, 1e3) if obs > 0 else 0.0
        high = solve(lambda w: tails(w)[1] < a, 1e-3, 1e3)
        return float(low), float(high)


def logistic(y, X):
    fit = sm.Logit(y, sm.add_constant(X)).fit(disp=0)
    ci = fit.conf_int().loc['carrier']
    return {'n': int(len(y)), 'odds_ratio': float(np.exp(fit.params.carrier)),
            'ci95': [float(np.exp(ci[0])), float(np.exp(ci[1]))],
            'p_one_sided': float(one_sided(fit.params.carrier, fit.pvalues.carrier)),
            'covariate_odds_ratios': {k: float(np.exp(v)) for k, v in fit.params.items() if k not in ('const', 'carrier')}}


def cox(df, time, event, strata=None):
    d = df[df[time].notna() & df[event].notna() & (df[time] > 0)]
    if d[event].sum() < 5 or d.carrier.sum() < 5:
        return None
    m = PHReg(d[time].values, d[['carrier']].values.astype(float), status=d[event].values.astype(float),
              strata=None if strata is None else d[strata].values).fit()
    lo, hi = np.exp(m.conf_int()[0])
    return {'hazard_ratio': float(np.exp(m.params[0])), 'ci95': [float(lo), float(hi)], 'p_two_sided': float(m.pvalues[0]),
            'tumours': int(len(d)), 'carriers': cell(d.carrier.sum()), 'events': int(d[event].sum()),
            'events_carriers': cell(d[d.carrier == 1][event].sum())}


# ------------------------------------------------------------------------------------------------ analyses
def h1(mss):
    r = mss.test(spec.CORE)
    lo, hi = mss.ratio_ci(spec.CORE)
    present = [g for g in spec.CORE if g in mss.layer.data]
    per_gene = []
    for g in present:
        pmf, mean = mss.layer._gene_null(g)
        per_gene.append({'gene': g, 'observed': int(mss.layer.observed(g)), 'expected': float(mean),
                         'cohort_driver': g in mss.drivers,
                         'dndscv_qtrunc': float(mss.sel.loc[g, 'qtrunc_cv']) if g in mss.sel.index else None})
    return {'observed': int(r['OBSERVED']), 'expected': float(r['EXPECTED']), 'ratio': float(r['RATIO']),
            'ratio_ci95': [lo, hi], 'p_one_sided': float(r['P_VALUE']), 'genes_tested': present,
            'genes_missing': [g for g in spec.CORE if g not in present], 'per_gene': per_gene,
            'background_fit': mss.fit, 'drivers_excluded_from_fit': len(mss.drivers)}


def stage_table(groups, mss):
    """MSS tumours analysed by dNdScv, with stage, carrier status and the H2 covariates"""
    g = groups[(groups.GROUP == 'MSS') & ~groups.index.isin(mss.excluded)].copy()
    missing = [c for c in ['STAGE', 'AGE', 'SEX'] if c not in g.columns]
    if missing:
        raise SystemExit(f'samples.tsv lacks {missing}, required for H2')
    sex = g.SEX.dropna().astype(str).str.upper().str.strip()
    if not set(sex) <= {'MALE', 'FEMALE', 'M', 'F'}:
        raise SystemExit(f'SEX values {sorted(set(sex))[:5]}: map them to MALE / FEMALE in samples.tsv')
    g['STAGE_GROUP'] = stage_group(g.STAGE)
    tr = truncating(mss.annot)
    core = tr[tr.gene.isin(spec.CORE)]
    g['carrier'] = g.index.isin(set(core.sampleID)).astype(float)
    g['carrier_substitution'] = g.index.isin(set(core[core.SUBSTITUTION].sampleID)).astype(float)
    g['IV'] = (g.STAGE_GROUP == 'IV').astype(float)
    g['age_10y'] = pd.to_numeric(g.AGE, errors='coerce') / 10
    g['male'] = g.SEX.astype(str).str.upper().str.strip().isin(['MALE', 'M']).astype(float)
    g.loc[g.SEX.isna(), 'male'] = np.nan
    g['rectum'] = (g.LOCATION == 'rectum').astype(float)
    g['proximal'] = (g.LOCATION == 'proximal').astype(float)
    g.loc[g.LOCATION == 'NA', ['rectum', 'proximal']] = np.nan
    g['log_coding'] = np.log(g.CODING.clip(lower=1))
    purity = pd.to_numeric(g['PURITY'], errors='coerce') if 'PURITY' in g.columns else pd.Series(np.nan, index=g.index)
    source = 'cohort purity' if purity.notna().mean() > 0.9 else 'median allele fraction'
    g['purity'] = purity if source == 'cohort purity' else g.MEDIAN_VAF
    for c in ['OS', 'OS_TIME', 'PFI', 'PFI_TIME']:
        if c in g.columns:
            g[c] = pd.to_numeric(g[c], errors='coerce')
    return g, tr, source


def h2(g, purity_source):
    d = g[g.STAGE_GROUP.notna()]
    cc = d[spec.H2_COVARIATES].notna().all(axis=1)
    model = logistic(d.IV[cc], d.loc[cc, ['carrier'] + spec.H2_COVARIATES])
    tab = pd.crosstab(d.carrier, d.IV).reindex(index=[0.0, 1.0], columns=[0.0, 1.0], fill_value=0)
    fisher = stats.fisher_exact(tab.values, alternative='greater')
    return {**model, 'purity_covariate': purity_source, 'tumours_with_stage': int(len(d)),
            'excluded_incomplete_covariates': int((~cc).sum()),
            'carriers': cell(tab.loc[1.0].sum()), 'carriers_stage_iv': cell(tab.loc[1.0, 1.0]),
            'others': cell(tab.loc[0.0].sum()), 'others_stage_iv': cell(tab.loc[0.0, 1.0]),
            'fraction_stage_iv_carriers': float(tab.loc[1.0, 1.0] / max(tab.loc[1.0].sum(), 1)),
            'fraction_stage_iv_others': float(tab.loc[0.0, 1.0] / max(tab.loc[0.0].sum(), 1)),
            'fisher_odds_ratio': float(fisher[0]), 'fisher_p_one_sided': float(fisher[1])}


def fixed_sequence(p1, p2):
    """H1 then H2, each one-sided at spec.ALPHA; testing stops at the first hypothesis not supported"""
    h1_ok = p1 <= spec.ALPHA
    return {'H1': {'p_one_sided': p1, 'threshold': spec.ALPHA, 'supported': bool(h1_ok)},
            'H2': {'p_one_sided': p2, 'threshold': spec.ALPHA, 'supported': bool(h1_ok and p2 <= spec.ALPHA),
                   'tested': bool(h1_ok)}}


def sets_secondary(mss, sets):
    rows = []
    for s in spec.SETS:
        genes = set(sets[sets.SET == s].SYMBOL) - mss.drivers
        a, b = mss.test(genes), mss.test(genes - set(spec.CORE))
        rows.append({'set': s, 'observed': a['OBSERVED'], 'expected': a['EXPECTED'], 'p': a['P_VALUE'],
                     'observed_without_core': b['OBSERVED'], 'expected_without_core': b['EXPECTED'], 'p_without_core': b['P_VALUE']})
    q = fdr_bh([r['p'] for r in rows])
    for r, qi in zip(rows, q):
        r['q'] = float(qi)
    return rows


def matched_null(g, tr, mss):
    d = g[g.STAGE_GROUP.notna()]
    t = tr[tr.sampleID.isin(d.index)]
    counts = t.groupby('gene').sampleID.nunique()
    pool = counts.drop([x for x in (mss.drivers | set(spec.CORE)) if x in counts.index])
    by_count = {k: list(v.index) for k, v in pool.groupby(pool)}
    carriers_of = t.groupby('gene').sampleID.apply(set)
    iv = set(d.index[d.IV == 1])
    core = set(t[t.gene.isin(spec.CORE)].sampleID)
    observed = len(core & iv) / max(len(core), 1)
    rng = np.random.default_rng(spec.SEED)
    null = []
    for _ in range(spec.NULL_DRAWS):
        pts = set()
        for gene in spec.CORE:
            k = counts.get(gene, 0)
            if k == 0:
                continue
            ks = [c for c in by_count if abs(c - k) <= max(1, k // 5)] or [min(by_count, key=lambda c: abs(c - k))]
            pts |= carriers_of[rng.choice(by_count[rng.choice(ks)])]
        null.append(len(pts & iv) / max(len(pts), 1))
    null = np.array(null)
    return {'observed_fraction_stage_iv': observed, 'null_median': float(np.median(null)),
            'null_95th_percentile': float(np.quantile(null, 0.95)),
            'p_empirical': float((np.sum(null >= observed) + 1) / (len(null) + 1)), 'draws': spec.NULL_DRAWS}


def per_gene_and_sensitivity(g, tr):
    d = g[g.STAGE_GROUP.notna()]
    cc = d[spec.H2_COVARIATES].notna().all(axis=1)
    rows = []
    for gene in spec.CORE:
        pts = set(tr[tr.gene == gene].sampleID) & set(d.index)
        loo = d.copy()
        loo['carrier'] = loo.index.isin(set(tr[tr.gene.isin([x for x in spec.CORE if x != gene])].sampleID)).astype(float)
        m = logistic(loo.IV[cc], loo.loc[cc, ['carrier'] + spec.H2_COVARIATES])
        rows.append({'gene': gene, 'carriers': cell(len(pts)), 'carriers_stage_iv': cell(d.loc[list(pts), 'IV'].sum()),
                     'leave_out_odds_ratio': m['odds_ratio'], 'leave_out_p_one_sided': m['p_one_sided']})
    unadjusted = logistic(d.IV, d[['carrier']])
    no_location = [c for c in spec.H2_COVARIATES if c not in ('rectum', 'proximal')]
    cc2 = d[no_location].notna().all(axis=1)
    without_location = logistic(d.IV[cc2], d.loc[cc2, ['carrier'] + no_location])
    sub = d.assign(carrier=d.carrier_substitution)
    substitution_only = logistic(sub.IV[cc], sub.loc[cc, ['carrier'] + spec.H2_COVARIATES])
    return rows, {'unadjusted': unadjusted, 'without_location': without_location, 'substitution_carriers_only': substitution_only}


def clonality(work, g, mss, tr):
    """Allele fraction relative to the tumour median (clonality proxy) of co-regulator truncations, APC/TP53
    truncations and synonymous mutations"""
    m = pd.read_csv(os.path.join(work, 'coding.tsv.gz'), sep='\t', dtype={'SAMPLE': str, 'CHROM': str})
    if not {'T_ALT', 'T_DEPTH'} <= set(m.columns):
        return None
    m = m[m.SAMPLE.isin(g.index)].copy()
    m['chr'] = m.CHROM.str.replace('^chr', '', regex=True)
    m['RVAF'] = m.T_ALT / m.T_DEPTH.where(m.T_DEPTH > 0) / m.SAMPLE.map(g.MEDIAN_VAF)
    key = ['SAMPLE', 'chr', 'POS']

    def rvaf(annot):
        a = annot.rename(columns={'sampleID': 'SAMPLE', 'pos': 'POS'})[key]
        return a.merge(m[key + ['RVAF']], on=key).drop_duplicates(key).RVAF.dropna()

    groups = {'core_truncating': rvaf(tr[tr.gene.isin(spec.CORE)]), 'apc_tp53_truncating': rvaf(tr[tr.gene.isin(['APC', 'TP53'])]),
              'synonymous': rvaf(mss.annot[mss.annot.impact == 'Synonymous'])}
    out = {k: {'n': int(len(v)), 'median_relative_vaf': float(v.median()) if len(v) else None,
               'clonal_fraction': float((v >= 0.8).mean()) if len(v) else None} for k, v in groups.items()}
    if len(groups['core_truncating']) and len(groups['synonymous']):
        out['p_core_vs_synonymous'] = float(stats.mannwhitneyu(groups['core_truncating'], groups['synonymous'],
                                                               alternative='greater').pvalue)
        out['p_core_vs_apc_tp53'] = float(stats.mannwhitneyu(groups['core_truncating'], groups['apc_tp53_truncating']).pvalue)
    return out


def feasibility(work, out):
    """Blinded feasibility: no observed co-regulator mutation and no carrier-by-stage table is looked at"""
    import power
    groups = pd.read_csv(os.path.join(work, 'groups.tsv'), sep='\t', dtype={'SAMPLE': str}).set_index('SAMPLE')
    mss = Stratum(work, 'MSS')
    expected = float(sum(mss.layer._gene_null(g)[1] for g in spec.CORE if g in mss.layer.data))
    g = groups[(groups.GROUP == 'MSS') & ~groups.index.isin(mss.excluded)]
    stage = stage_group(g.STAGE) if 'STAGE' in g.columns else pd.Series(dtype=str)
    n, p_iv = int(stage.notna().sum()), float((stage == 'IV').sum() / max(stage.notna().sum(), 1))
    power.EXPECTED_PER_TUMOUR = expected / max(len(g), 1)
    res = {'mss_tumours': int(len(g)), 'mss_with_stage': n, 'stage_iv_prevalence': p_iv, 'neutral_expected_h1': expected,
           'power_h1': {str(r): power.power_h1(len(g), r) for r in [1.5, 2, 3]},
           'power_h2': {f'carriers {c:.0%}': {str(o): power.power_h2(n, c, p_iv, o) for o in [1.5, 2.0, 2.8]} for c in [0.06, 0.09]},
           'note': 'power at one-sided alpha {} for the observed size and stage mix; carrier rates assumed'.format(spec.ALPHA)}
    json.dump(clean(res), open(os.path.join(out, 'feasibility.json'), 'w'), indent=1)
    print(json.dumps(clean(res), indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--work', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--cornish', default=None, help='driver table of Cornish et al. (column gene, Gene, SYMBOL or Hugo_Symbol)')
    ap.add_argument('--feasibility', action='store_true',
                    help='power from the cohort size, neutral expectation and stage mix only, before any test (blinded)')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.feasibility:
        return feasibility(a.work, a.out)
    groups = pd.read_csv(os.path.join(a.work, 'groups.tsv'), sep='\t', dtype={'SAMPLE': str}).set_index('SAMPLE')
    sets = pd.read_csv(os.path.join(HERE, 'gene_sets.tsv'), sep='\t')
    mss = Stratum(a.work, 'MSS')
    msi = Stratum(a.work, 'MSI')
    res = {'cohort': {'tumours': int(len(groups)), 'groups': {k: cell(v) for k, v in groups.GROUP.value_counts().items()},
                      'group_source': groups.GROUP_SOURCE.iloc[0], 'mss_excluded_by_dndscv': cell(len(mss.excluded)),
                      'mss_drivers': sorted(mss.drivers)}}

    # primary
    res['H1'] = h1(mss)
    g, tr, purity_source = stage_table(groups, mss)
    res['H2'] = h2(g, purity_source)
    res['decision'] = fixed_sequence(res['H1']['p_one_sided'], res['H2']['p_one_sided'])

    # secondary
    res['S1_reactome_sets'] = sets_secondary(mss, sets)
    res['S2_matched_null'] = matched_null(g, tr, mss)
    res['S3_per_gene'], res['S7_sensitivity'] = per_gene_and_sensitivity(g, tr)
    res['S4_clonality'] = clonality(a.work, g, mss, tr)
    survival = {}
    for end in ['OS', 'PFI']:
        if end in g.columns and f'{end}_TIME' in g.columns:
            d = g[g.STAGE_GROUP.notna()]
            survival[end] = {'unadjusted': cox(d, f'{end}_TIME', end), 'stage_stratified': cox(d, f'{end}_TIME', end, 'STAGE_GROUP'),
                             'stages_i_iii': cox(d[d.STAGE_GROUP != 'IV'], f'{end}_TIME', end, 'STAGE_GROUP')}
    res['S5_survival'] = survival
    if msi.ok:
        r = msi.test(spec.CORE)
        msi_groups = groups[groups.GROUP == 'MSI']
        res['S6_msi'] = {'tumours': cell(len(msi_groups)), 'observed': int(r['OBSERVED']), 'expected': float(r['EXPECTED']),
                         'p_one_sided': float(r['P_VALUE'])}
    if a.cornish:
        t = pd.read_csv(a.cornish, sep=None, engine='python')
        col = next(c for c in t.columns if c in ('gene', 'Gene', 'SYMBOL', 'Hugo_Symbol', 'gene_name', 'Gene name'))
        cg = set(t[col].astype(str).str.strip())
        tcga_less = ['TGIF1', 'PCBP1', 'ZFP36L2', 'ELF3', 'RGMB', 'PRKCI', 'MBD6', 'HNRNPAB', 'ING1', 'RBM10', 'TNRC6B', 'KMT2B', 'BMPR2']
        res['S8_cornish'] = {'drivers_in_table': len(cg), 'core_in_cornish': sorted(set(spec.CORE) & cg),
                             'core_not_in_cornish': sorted(set(spec.CORE) - cg),
                             'tcga_less_established_in_cornish': sorted(set(tcga_less) & cg),
                             'tcga_less_established_not_in_cornish': sorted(set(tcga_less) - cg),
                             'cohort_mss_drivers_in_cornish': f'{len(mss.drivers & cg)} of {len(mss.drivers)}'}

    # provenance
    try:
        commit = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=HERE, capture_output=True, text=True).stdout.strip()
        # uncommitted changes to the registered code make the commit an incomplete record
        dirty = bool(subprocess.run(['git', 'status', '--porcelain', '--', '.', '../../core'], cwd=HERE,
                                    capture_output=True, text=True).stdout.strip())
    except OSError:
        commit, dirty = None, None
    version = os.path.join(a.work, 'dndscv_version.txt')
    res['provenance'] = {'commit': commit, 'uncommitted_changes': dirty, 'spec_sha256': sha256(os.path.join(HERE, 'spec.py')),
                         'gene_sets_sha256': sha256(os.path.join(HERE, 'gene_sets.tsv')),
                         'dndscv': open(version).read().strip() if os.path.exists(version) else None,
                         'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__,
                         'scipy': scipy.__version__, 'statsmodels': statsmodels.__version__,
                         'run_at': datetime.now(timezone.utc).isoformat(timespec='seconds')}
    res = clean(res)
    json.dump(res, open(os.path.join(a.out, 'replication.json'), 'w'), indent=1)
    write_summary(res, os.path.join(a.out, 'summary.md'))
    print(open(os.path.join(a.out, 'summary.md')).read())


def write_summary(r, path):
    h1, h2, dec = r['H1'], r['H2'], r['decision']
    lines = [
        '# Replication of the co-regulator findings', '',
        f"Cohort: {r['cohort']['tumours']} tumours ({r['cohort']['groups']}), groups from {r['cohort']['group_source']}.", '',
        '## Primary hypotheses (fixed sequence H1 then H2, one-sided, family-wise alpha {:.2f})'.format(spec.ALPHA), '',
        '| Hypothesis | Estimate | 95% CI | p (one-sided) | Threshold | Supported |', '|---|---|---|---|---|---|',
        f"| H1 truncating excess in the nine co-regulators (MSS) | {h1['observed']} vs {h1['expected']:.1f} expected, "
        f"ratio {h1['ratio']:.2f} | {h1['ratio_ci95'][0]:.2f}-{h1['ratio_ci95'][1]:.2f} | {h1['p_one_sided']:.3g} | "
        f"{dec['H1']['threshold']:.3f} | {'yes' if dec['H1']['supported'] else 'no'} |",
        f"| H2 carriers more often stage IV (MSS, adjusted) | OR {h2['odds_ratio']:.2f} (n = {h2['n']}) | "
        f"{h2['ci95'][0]:.2f}-{h2['ci95'][1]:.2f} | {h2['p_one_sided']:.3g} | {dec['H2']['threshold']:.3f} | "
        f"{'yes' if dec['H2']['supported'] else ('no' if dec['H2']['tested'] else 'not tested (H1 not supported)')} |", '',
        f"Stage IV: carriers {h2['carriers_stage_iv']} of {h2['carriers']} ({h2['fraction_stage_iv_carriers']:.0%}), "
        f"others {h2['others_stage_iv']} of {h2['others']} ({h2['fraction_stage_iv_others']:.0%}); "
        f"purity covariate: {h2['purity_covariate']}.", '',
        '## Secondary', '',
        '| Reactome set | Observed / expected | q | Without the nine genes | p |', '|---|---|---|---|---|']
    for s in r['S1_reactome_sets']:
        lines.append(f"| {s['set']} | {s['observed']} / {s['expected']:.1f} | {s['q']:.3g} | "
                     f"{s['observed_without_core']} / {s['expected_without_core']:.1f} | {s['p_without_core']:.3g} |")
    n = r['S2_matched_null']
    lines += ['', f"Matched random genes: stage IV fraction among carriers {n['observed_fraction_stage_iv']:.2f}, null median "
              f"{n['null_median']:.2f}, 95th percentile {n['null_95th_percentile']:.2f}, empirical p = {n['p_empirical']:.3g}.", '',
              '| Gene | Truncating obs / exp (substitutions) | Carriers | Stage IV | OR without it | p without it |', '|---|---|---|---|---|---|']
    obs = {x['gene']: x for x in h1['per_gene']}
    for x in r['S3_per_gene']:
        o = obs.get(x['gene'])
        lines.append(f"| {x['gene']} | {o['observed']} / {o['expected']:.2f} |" if o else f"| {x['gene']} | not in RefCDS |")
        lines[-1] += f" {x['carriers']} | {x['carriers_stage_iv']} | {x['leave_out_odds_ratio']:.2f} | {x['leave_out_p_one_sided']:.3g} |"
    s7 = r['S7_sensitivity']
    lines += ['', 'Sensitivity of H2: ' + '; '.join(f"{k.replace('_', ' ')} OR {v['odds_ratio']:.2f} (p = {v['p_one_sided']:.3g})"
                                                    for k, v in s7.items()) + '.']
    if r.get('S4_clonality'):
        c = r['S4_clonality']
        lines += ['', 'Relative allele fraction (median, clonal fraction): ' + '; '.join(
            f"{k.replace('_', ' ')} {v['median_relative_vaf']:.2f}, {v['clonal_fraction']:.0%} (n = {v['n']})"
            for k, v in c.items() if isinstance(v, dict) and v.get('median_relative_vaf') is not None) + '.']
    for end, v in r.get('S5_survival', {}).items():
        parts = [f"{k.replace('_', ' ')} HR {x['hazard_ratio']:.2f} ({x['ci95'][0]:.2f}-{x['ci95'][1]:.2f}), p = {x['p_two_sided']:.3g}"
                 for k, x in v.items() if x]
        lines += ['', f'{end}: ' + ('; '.join(parts) if parts else 'too few events') + '.']
    if r.get('S6_msi'):
        m = r['S6_msi']
        lines += ['', f"MSI tumours ({m['tumours']}): co-regulator truncating {m['observed']} vs {m['expected']:.1f} (p = {m['p_one_sided']:.3g})."]
    if r.get('S8_cornish'):
        c = r['S8_cornish']
        lines += ['', f"Cornish et al. drivers ({c['drivers_in_table']}): co-regulators in their list {c['core_in_cornish']}, "
                  f"not in it {c['core_not_in_cornish']}."]
    p = r['provenance']
    lines += ['', f"Commit {p['commit']}{' (with uncommitted changes)' if p['uncommitted_changes'] else ''}, spec.py sha256 {p['spec_sha256'][:12]}, gene_sets.tsv sha256 {p['gene_sets_sha256'][:12]}, "
              f"dNdScv {p['dndscv']}, run {p['run_at']}."]
    open(path, 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
