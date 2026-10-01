"""
Stage IV and outcome of MSS tumours with a truncating mutation in one of the nine co-regulators: adjusted model,
per gene, a null of carrier sets of matched rare genes, and progression-free interval (TCGA CDR).
Needs mutations_reads.tsv.gz and the CDR columns of clinical.tsv (mutations.R).
"""
import json

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.duration.hazard_regression import PHReg

rng = np.random.default_rng(11)
CORE = ['EP300', 'CREBBP', 'KMT2C', 'KMT2B', 'SIN3A', 'ASXL1', 'ASXL2', 'USP9X', 'TNRC6B']
TRUNC = ['Nonsense_Mutation', 'Splice_Site', 'Frame_Shift_Del', 'Frame_Shift_Ins']
out = {}

groups = pd.read_csv('tumour_groups.tsv', sep='\t').set_index('PATIENT')
mss = groups.index[groups.GROUP == 'MSS']
clin = pd.read_csv('clinical.tsv', sep='\t').drop_duplicates('PATIENT').set_index('PATIENT')
reads = pd.read_csv('mutations_reads.tsv.gz', sep='\t', dtype={'chr': str})
reads['PATIENT'] = reads.SAMPLE.str[:12]
reads = reads[reads.PATIENT.isin(mss)].drop_duplicates(['SAMPLE', 'chr', 'pos', 'gene', 'class'])
reads['VAF'] = reads.t_alt / (reads.t_alt + reads.t_ref)
drivers = set(pd.read_csv('results/MSS/vet.tsv', sep='\t').SYMBOL)

d = groups.loc[mss, ['CODING', 'LOCATION', 'PROJECT']].join(clin[['AGE', 'SEX', 'CDR_STAGE', 'OS', 'OS_TIME', 'PFI', 'PFI_TIME']])
d['STAGE'] = d.CDR_STAGE.str.extract(r'Stage (IV|III|II|I)')[0]
d['IV'] = (d.STAGE == 'IV').astype(float)
d['VAF_MEDIAN'] = reads.groupby('PATIENT').VAF.median()          # tumour purity proxy
trunc = reads[reads['class'].isin(TRUNC)]
carriers = set(trunc[trunc.gene.isin(CORE)].PATIENT)
d['CARRIER'] = d.index.isin(carriers).astype(float)
d = d[d.STAGE.notna()]
tab = pd.crosstab(d.CARRIER, d.IV)
odds, p = stats.fisher_exact(tab.values)
print(tab, f'\nstage IV: carriers {tab.loc[1.0, 1.0]}/{tab.loc[1.0].sum()} vs {tab.loc[0.0, 1.0]}/{tab.loc[0.0].sum()}; '
      f'OR {odds:.2f}, Fisher p = {p:.3g}')
out['fisher'] = {'carriers_iv': int(tab.loc[1.0, 1.0]), 'carriers': int(tab.loc[1.0].sum()),
                 'others_iv': int(tab.loc[0.0, 1.0]), 'others': int(tab.loc[0.0].sum()), 'odds': odds, 'p': p}

# adjusted logistic model
X = pd.DataFrame({'carrier': d.CARRIER, 'age_10y': d.AGE / 10, 'male': (d.SEX == 'MALE').astype(float),
                  'rectal_project': (d.PROJECT == 'READ').astype(float), 'proximal': (d.LOCATION == 'proximal').astype(float),
                  'log_coding': np.log(d.CODING), 'vaf_median': d.VAF_MEDIAN})
ok = X.notna().all(axis=1)
fit = sm.Logit(d.IV[ok], sm.add_constant(X[ok])).fit(disp=0)
ci = fit.conf_int().loc['carrier']
print(f'adjusted (n = {ok.sum()}): carrier OR {np.exp(fit.params.carrier):.2f} [{np.exp(ci[0]):.2f}-{np.exp(ci[1]):.2f}], '
      f'p = {fit.pvalues.carrier:.3g}')
print(np.exp(fit.params).round(2).to_dict())
out['adjusted'] = {'n': int(ok.sum()), 'odds': float(np.exp(fit.params.carrier)), 'ci': [float(np.exp(ci[0])), float(np.exp(ci[1]))],
                   'p': float(fit.pvalues.carrier), 'covariates': list(X.columns[1:])}

# per gene and per definition
per_gene = []
for g in CORE:
    pts = set(trunc[trunc.gene == g].PATIENT) & set(d.index)
    per_gene.append({'gene': g, 'carriers': len(pts), 'stage_iv': int(d.loc[list(pts), 'IV'].sum())})
per_gene = pd.DataFrame(per_gene)
print(per_gene.to_string(index=False))
out['per_gene'] = per_gene.to_dict('records')
subst = set(trunc[trunc.gene.isin(CORE) & trunc['class'].isin(['Nonsense_Mutation', 'Splice_Site'])].PATIENT) & set(d.index)
x = d.loc[list(subst), 'IV']
out['substitution_carriers'] = {'carriers': len(subst), 'stage_iv': int(x.sum()),
                                'p': float(stats.fisher_exact(pd.crosstab(d.index.isin(subst), d.IV).values)[1])}
print('substitution carriers:', out['substitution_carriers'])
loo = {g: stats.fisher_exact(pd.crosstab(d.index.isin(set(trunc[trunc.gene.isin([c for c in CORE if c != g])].PATIENT)), d.IV).values)
       for g in CORE}
out['leave_one_gene_out'] = {g: {'odds': float(v[0]), 'p': float(v[1])} for g, v in loo.items()}
print('leave one gene out p:', {g: round(v[1], 3) for g, v in loo.items()})

# null: nine random non-driver genes with the same number of truncation carriers as each co-regulator
counts = trunc[trunc.PATIENT.isin(d.index)].groupby('gene').PATIENT.nunique()
pool = counts.drop(list((drivers | set(CORE)) & set(counts.index)))
by_count = {k: list(v.index) for k, v in pool.groupby(pool)}
carriers_of = trunc[trunc.PATIENT.isin(d.index)].groupby('gene').PATIENT.apply(set)
iv = set(d.index[d.IV == 1])
obs_frac = tab.loc[1.0, 1.0] / tab.loc[1.0].sum()
null = []
for _ in range(5000):
    pts = set()
    for g in CORE:
        k = counts.get(g, 0)
        ks = [c for c in by_count if abs(c - k) <= max(1, k // 5)]
        pts |= carriers_of[rng.choice(by_count[rng.choice(ks)])]
    null.append(len(pts & iv) / len(pts))
null = np.array(null)
emp = float((np.sum(null >= obs_frac) + 1) / (len(null) + 1))
print(f'stage IV fraction among carriers {obs_frac:.3f}; matched random gene sets: median {np.median(null):.3f}, '
      f'95th percentile {np.quantile(null, 0.95):.3f}, empirical p = {emp:.3g}')
out['matched_null'] = {'observed_fraction': float(obs_frac), 'null_median': float(np.median(null)),
                       'null_95': float(np.quantile(null, 0.95)), 'p': emp, 'draws': len(null)}

# progression-free interval and overall survival (TCGA CDR), with and without stage strata
for end in ['PFI', 'OS']:
    s = d[d[f'{end}_TIME'].notna() & d[end].notna() & (d[f'{end}_TIME'] > 0)]
    res = {}
    for name, strata in [('unadjusted', None), ('stage_stratified', s.STAGE.values)]:
        m = PHReg(s[f'{end}_TIME'].values, s[['CARRIER']].values, status=s[end].values, strata=strata).fit()
        lo, hi = np.exp(m.conf_int()[0])
        res[name] = {'hr': float(np.exp(m.params[0])), 'ci': [float(lo), float(hi)], 'p': float(m.pvalues[0])}
    res['events_carriers'] = int(s[s.CARRIER == 1][end].sum())
    res['n_carriers'] = int((s.CARRIER == 1).sum())
    res['events_others'] = int(s[s.CARRIER == 0][end].sum())
    res['n_others'] = int((s.CARRIER == 0).sum())
    # within stages I-III only
    s3 = s[s.STAGE != 'IV']
    m = PHReg(s3[f'{end}_TIME'].values, s3[['CARRIER']].values, status=s3[end].values, strata=s3.STAGE.values).fit()
    res['stage_i_iii'] = {'hr': float(np.exp(m.params[0])), 'ci': [float(x) for x in np.exp(m.conf_int()[0])], 'p': float(m.pvalues[0]),
                          'n_carriers': int((s3.CARRIER == 1).sum()), 'events_carriers': int(s3[s3.CARRIER == 1][end].sum())}
    out[end] = res
    print(end, json.dumps(res, default=float))

json.dump(out, open('results/coregulators_clinical.json', 'w'), indent=1, default=float)
