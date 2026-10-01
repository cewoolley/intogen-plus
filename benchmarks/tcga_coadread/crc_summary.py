"""Consolidated results of the TCGA COADREAD analysis for the report"""
import json
import pandas as pd

S = {}
groups = pd.read_csv('tumour_groups.tsv', sep='\t')
S['cohort'] = {'patients': int(len(groups)), 'groups': groups.GROUP.value_counts().to_dict(),
               'project': groups.PROJECT.value_counts().to_dict(),
               'msi_location': groups[groups.GROUP == 'MSI'].LOCATION.value_counts().to_dict(),
               'median_coding': groups.groupby('GROUP').CODING.median().to_dict(),
               'indel_fraction': groups.groupby('GROUP').INDEL_FRACTION.median().to_dict(),
               'ihc_vs_group': pd.crosstab(groups.GROUP, groups.MMR_IHC_LOSS.fillna('NA')).to_dict()}
drivers = {}
for s in ['ALL', 'MSS', 'MSI']:
    sel = pd.read_csv(f'data/{s}/sel_cv.tsv', sep='\t')
    d = sel[sel.qglobal_cv < 0.1].sort_values('qglobal_cv')
    drivers[s] = [{'gene': g, 'q': float(q), 'n_mis': int(a), 'n_trunc': int(b + c), 'n_ind': int(e)}
                  for g, q, a, b, c, e in zip(d.gene_name, d.qglobal_cv, d.n_mis, d.n_non, d.n_spl, d.n_ind)]
S['drivers'] = drivers
lt = {}
for s in ['ALL', 'MSS', 'MSI']:
    df = pd.read_csv(f'results/{s}/pathways.tsv.gz', sep='\t')
    comb = df[(df.SCOPE == 'long_tail') & (df.LAYER == 'combined')]
    sig = df[(df.SCOPE == 'long_tail') & (df.LAYER != 'combined')]
    rows = []
    for _, r in comb[comb.Q_VALUE < 0.1].sort_values('Q_VALUE').iterrows():
        b = sig[sig.SET == r.SET].sort_values('P_VALUE').iloc[0]
        rows.append({'set': r.SET, 'q': float(r.Q_VALUE), 'n_genes': int(r.N_GENES), 'layer': b.LAYER,
                     'observed': int(b.OBSERVED), 'expected': float(b.EXPECTED), 'ratio': float(b.RATIO), 'top': b.TOP_GENES})
    lt[s] = rows
S['long_tail'] = lt
pairs = {}
for s in ['MSS', 'MSI']:
    p = pd.read_csv(f'results/{s}/network_pairs.tsv.gz', sep='\t')
    pairs[s] = {'events': int(len(set(p.EVENT_1) | set(p.EVENT_2))), 'pairs': int(len(p)),
                'exclusive': p[p.Q_VALUE_EXCLUSIVITY < 0.1][['A', 'B', 'TUMOURS_BOTH', 'EXPECTED_BOTH', 'Q_VALUE_EXCLUSIVITY', 'KNOWN']].to_dict('records'),
                'cooccurring': p[p.Q_VALUE < 0.1][['A', 'B', 'TUMOURS_BOTH', 'EXPECTED_BOTH', 'Q_VALUE', 'KNOWN']].to_dict('records')}
for s in ['ALL']:
    p = pd.read_csv(f'results/{s}/cooccurrence.tsv.gz', sep='\t')
    p['A'] = p.EVENT_1.str.split(':', n=1).str[1]
    p['B'] = p.EVENT_2.str.split(':', n=1).str[1]
    pairs[s] = {'pairs': int(len(p)),
                'exclusive': p[p.Q_VALUE_EXCLUSIVITY < 0.1][['A', 'B', 'TUMOURS_BOTH', 'EXPECTED_BOTH', 'Q_VALUE_EXCLUSIVITY']].to_dict('records'),
                'cooccurring': p[p.Q_VALUE < 0.1][['A', 'B', 'TUMOURS_BOTH', 'EXPECTED_BOTH', 'Q_VALUE']].to_dict('records')}
S['pairs'] = pairs
S['robustness'] = json.load(open('results/robustness.json'))
S['immune'] = json.load(open('results/immune_escape.json'))
assoc = pd.read_csv('results/site_msi_associations.tsv', sep='\t')
S['location'] = assoc[(assoc.TEST == 'location (MSS)')].sort_values('P')[['GENE', 'MUTATED', 'FREQ_PROXIMAL', 'FREQ_DISTAL', 'FREQ_RECTUM', 'P', 'Q']].to_dict('records')
S['msi_vs_mss'] = assoc[(assoc.TEST == 'MSI vs MSS')].sort_values('P')[['GENE', 'FREQ_MSS', 'FREQ_MSI', 'ODDS', 'Q']].to_dict('records')
S['power'] = pd.read_csv('results/power.tsv', sep='\t').to_dict('records')
json.dump(S, open('results/crc_summary.json', 'w'), indent=1, default=float)
print(json.dumps({k: v for k, v in S['cohort'].items() if k != 'ihc_vs_group'}, default=float))
print({s: len(v) for s, v in drivers.items()}, {s: len(v) for s, v in lt.items()})
print(S['pairs']['ALL']['exclusive'])
