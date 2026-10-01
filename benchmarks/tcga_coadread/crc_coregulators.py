"""
What the MSS long-tail signal of transcriptional co-regulators is made of: genes, overlap of the selected gene sets,
mutation-level checks (MC3 filters, allele fractions, synonymous and missense rates), co-occurrence and clinical features.
Needs mutations_reads.tsv.gz (MC3 read counts and filters, see mutations.R).
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'core'))
from intogen_core.omics.stats import upper_tail                                        # noqa: E402
from intogen_core.pathways import analysis                                             # noqa: E402
from intogen_core.pathways.selection import MutationLayer, background_omega            # noqa: E402

pd.set_option('display.width', 250)
out = {}
groups = pd.read_csv('tumour_groups.tsv', sep='\t').set_index('PATIENT')
clinical = pd.read_csv('clinical.tsv', sep='\t').drop_duplicates('PATIENT').set_index('PATIENT')
mss = groups.index[groups.GROUP == 'MSS']

# ---------------------------------------------------------------- the fork's null, as in the pathway analysis of MSS
gm = pd.read_csv('data/MSS/genemuts.tsv', sep='\t')
drivers = set(pd.read_csv('results/MSS/vet.tsv', sep='\t').SYMBOL)
covs = analysis.read_covariates('data/covariates_pcs.tsv')
adjusted, _ = background_omega(gm, exclude=drivers, covariates=covs)
trunc = MutationLayer(adjusted, 'mutation_truncating')
mis = MutationLayer(adjusted, 'mutation_missense')
sets = pd.read_csv('gene_sets_crc.tsv.gz', sep='\t')
res = pd.read_csv('results/MSS/pathways.tsv.gz', sep='\t')
lt = res[(res.SCOPE == 'long_tail')]
comb = lt[lt.LAYER == 'combined'].set_index('SET')
sig = comb.index[comb.Q_VALUE < 0.1]
tr = lt[lt.LAYER == 'mutation_truncating'].set_index('SET')
sig_trunc = [s for s in sig if tr.loc[s, 'P_VALUE'] <= lt[(lt.SET == s) & (lt.LAYER != 'combined')].P_VALUE.min()]
print('significant MSS long tails:', len(sig), '| best layer truncating:', len(sig_trunc))


def gene_row(g):
    pmf, mean = trunc._gene_null(g)
    x = trunc.observed(g)
    pm, mm = mis._gene_null(g)
    xs, s, _, es, mu = trunc.data[g]
    return {'GENE': g, 'TRUNC_OBS': int(x), 'TRUNC_EXP': mean, 'TRUNC_P': upper_tail(pmf, x),
            'MIS_OBS': int(mis.observed(g)), 'MIS_EXP': mm, 'MIS_P': upper_tail(pm, mis.observed(g)),
            'SYN_OBS': int(s), 'SYN_EXP': float(mu)}


# ---------------------------------------------------------------- which genes carry the excess
# transcriptional co-regulators among the top contributors (named after inspection: descriptive, not a new test)
CO_REG = {'EP300': 'KAT3 acetyltransferase co-activator', 'CREBBP': 'KAT3 acetyltransferase co-activator',
          'KMT2C': 'enhancer H3K4 methyltransferase (COMPASS-like)', 'KMT2B': 'H3K4 methyltransferase (COMPASS)',
          'SIN3A': 'SIN3-HDAC co-repressor scaffold', 'ASXL1': 'PR-DUB (BAP1 complex) subunit',
          'ASXL2': 'PR-DUB (BAP1 complex) subunit', 'USP9X': 'deubiquitinase', 'TNRC6B': 'miRNA silencing (GW182)'}
core = list(CO_REG)
union = sorted(set(sets[sets.SET.isin(sig_trunc)].SYMBOL) - drivers)
union = [g for g in union if g in trunc.data]
genes = pd.DataFrame([gene_row(g) for g in union])
genes['N_SETS'] = [int(sets[(sets.SYMBOL == g) & sets.SET.isin(sig_trunc)].shape[0]) for g in genes.GENE]
genes['EXCESS'] = genes.TRUNC_OBS - genes.TRUNC_EXP
genes['ROLE'] = genes.GENE.map(CO_REG).fillna('')
genes = genes.sort_values('EXCESS', ascending=False)
total, rest = trunc.test(union), trunc.test([g for g in union if g not in CO_REG])
excess_core = float(genes[genes.GENE.isin(core)].EXCESS.sum())
print(f'\nunion of the {len(sig_trunc)} truncating sets: {len(union)} long-tail genes, truncating {total["OBSERVED"]} vs '
      f'{total["EXPECTED"]:.1f} (p = {total["P_VALUE"]:.2g}); without the co-regulators {rest["OBSERVED"]} vs '
      f'{rest["EXPECTED"]:.1f} (p = {rest["P_VALUE"]:.2g}); co-regulators carry {excess_core:.1f} of '
      f'{total["OBSERVED"] - total["EXPECTED"]:.1f} excess mutations')
print(genes.head(20).round(4).to_string(index=False))
out['union'] = {'sets': len(sig_trunc), 'genes': len(union), 'observed': total['OBSERVED'], 'expected': total['EXPECTED'],
                'p': total['P_VALUE'], 'rest_observed': rest['OBSERVED'], 'rest_expected': rest['EXPECTED'],
                'rest_p': rest['P_VALUE'], 'excess_core': excess_core}
out['genes'] = genes.head(16).to_dict('records')

# ---------------------------------------------------------------- every selected set with and without the co-regulators
rows = []
for s in sig:
    g = set(sets[sets.SET == s].SYMBOL) - drivers
    a, b = trunc.test(g), trunc.test(g - set(core))
    rows.append({'SET': s, 'LAYER': 'truncating' if s in sig_trunc else 'other', 'GENES': len(g), 'CORE_IN_SET': sorted(g & set(core)),
                 'OBS': a['OBSERVED'], 'EXP': a['EXPECTED'], 'P': a['P_VALUE'],
                 'OBS_NO_CORE': b['OBSERVED'], 'EXP_NO_CORE': b['EXPECTED'], 'P_NO_CORE': b['P_VALUE']})
by_set = pd.DataFrame(rows).sort_values('P')
print('\n', by_set.drop(columns='SET').assign(SET=by_set.SET.str[:45]).round(4).to_string(index=False))
out['sets'] = by_set.to_dict('records')

# leave one gene out of the most significant (pre-specified) set
top_set = by_set.SET.iloc[0]
g0 = sorted(set(sets[sets.SET == top_set].SYMBOL) - drivers)
loo = []
for g in [x for x in g0 if trunc.observed(x) > 0]:
    r = trunc.test([x for x in g0 if x != g])
    loo.append({'LEFT_OUT': g, 'OBS': r['OBSERVED'], 'EXP': r['EXPECTED'], 'P': r['P_VALUE']})
loo = pd.DataFrame(loo).sort_values('P', ascending=False)
print(f'\nleave one gene out of {top_set[:50]}:\n', loo.round(6).to_string(index=False))
out['leave_one_out'] = {'set': top_set, 'max_p': float(loo.P.max()), 'worst': loo.iloc[0].LEFT_OUT,
                        'rows': loo.to_dict('records')}
# missense and synonymous of the co-regulators: a higher local mutation rate would raise them too
ct, cm = trunc.test(core), mis.test(core)
cs_o = int(sum(trunc.data[g][1] for g in core))
cs_e = float(sum(trunc.data[g][4] for g in core))
print(f'co-regulators truncating {ct["OBSERVED"]} vs {ct["EXPECTED"]:.1f}; missense {cm["OBSERVED"]} vs {cm["EXPECTED"]:.1f} '
      f'(p = {cm["P_VALUE"]:.2g}); synonymous {cs_o} vs {cs_e:.1f} (Poisson p = {stats.poisson.sf(cs_o - 1, cs_e):.2g})')
out['core_truncating'] = {'observed': ct['OBSERVED'], 'expected': ct['EXPECTED']}
out['core_missense'] = {'observed': cm['OBSERVED'], 'expected': cm['EXPECTED'], 'p': cm['P_VALUE']}
out['core_synonymous'] = {'observed': cs_o, 'expected': cs_e, 'p': float(stats.poisson.sf(cs_o - 1, cs_e))}
sel_cv = pd.read_csv('data/MSS/sel_cv.tsv', sep='\t').set_index('gene_name')
print(sel_cv.loc[core, ['n_non', 'n_spl', 'wnon_cv', 'ptrunc_cv', 'qtrunc_cv', 'qglobal_cv']].round(4).to_string())
out['dndscv_genes'] = sel_cv.loc[core, ['ptrunc_cv', 'qtrunc_cv', 'qglobal_cv']].reset_index().to_dict('records')

# ---------------------------------------------------------------- the mutations: filters and allele fractions
reads = pd.read_csv('mutations_reads.tsv.gz', sep='\t', dtype={'chr': str})
reads['PATIENT'] = reads.SAMPLE.str[:12]
first = pd.read_csv('coadread_mutations.tsv.gz', sep='\t', usecols=['SAMPLE']).SAMPLE.unique()
reads = reads[reads.SAMPLE.isin(first) & reads.PATIENT.isin(mss)].drop_duplicates(['SAMPLE', 'chr', 'pos', 'gene', 'class']).copy()
reads['VAF'] = reads.t_alt / (reads.t_alt + reads.t_ref)
med = reads[reads['filter'] == 'PASS'].groupby('PATIENT').VAF.median()
reads['RVAF'] = reads.VAF / reads.PATIENT.map(med)
TRUNC = ['Nonsense_Mutation', 'Splice_Site', 'Frame_Shift_Del', 'Frame_Shift_Ins']
cr = reads[reads.gene.isin(core) & reads['class'].isin(TRUNC)].copy()
cr = cr.join(groups[['CODING', 'INDEL_FRACTION', 'LOCATION']], on='PATIENT')
print('\n', cr[['PATIENT', 'gene', 'class', 'protein', 'exon', 'VAF', 'RVAF', 'filter', 'CODING', 'LOCATION']]
      .sort_values(['gene', 'protein']).round(2).to_string(index=False))


def describe(df):
    return {'n': int(len(df)), 'pass': float((df['filter'] == 'PASS').mean()),
            'wga': float(df['filter'].str.contains('wga', na=False).mean()),
            'oxog': float(df['filter'].str.contains('oxog', na=False).mean()),
            'rvaf_median': float(df.RVAF.median()), 'clonal': float((df.RVAF >= 0.8).mean())}


comparisons = {
    'core truncating': cr,
    'core truncating substitutions': cr[cr['class'].isin(['Nonsense_Mutation', 'Splice_Site'])],
    'driver truncating (APC, TP53)': reads[reads.gene.isin(['APC', 'TP53']) & reads['class'].isin(TRUNC)],
    'driver hotspots (KRAS, PIK3CA, BRAF)': reads[reads.gene.isin(['KRAS', 'PIK3CA', 'BRAF']) & (reads['class'] == 'Missense_Mutation')],
    'other long-tail truncating': reads[~reads.gene.isin(set(core) | drivers) & reads['class'].isin(TRUNC)],
    'synonymous': reads[reads['class'] == 'Silent'],
}
cmp_rows = {k: describe(v) for k, v in comparisons.items()}
print('\n', pd.DataFrame(cmp_rows).T.round(3).to_string())
u = stats.mannwhitneyu(cr.RVAF.dropna(), comparisons['synonymous'].RVAF.dropna())
v = stats.mannwhitneyu(cr.RVAF.dropna(), comparisons['driver truncating (APC, TP53)'].RVAF.dropna())
print(f'relative VAF core vs synonymous p = {u.pvalue:.2g}; core vs APC/TP53 truncating p = {v.pvalue:.2g}')
out['mutations'] = {'table': cr[['PATIENT', 'gene', 'class', 'protein', 'VAF', 'RVAF', 'filter', 'CODING', 'LOCATION']]
                    .sort_values(['gene', 'protein']).to_dict('records'),
                    'comparison': cmp_rows, 'p_rvaf_vs_synonymous': float(u.pvalue), 'p_rvaf_vs_apc_tp53': float(v.pvalue)}
# substitution spectrum of the core nonsense mutations (G>T / C>A would point to oxidative artefacts)
PYRIMIDINE = {'G>A': 'C>T', 'G>T': 'C>A', 'G>C': 'C>G', 'T>A': 'A>T', 'T>C': 'A>G', 'T>G': 'A>C'}
annot = pd.read_csv('data/MSS/annotmuts.tsv.gz', sep='\t', dtype={'chr': str})
ann = annot[annot.gene.isin(core) & (annot.impact == 'Nonsense')].copy()
ann['SUB'] = (ann.ref + '>' + ann.mut).replace(PYRIMIDINE)
spec = ann.SUB.value_counts()
bg = annot[annot.impact == 'Nonsense']
bspec = (bg.ref + '>' + bg.mut).replace(PYRIMIDINE).value_counts(normalize=True)
# arginine CGA codons to TGA (C>T at CpG, the clock-like process)
ann = ann.merge(cr[['PATIENT', 'chr', 'pos', 'protein']].astype({'pos': int}), left_on=['sampleID', 'chr', 'pos'],
                right_on=['PATIENT', 'chr', 'pos'])
cga = int((ann.protein.str.startswith('p.R') & (ann.SUB == 'C>T')).sum())
print('\ncore nonsense spectrum:', spec.to_dict(), f'(CGA>TGA {cga}) | all MSS nonsense:', bspec.round(3).to_dict())
out['spectrum'] = {'core': spec.to_dict(), 'background_fraction': bspec.to_dict(), 'cga_tga': cga}

# ---------------------------------------------------------------- the tumours
carriers = sorted(set(cr.PATIENT))
sub_carriers = sorted(set(cr[cr['class'].isin(['Nonsense_Mutation', 'Splice_Site'])].PATIENT))
g = groups.loc[mss].join(clinical[['STAGE', 'AGE', 'SEX']])
g['CARRIER'] = g.index.isin(carriers)
per_tumour = cr.groupby('PATIENT').gene.nunique()
# calls within 10 bp in the same gene and tumour are one event (e.g. adjacent splice-site changes)
n_events = int(cr.sort_values('pos').groupby(['PATIENT', 'gene']).pos.apply(lambda x: (x.diff().fillna(99) > 10).sum()).sum())
stage = g.STAGE.str.extract(r'Stage (IV|III|II|I)')[0]
# tumours with two or more co-regulator genes truncated if the genes were hit independently
f = cr.groupby('gene').PATIENT.nunique() / len(mss)
p_none = np.prod(1 - f.values)
p_one = sum(fi * p_none / (1 - fi) for fi in f.values)
expected_two = float(len(mss) * (1 - p_none - p_one))
usp9x = cr[cr.gene == 'USP9X'].join(clinical[['SEX']], on='PATIENT')[['PATIENT', 'VAF', 'SEX']]
print('USP9X carriers:', usp9x.round(2).to_dict('records'))
stage4 = stats.fisher_exact(pd.crosstab(g.CARRIER, stage == 'IV').values)
print(f'\ncarriers: {len(carriers)} tumours ({len(sub_carriers)} with substitutions), {len(cr)} calls, {n_events} events; '
      f'tumours with 2+ core genes hit: {(per_tumour > 1).sum()} (expected if independent {expected_two:.1f}); '
      f'stage IV odds ratio {stage4[0]:.2f} (Fisher p = {stage4[1]:.2g})')
loc = pd.crosstab(g.CARRIER, g.LOCATION)
stg = pd.crosstab(g.CARRIER, stage)
print(loc, '\n', stg)
burden_p = stats.mannwhitneyu(g[g.CARRIER].CODING, g[~g.CARRIER].CODING).pvalue
age_p = stats.mannwhitneyu(g[g.CARRIER].AGE.dropna(), g[~g.CARRIER].AGE.dropna()).pvalue
loc_p = stats.chi2_contingency(loc.values)[1]
stage_p = stats.chi2_contingency(stg.values)[1]
print(f'burden {g[g.CARRIER].CODING.median():.0f} vs {g[~g.CARRIER].CODING.median():.0f} (p = {burden_p:.2g}); '
      f'age {g[g.CARRIER].AGE.median():.0f} vs {g[~g.CARRIER].AGE.median():.0f} (p = {age_p:.2g}); '
      f'location p = {loc_p:.2g}; stage p = {stage_p:.2g}')
out['tumours'] = {'carriers': len(carriers), 'carriers_substitution': len(sub_carriers), 'mutations': int(len(cr)), 'events': n_events,
                  'two_or_more_genes': int((per_tumour > 1).sum()), 'two_or_more_expected': expected_two,
                  'stage_iv_odds': float(stage4[0]), 'p_stage_iv': float(stage4[1]), 'usp9x': usp9x.to_dict('records'),
                  'frameshift_fraction': float(cr['class'].str.startswith('Frame').mean()),
                  'location': {str(k): {str(c): int(v) for c, v in r.items()} for k, r in loc.iterrows()},
                  'stage': {str(k): {str(c): int(v) for c, v in r.items()} for k, r in stg.iterrows()},
                  'burden_carriers': float(g[g.CARRIER].CODING.median()), 'burden_others': float(g[~g.CARRIER].CODING.median()),
                  'p_burden': float(burden_p), 'p_age': float(age_p), 'p_location': float(loc_p), 'p_stage': float(stage_p)}

# co-occurrence and exclusivity of the co-regulator events with the drivers (exhaustive MSS pairs)
pairs = pd.read_csv('results/MSS/network_pairs.tsv.gz', sep='\t')
events = [f'mutation_truncating:{s}' for s in sig_trunc]
sel = pairs[pairs.EVENT_1.isin(events) | pairs.EVENT_2.isin(events)].copy()
sel['OTHER'] = np.where(sel.EVENT_1.isin(events), sel.EVENT_2, sel.EVENT_1)
sel = sel[~sel.OTHER.isin(events) & ~sel.OTHER.str.contains('REACTOME')]
sel['SET'] = np.where(sel.EVENT_1.isin(events), sel.EVENT_1, sel.EVENT_2)
best = sel.loc[sel.groupby('OTHER').P_VALUE.idxmin()].sort_values('P_VALUE')
bestx = sel.loc[sel.groupby('OTHER').P_VALUE_EXCLUSIVITY.idxmin()].sort_values('P_VALUE_EXCLUSIVITY')
cols = ['OTHER', 'SET', 'TUMOURS_1', 'TUMOURS_2', 'TUMOURS_BOTH', 'EXPECTED_BOTH']
print('\nco-occurrence with co-regulator events:\n', best[cols + ['P_VALUE', 'Q_VALUE']].head(8).round(4).to_string(index=False))
print('exclusivity:\n', bestx[cols + ['P_VALUE_EXCLUSIVITY', 'Q_VALUE_EXCLUSIVITY']].head(6).round(4).to_string(index=False))
out['cooccurrence'] = best[cols + ['P_VALUE', 'Q_VALUE']].head(6).to_dict('records')
out['exclusivity'] = bestx[cols + ['P_VALUE_EXCLUSIVITY', 'Q_VALUE_EXCLUSIVITY']].head(4).to_dict('records')

# the same genes in MSI tumours and in all tumours
for stratum in ['MSI', 'ALL']:
    gm2 = pd.read_csv(f'data/{stratum}/genemuts.tsv', sep='\t')
    drv2 = set(pd.read_csv(f'results/{stratum}/vet.tsv', sep='\t').SYMBOL)
    adj2, _ = background_omega(gm2, exclude=drv2, covariates=covs)
    r = MutationLayer(adj2, 'mutation_truncating').test([x for x in core if x not in drv2])
    print(f'{stratum}: core truncating {r["OBSERVED"]} vs {r["EXPECTED"]:.1f} (p = {r["P_VALUE"]:.2g}); drivers among core: {sorted(drv2 & set(core))}')
    out[f'core_{stratum}'] = {'observed': r['OBSERVED'], 'expected': r['EXPECTED'], 'p': r['P_VALUE'], 'drivers': sorted(drv2 & set(core))}


def clean(x):
    """JSON-safe copy: numpy scalars to Python, NaN to null"""
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, (np.integer, np.bool_)):
        return x.item()
    if isinstance(x, (float, np.floating)):
        return None if np.isnan(x) else float(x)
    return x


json.dump(clean(out), open('results/coregulators.json', 'w'), indent=1)
