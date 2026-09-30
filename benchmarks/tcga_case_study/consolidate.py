import json
import numpy as np, pandas as pd
C = ['KIRC', 'LUAD', 'BRCA', 'COAD']
out = {'cohorts': {}}
for c in C:
    info = dict(l.strip().split('\t') for l in open(f'data/{c}/info.tsv'))
    drivers = pd.read_csv(f'results/{c}/vet.tsv', sep='\t')
    gd = pd.read_csv(f'data/{c}/globaldnds.tsv', sep='\t').set_index('name')
    muts = pd.read_csv(f'results/{c}/mutations.tsv.gz', sep='\t', usecols=['Consequence'])
    out['cohorts'][c] = {'tumours': int(info['tumours']), 'theta_dndscv': float(info['theta']),
                         'coding_mutations': int(len(muts)), 'drivers_dndscv_q10': int(len(drivers)),
                         'global_wmis': float(gd.loc['wmis', 'mle']), 'global_wtru': float(gd.loc['wtru', 'mle'])}
# calibration (random sets of non-driver genes), with genesetdnds for all cohorts
cal = {}
for c in C:
    df = pd.read_csv(f'results/{c}/calibration_final.tsv', sep='\t')
    gsd = pd.read_csv(f'results/{c}/gsd.tsv', sep='\t')
    gsd = gsd[gsd.UNIT.str.startswith('RND|')].copy()
    gsd['SET'] = gsd.UNIT.str.split('|').str[1].astype(int)
    gsd['gsd'] = np.where(gsd.wall > 1, gsd.p_wall / 2, 1 - gsd.p_wall / 2)
    df = df.drop(columns=[x for x in ['genesetdnds'] if x in df.columns]).merge(gsd[['SET', 'gsd']], on='SET')
    big = df.N_GENES >= df.N_GENES.quantile(2 / 3)
    cal[c] = {}
    for m, label in [('naive_poisson', 'Naive Poisson'), ('gsd', 'dNdScv genesetdnds'), ('conditional', 'Conditional binomial'),
                     ('shipped_combined', 'Fork as shipped'), ('final_combined', 'Fork corrected')]:
        cal[c][label] = {'p05': float((df[m] < 0.05).mean()), 'p01': float((df[m] < 0.01).mean()),
                         'p05_large': float((df.loc[big, m] < 0.05).mean())}
    cal[c]['n_sets'] = int(len(df))
out['calibration'] = cal
# spike-in power
sp = []
for c in C:
    df = pd.read_csv(f'results/{c}/spikein_final.tsv', sep='\t')
    for (k, w), g in df.groupby(['KIND', 'OMEGA']):
        sp.append({'cohort': c, 'kind': k, 'omega': float(w), 'genes_significant_mean': float(g.GENES_SIGNIFICANT.mean()),
                   'prefork_ora': float(g.PREFORK_ORA.mean()), 'fork_shipped': float(g.FORK.mean()),
                   'fork_corrected': float(g.FINAL.mean()), 'reps': int(len(g))})
out['spikein'] = sp
out['sets'] = {c: {k: {'n': v['n'], 'artefact_like': v['neuronal_muscle_ecm_channel']} for k, v in d.items()}
               for c, d in json.load(open('results/sets_compare.json')).items()}
out['cooccurrence'] = json.load(open('results/cooccurrence.json'))
# final pathway analysis
fin = {}
for c in C:
    st = json.load(open(f'results/{c}/final.pathways.tsv.gz.stats.json'))
    df = pd.read_csv(f'results/{c}/final.pathways.tsv.gz', sep='\t')
    lt = df[(df.SCOPE == 'long_tail') & (df.LAYER == 'combined')].sort_values('Q_VALUE')
    cg = pd.read_csv(f'results/{c}/final_cgc.pathways.tsv.gz', sep='\t')
    ltc = cg[(cg.SCOPE == 'long_tail') & (cg.LAYER == 'combined')]
    sig = df[(df.SCOPE == 'long_tail') & (df.LAYER != 'combined') & (df.Q_VALUE < 0.1)]
    top = []
    for _, r in lt[lt.Q_VALUE < 0.1].iterrows():
        best = sig[sig.SET == r.SET].sort_values('Q_VALUE').iloc[0]
        top.append({'set': r.SET, 'q': float(r.Q_VALUE), 'n_genes': int(r.N_GENES), 'layer': best.LAYER,
                    'ratio': float(best.RATIO), 'top_genes': best.TOP_GENES})
    pairs = pd.read_csv(f'results/{c}/final.cooccurrence.tsv.gz', sep='\t')
    def lst(df):
        return [{'a': a.split(':', 1)[1], 'b': b.split(':', 1)[1], 'both': int(n), 'expected': float(e), 'q': float(q)}
                for a, b, n, e, q in zip(df.EVENT_1, df.EVENT_2, df.TUMOURS_BOTH, df.EXPECTED_BOTH, df.QQ)]
    ex = pairs[pairs.Q_VALUE_EXCLUSIVITY < 0.1].assign(QQ=lambda d: d.Q_VALUE_EXCLUSIVITY)
    co = pairs[pairs.Q_VALUE < 0.1].assign(QQ=lambda d: d.Q_VALUE)
    fin[c] = {'long_tail_sets': int((lt.Q_VALUE < 0.1).sum()), 'long_tail_sets_without_cgc': int((ltc.Q_VALUE < 0.1).sum()),
              'long_tail': top, 'events': st['cooccurrence_events'], 'pairs': st['cooccurrence_pairs'],
              'exclusive': lst(ex), 'cooccurring': lst(co), 'background': st['mutation_background'],
              'theta_used': st['mutation_theta'], 'theta_mle': st['mutation_theta_mle']}
out['final'] = fin
out['timing'] = json.load(open('results/timing.json'))
json.dump(out, open('results/summary.json', 'w'), indent=1, default=float)
for c in C:
    print(c, {k: (round(v['p05'], 3), round(v['p05_large'], 3)) for k, v in cal[c].items() if k != 'n_sets'})
