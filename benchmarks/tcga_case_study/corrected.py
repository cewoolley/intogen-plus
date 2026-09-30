"""Case study with the final code of the fork (covariate-matched background, burden-elastic co-occurrence)"""
import json, sys, time
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'core'))
from intogen_core.pathways import analysis
from intogen_core.pathways.selection import background_omega

COVS = 'data/covariates_pcs.tsv'
for cohort in ['KIRC', 'LUAD', 'BRCA', 'COAD']:
    o = f'results/{cohort}'
    t0 = time.time()
    analysis.run(genemuts=f'data/{cohort}/genemuts.tsv', mutations=f'{o}/mutations.tsv.gz', vet=f'{o}/vet.tsv',
                 gene_sets='data/gene_sets.tsv.gz', output=f'{o}/final.pathways.tsv.gz', covariates=COVS,
                 cooccurrence_output=f'{o}/final.cooccurrence.tsv.gz', genes_output=f'{o}/final.genes.tsv.gz')
    seconds = time.time() - t0
    gm = pd.read_csv(f'data/{cohort}/genemuts.tsv', sep='\t')
    drivers = set(pd.read_csv(f'{o}/vet.tsv', sep='\t')['SYMBOL'])
    adjusted, info = background_omega(gm, exclude=drivers, covariates=analysis.read_covariates(COVS))
    adjusted.to_csv(f'{o}/genemuts_final.tsv', sep='\t', index=False)
    st = json.load(open(f'{o}/final.pathways.tsv.gz.stats.json'))
    df = pd.read_csv(f'{o}/final.pathways.tsv.gz', sep='\t')
    lt = df[(df.SCOPE == 'long_tail') & (df.LAYER == 'combined')].sort_values('Q_VALUE')
    pairs = pd.read_csv(f'{o}/final.cooccurrence.tsv.gz', sep='\t')
    print(f'== {cohort} ({seconds:.0f} s) background: {st["mutation_background"]["mutation_missense"]["omega_quantiles"]}')
    print('  long tail q<0.1:', (lt.Q_VALUE < 0.1).sum(), '|', '; '.join(f'{s}({q:.1e})' for s, q in zip(lt.SET.head(8), lt.Q_VALUE.head(8))))
    print('  events', st['cooccurrence_events'], 'pairs', st['cooccurrence_pairs'], 'co-occurring', st['cooccurrence_significant_pairs'],
          'exclusive', st['exclusivity_significant_pairs'])
    ex = pairs[pairs.Q_VALUE_EXCLUSIVITY < 0.1].head(12)
    print('  exclusive:', '; '.join(f'{a.split(":")[1][:28]}|{b.split(":")[1][:28]} {n}/{e:.1f}' for a, b, n, e in zip(ex.EVENT_1, ex.EVENT_2, ex.TUMOURS_BOTH, ex.EXPECTED_BOTH)))
    co = pairs[pairs.Q_VALUE < 0.1].head(12)
    print('  co-occurring:', '; '.join(f'{a.split(":")[1][:28]}|{b.split(":")[1][:28]} {n}/{e:.1f}' for a, b, n, e in zip(co.EVENT_1, co.EVENT_2, co.TUMOURS_BOTH, co.EXPECTED_BOTH)), flush=True)

# sensitivity: the long tail also excludes all the Cancer Gene Census genes
for cohort in ['KIRC', 'LUAD', 'BRCA', 'COAD']:
    o = f'results/{cohort}'
    analysis.run(genemuts=f'data/{cohort}/genemuts.tsv', mutations=f'{o}/mutations.tsv.gz', vet=f'{o}/vet_cgc.tsv',
                 gene_sets='data/gene_sets.tsv.gz', output=f'{o}/final_cgc.pathways.tsv.gz', covariates=COVS,
                 cooccurrence_output=f'{o}/final_cgc.cooccurrence.tsv.gz', genes_output=f'{o}/final_cgc.genes.tsv.gz')
    df = pd.read_csv(f'{o}/final_cgc.pathways.tsv.gz', sep='\t')
    lt = df[(df.SCOPE == 'long_tail') & (df.LAYER == 'combined')].sort_values('Q_VALUE')
    print(f'== {cohort} long tail without CGC genes q<0.1:', (lt.Q_VALUE < 0.1).sum(), '|',
          '; '.join(f'{s}({q:.1e})' for s, q in zip(lt.SET.head(5), lt.Q_VALUE.head(5))), flush=True)
