"""
Case study, step 1: inputs of the pathway analysis from the dNdScv results of the TCGA cohorts,
the gene-set tests against the neutral model (as first implemented, --no-background) and the
gene sets to be tested with dNdScv genesetdnds (long tails and random sets).
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'core'))
from intogen_core.pathways import analysis                     # noqa: E402
from intogen_core.pathways.genesets import GeneSets, read_gene_sets   # noqa: E402

DATA = 'data'
OUT = 'results'
COHORTS = ['KIRC', 'LUAD', 'BRCA', 'COAD']
IMPACT = {'Synonymous': 'synonymous_variant', 'Missense': 'missense_variant', 'Nonsense': 'stop_gained',
          'Essential_Splice': 'splice_donor_variant', 'Stop_loss': 'stop_lost'}
CGC = set(open(f'{DATA}/cgc81.txt').read().split())


def consequence(row):
    if row['impact'] != 'no-SNV':
        return IMPACT.get(row['impact'], 'other')
    ref, mut = row['ref'].replace('-', ''), row['mut'].replace('-', '')
    diff = len(ref) - len(mut)
    if diff == 0:
        return 'other'                      # dinucleotide / multi-nucleotide substitutions
    if diff % 3:
        return 'frameshift_variant'
    return 'inframe_deletion' if diff > 0 else 'inframe_insertion'


def prepare(cohort):
    d, o = f'{DATA}/{cohort}', f'{OUT}/{cohort}'
    os.makedirs(o, exist_ok=True)
    sel = pd.read_csv(f'{d}/sel_cv.tsv', sep='\t')
    drivers = sel[sel['qglobal_cv'] < 0.1][['gene_name', 'qglobal_cv']]
    annot = pd.read_csv(f'{d}/annotmuts.tsv.gz', sep='\t', dtype=str)
    annot['Consequence'] = annot.apply(consequence, axis=1)
    annot['#Uploaded_variation'] = (annot['chr'] + '_' + annot['pos'] + '__' + annot['sampleID'] + '__' +
                                    annot['ref'] + '__' + annot['mut'])
    annot = annot.rename(columns={'gene': 'SYMBOL'})
    annot[['#Uploaded_variation', 'Consequence', 'SYMBOL']].to_csv(f'{o}/mutations.tsv.gz', sep='\t', index=False)
    # pre-fork IntOGen proxy: genes significant in dNdScv are the drivers (tier 1)
    vet = pd.DataFrame({'SYMBOL': drivers['gene_name'], 'TIER': 1, 'FILTER': 'PASS',
                        'QVALUE_COMBINATION': drivers['qglobal_cv']})
    vet.to_csv(f'{o}/vet.tsv', sep='\t', index=False)
    # sensitivity analysis: the long tail also excludes all the Cancer Gene Census genes
    extra = sorted(CGC - set(vet['SYMBOL']))
    vet_cgc = pd.concat([vet, pd.DataFrame({'SYMBOL': extra, 'TIER': 2, 'FILTER': 'Known cancer gene',
                                            'QVALUE_COMBINATION': 1.0})])
    vet_cgc.to_csv(f'{o}/vet_cgc.tsv', sep='\t', index=False)
    return drivers


def run_fork(cohort, vet_name, tag):
    o = f'{OUT}/{cohort}'
    files = {k: f'{o}/{tag}.{k}.tsv.gz' for k in ['pathways', 'cooccurrence', 'genes']}
    t0 = time.time()
    analysis.run(genemuts=f'{DATA}/{cohort}/genemuts.tsv', mutations=f'{o}/mutations.tsv.gz', vet=f'{o}/{vet_name}',
                 gene_sets=f'{DATA}/gene_sets.tsv.gz', output=files['pathways'],
                 cooccurrence_output=files['cooccurrence'], genes_output=files['genes'], background=False)
    return time.time() - t0


def units(cohort, drivers, rng):
    """Long tails of the gene sets and random gene sets, for genesetdnds"""
    genemuts = pd.read_csv(f'{DATA}/{cohort}/genemuts.tsv', sep='\t')
    universe = set(genemuts['gene_name'])
    sets = GeneSets(read_gene_sets(f'{DATA}/gene_sets.tsv.gz'), universe, min_size=10, max_size=500)
    rows = []
    for s in sets.ids:
        tail = sorted(sets.genes[s] - set(drivers['gene_name']))
        if len(tail) >= 2:
            rows.append((f'LT|{s}', ';'.join(tail)))
    pool = np.array(sorted(universe - set(drivers['gene_name'])))
    sizes = [len(sets.genes[s]) for s in sets.ids]
    random_sets = []
    for i in range(400):
        genes = sorted(rng.choice(pool, int(rng.choice(sizes)), replace=False))
        random_sets.append(genes)
        rows.append((f'RND|{i}', ';'.join(genes)))
    pd.DataFrame(rows, columns=['UNIT', 'GENES']).to_csv(f'{OUT}/{cohort}/gsd_units.tsv', sep='\t', index=False)
    json.dump(random_sets, open(f'{OUT}/{cohort}/random_sets.json', 'w'))
    return len(rows)


if __name__ == '__main__':
    timing = {}
    for cohort in COHORTS:
        drivers = prepare(cohort)
        timing[cohort] = {'drivers': len(drivers)}
        timing[cohort]['fork_seconds'] = run_fork(cohort, 'vet.tsv', 'fork')
        timing[cohort]['fork_cgc_seconds'] = run_fork(cohort, 'vet_cgc.tsv', 'fork_cgc')
        timing[cohort]['gsd_units'] = units(cohort, drivers, np.random.default_rng(1))
        print(cohort, timing[cohort], flush=True)
    json.dump(timing, open(f'{OUT}/timing.json', 'w'), indent=2)
