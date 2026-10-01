"""Pathway selection and co-occurrence/exclusivity of the fork on TCGA COADREAD strata (ALL, MSS, MSI)"""
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'core'))
from intogen_core.pathways import analysis      # noqa: E402

IMPACT = {'Synonymous': 'synonymous_variant', 'Missense': 'missense_variant', 'Nonsense': 'stop_gained',
          'Essential_Splice': 'splice_donor_variant', 'Stop_loss': 'stop_lost'}
COVARIATES = 'data/covariates_pcs.tsv'
GENE_SETS = 'gene_sets_crc.tsv.gz'


def consequence(impact, ref, mut):
    if impact != 'no-SNV':
        return IMPACT.get(impact, 'other')
    diff = len(ref.replace('-', '')) - len(mut.replace('-', ''))
    if diff == 0:
        return 'other'
    if diff % 3:
        return 'frameshift_variant'
    return 'inframe_deletion' if diff > 0 else 'inframe_insertion'


def prepare(stratum):
    d, o = f'data/{stratum}', f'results/{stratum}'
    os.makedirs(o, exist_ok=True)
    annot = pd.read_csv(f'{d}/annotmuts.tsv.gz', sep='\t', dtype=str)
    annot['Consequence'] = [consequence(i, r, m) for i, r, m in zip(annot['impact'], annot['ref'], annot['mut'])]
    annot['#Uploaded_variation'] = annot['chr'] + '_' + annot['pos'] + '__' + annot['sampleID'] + '__x'
    annot.rename(columns={'gene': 'SYMBOL'})[['#Uploaded_variation', 'Consequence', 'SYMBOL']].to_csv(
        f'{o}/mutations.tsv.gz', sep='\t', index=False)
    sel = pd.read_csv(f'{d}/sel_cv.tsv', sep='\t')
    drivers = sel[sel['qglobal_cv'] < 0.1]
    pd.DataFrame({'SYMBOL': drivers['gene_name'], 'TIER': 1, 'FILTER': 'PASS',
                  'QVALUE_COMBINATION': drivers['qglobal_cv']}).to_csv(f'{o}/vet.tsv', sep='\t', index=False)
    return drivers


def run(stratum):
    o = f'results/{stratum}'
    analysis.run(genemuts=f'data/{stratum}/genemuts.tsv', mutations=f'{o}/mutations.tsv.gz', vet=f'{o}/vet.tsv',
                 gene_sets=GENE_SETS, covariates=COVARIATES, min_size=5, max_events=150,
                 output=f'{o}/pathways.tsv.gz', cooccurrence_output=f'{o}/cooccurrence.tsv.gz',
                 genes_output=f'{o}/genes.tsv.gz')
    st = json.load(open(f'{o}/pathways.tsv.gz.stats.json'))
    df = pd.read_csv(f'{o}/pathways.tsv.gz', sep='\t')
    lt = df[(df.SCOPE == 'long_tail') & (df.LAYER == 'combined')].sort_values('Q_VALUE')
    sig = df[(df.SCOPE == 'long_tail') & (df.LAYER != 'combined') & (df.Q_VALUE < 0.1)]
    print(f'== {stratum}: long-tail sets q<0.1: {(lt.Q_VALUE < 0.1).sum()} | events {st["cooccurrence_events"]}, '
          f'pairs {st["cooccurrence_pairs"]}, co-occurring {st["cooccurrence_significant_pairs"]}, '
          f'exclusive {st["exclusivity_significant_pairs"]}')
    for _, r in lt[lt.Q_VALUE < 0.25].head(15).iterrows():
        best = sig[sig.SET == r.SET].sort_values('Q_VALUE')
        extra = '' if best.empty else f' {best.iloc[0].LAYER} ratio={best.iloc[0].RATIO:.2f} top={best.iloc[0].TOP_GENES}'
        print(f'   q={r.Q_VALUE:.3g} {r.SET} (n={r.N_GENES}){extra}')
    pairs = pd.read_csv(f'{o}/cooccurrence.tsv.gz', sep='\t')
    for label, col in [('exclusive', 'Q_VALUE_EXCLUSIVITY'), ('co-occurring', 'Q_VALUE')]:
        sub = pairs[pairs[col] < 0.1].sort_values(col)
        print(f'   {label}:', '; '.join(f'{a.split(":", 1)[1]}|{b.split(":", 1)[1]} {n}/{e:.1f} q={q:.1g}'
                                       for a, b, n, e, q in zip(sub.EVENT_1, sub.EVENT_2, sub.TUMOURS_BOTH,
                                                                sub.EXPECTED_BOTH, sub[col])))


if __name__ == '__main__':
    for stratum in sys.argv[1:] or ['ALL', 'MSS', 'MSI']:
        drivers = prepare(stratum)
        print(stratum, 'drivers:', ', '.join(drivers.sort_values('qglobal_cv')['gene_name']), flush=True)
        run(stratum)
