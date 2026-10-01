"""Robustness of the CRC candidates: location-stratified exclusivity, composition of the co-regulator long tail"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'core'))
from intogen_core.omics.stats import poisson_binomial_tail       # noqa: E402
from intogen_core.pathways import analysis, cooccurrence          # noqa: E402
from intogen_core.pathways.selection import MUTATION_LAYERS, MutationLayer, background_omega   # noqa: E402
from intogen_core.omics.stats import simes_combine               # noqa: E402

groups = pd.read_csv('tumour_groups.tsv', sep='\t').set_index('PATIENT')
out = {}


def stratified_test(layer, g1, g2, strata):
    """Burden-elastic null fitted within each stratum (e.g. location); both tails"""
    h1, _ = layer.event([g1]) if isinstance(g1, str) else layer.event(g1)
    h2, _ = layer.event([g2]) if isinstance(g2, str) else layer.event(g2)
    q = np.zeros(len(layer.tumours))
    burden = layer.burden - layer.hits[[layer.gene_index[g] for g in ([g1] if isinstance(g1, str) else g1) if g in layer.gene_index]].sum(axis=0) \
        - layer.hits[[layer.gene_index[g] for g in ([g2] if isinstance(g2, str) else g2) if g in layer.gene_index]].sum(axis=0)
    for s in np.unique(strata):
        idx = strata == s
        p1 = cooccurrence.event_probabilities(h1[idx], np.log1p(burden[idx]))
        p2 = cooccurrence.event_probabilities(h2[idx], np.log1p(burden[idx]))
        q[idx] = p1 * p2
    tail = np.append(poisson_binomial_tail(q), 0.0)
    both = int((h1 & h2).sum())
    return {'both': both, 'expected': float(q.sum()), 'p_cooccurrence': float(tail[both]),
            'p_exclusivity': float(1 - tail[both + 1]), 'n1': int(h1.sum()), 'n2': int(h2.sum())}


tumours, hits = analysis.mutation_hits('results/MSS/mutations.tsv.gz')
layer = cooccurrence.Layer('mutation', hits['mutation'][['SYMBOL', 'SAMPLE']], tumours)
loc = groups.reindex(tumours)['LOCATION'].fillna('NA').values
side = np.where(loc == 'proximal', 'proximal', np.where(loc == 'NA', 'NA', 'distal_rectum'))
for g1, g2 in [('PIK3CA', 'TP53'), ('KRAS', 'NRAS'), ('KRAS', 'BRAF'), ('KRAS', 'TP53'), ('APC', 'TP53'), ('APC', 'KRAS')]:
    res = {'unstratified': stratified_test(layer, g1, g2, np.zeros(len(tumours))),
           'by_location': stratified_test(layer, g1, g2, loc),
           'by_side': stratified_test(layer, g1, g2, side)}
    out[f'MSS {g1}-{g2}'] = res
    print(f'MSS {g1}-{g2}:', {k: f"{v['both']} vs {v['expected']:.1f} p_ex={v['p_exclusivity']:.2g} p_co={v['p_cooccurrence']:.2g}"
                               for k, v in res.items()})

# strongest sub-threshold pairs in MSS (hypotheses), re-tested with a location-stratified null
pairs = pd.read_csv('results/MSS/network_pairs.tsv.gz', sep='\t')
top = []
for direction, col in [('co-occurrence', 'P_VALUE'), ('exclusivity', 'P_VALUE_EXCLUSIVITY')]:
    for _, r in pairs.sort_values(col).head(8).iterrows():
        e1 = r.EVENT_1.split(':', 1)[1]
        e2 = r.EVENT_2.split(':', 1)[1]
        if 'LONG_TAIL' in e1 or 'LONG_TAIL' in e2 or 'REACTOME' in e1 or 'REACTOME' in e2:
            continue
        s = stratified_test(layer, e1, e2, loc)
        top.append({'pair': f'{e1}-{e2}', 'direction': direction, 'both': int(r.TUMOURS_BOTH),
                    'expected': float(r.EXPECTED_BOTH), 'p': float(r[col]),
                    'q': float(r['Q_VALUE' if direction == 'co-occurrence' else 'Q_VALUE_EXCLUSIVITY']),
                    'expected_location': s['expected'],
                    'p_location': s['p_cooccurrence' if direction == 'co-occurrence' else 'p_exclusivity'],
                    'known': r.KNOWN})
out['top_mss_pairs'] = top
print(pd.DataFrame(top).round(4).to_string(index=False))

# co-regulator long tail in MSS: who carries the truncating mutations
CO_REG = ['EP300', 'CREBBP', 'KMT2C', 'KMT2B', 'SIN3A', 'TNRC6B', 'HDAC3', 'ASXL1', 'ASXL2', 'USP9X', 'KMT2D', 'NCOR1']
m = pd.read_csv('coadread_mutations.tsv.gz', sep='\t', dtype=str)
trunc = m[m.gene.isin(CO_REG) & m['class'].isin(['Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'Splice_Site'])
          & m.PATIENT.isin(groups.index[groups.GROUP == 'MSS'])].copy()
trunc = trunc.join(groups[['CODING', 'INDEL_FRACTION', 'LOCATION']], on='PATIENT')
print(trunc[['PATIENT', 'gene', 'class', 'protein', 'CODING', 'INDEL_FRACTION', 'LOCATION']].sort_values(['gene', 'CODING']).to_string(index=False))
out['coregulator_truncating'] = {'tumours': int(trunc.PATIENT.nunique()), 'mutations': int(len(trunc)),
                                 'frameshift_fraction': float(trunc['class'].str.startswith('Frame').mean()),
                                 'median_coding_mutations': float(trunc.CODING.median()),
                                 'high_indel_tumours': int(trunc[trunc.INDEL_FRACTION > 0.08].PATIENT.nunique())}
print(out['coregulator_truncating'])

# long-tail test of the co-regulator signal without the MSS tumours that look MSI-like (high indel fraction)
gm = pd.read_csv('data/MSS/genemuts.tsv', sep='\t')
covs = analysis.read_covariates('data/covariates_pcs.tsv')
drivers = set(pd.read_csv('results/MSS/vet.tsv', sep='\t').SYMBOL)
adjusted, _ = background_omega(gm, exclude=drivers, covariates=covs)
layers = {n: MutationLayer(adjusted, n) for n in MUTATION_LAYERS}
sets = pd.read_csv('gene_sets_crc.tsv.gz', sep='\t')
for s in ['REACTOME_RUNX1_REGULATES_GENES_INVOLVED_IN_MEGAKARYOCYTE_DIFFERENTIATION_AND_PLATELET_FUNCTION',
          'REACTOME_STAT3_NUCLEAR_EVENTS_DOWNSTREAM_OF_ALK_SIGNALING', 'REACTOME_SUMOYLATION_OF_TRANSCRIPTION_COFACTORS',
          'CRC_RTK_RAS', 'CRC_TGFB_BMP', 'CRC_SWI_SNF', 'CRC_IMMUNE_ESCAPE']:
    genes = set(sets[sets.SET == s].SYMBOL) - drivers
    r = {n: L.test(genes) for n, L in layers.items()}
    print(f'{s[:60]:60s} trunc obs/exp {r["mutation_truncating"]["OBSERVED"]}/{r["mutation_truncating"]["EXPECTED"]:.1f} '
          f'p={r["mutation_truncating"]["P_VALUE"]:.2g} | all {r["mutation"]["OBSERVED"]}/{r["mutation"]["EXPECTED"]:.1f} '
          f'p={r["mutation"]["P_VALUE"]:.2g}')
json.dump(out, open('results/robustness.json', 'w'), indent=1, default=float)
