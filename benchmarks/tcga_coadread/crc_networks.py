"""
Mining of driver networks in TCGA COADREAD with the fork's burden-aware co-occurrence / exclusivity test.

Events tested in each stratum (MSS, MSI):
- every gene that is a driver (dNdScv q < 0.1) in any stratum, altered in >= 3 tumours of the stratum
- the long tails (genes other than the stratum's drivers) of the curated CRC pathways
- the long tails of gene sets selected by the pathway analysis of the stratum (layer where most enriched)
Associations with the anatomical site (MSS) and with MSI status are tested for every driver.
"""
import itertools
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats as sps

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'core'))
from intogen_core.omics.stats import fdr_bh                            # noqa: E402
from intogen_core.pathways import analysis, cooccurrence              # noqa: E402
from intogen_core.pathways.genesets import read_gene_sets             # noqa: E402

STRATA = ['MSS', 'MSI']
CURATED = ['CRC_WNT', 'CRC_TGFB_BMP', 'CRC_PI3K', 'CRC_RTK_RAS', 'CRC_TP53_DDR', 'CRC_IMMUNE_ESCAPE', 'CRC_SWI_SNF']

# relationships reported in the CRC literature (direction, short reference)
KNOWN = {
    frozenset(['KRAS', 'BRAF']): ('exclusive', 'RAS/RAF alternatives (TCGA 2012 and many others)'),
    frozenset(['KRAS', 'NRAS']): ('exclusive', 'RAS alternatives'),
    frozenset(['NRAS', 'BRAF']): ('exclusive', 'RAS/RAF alternatives'),
    frozenset(['APC', 'BRAF']): ('exclusive', 'serrated (BRAF) vs conventional (APC) pathway'),
    frozenset(['APC', 'CTNNB1']): ('exclusive', 'WNT alternatives'),
    frozenset(['APC', 'RNF43']): ('exclusive', 'WNT alternatives; RNF43 in APC-wild-type, serrated/MSI'),
    frozenset(['BRAF', 'RNF43']): ('co-occurrence', 'serrated / MSI pathway'),
    frozenset(['KRAS', 'PIK3CA']): ('co-occurrence', 'PIK3CA mutations enriched in KRAS-mutant CRC'),
    frozenset(['APC', 'KRAS']): ('co-occurrence', 'adenoma-carcinoma sequence'),
    frozenset(['APC', 'TP53']): ('co-occurrence', 'adenoma-carcinoma sequence'),
    frozenset(['ACVR2A', 'TGFBR2']): ('co-occurrence', 'MSI target genes'),
    frozenset(['PIK3CA', 'PTEN']): ('exclusive', 'PI3K alternatives'),
    frozenset(['BRAF', 'TP53']): ('exclusive', 'BRAF/serrated tumours less often TP53-mutant'),
}


def load(stratum):
    tumours, hits = analysis.mutation_hits(f'results/{stratum}/mutations.tsv.gz')
    layers = {name: cooccurrence.Layer(name, df[['SYMBOL', 'SAMPLE']], tumours) for name, df in hits.items()}
    return tumours, layers


def drivers(stratum):
    sel = pd.read_csv(f'data/{stratum}/sel_cv.tsv', sep='\t')
    return dict(zip(sel.loc[sel.qglobal_cv < 0.1, 'gene_name'], sel.loc[sel.qglobal_cv < 0.1, 'qglobal_cv']))


def candidate_events(stratum, all_drivers, layers, sets):
    own = drivers(stratum)
    events = []
    for gene in sorted(all_drivers):
        if layers['mutation'].altered_tumours([gene]) >= 3:
            events.append(cooccurrence.Event('mutation', 'gene', gene, gene, [gene], own.get(gene, 1.0)))
    for s in CURATED:
        genes = set(sets.get(s, [])) - set(own)
        if layers['mutation'].altered_tumours(genes) >= 3:
            events.append(cooccurrence.Event('mutation', 'pathway', s + '_LONG_TAIL', s + ' (long tail)', genes, 0.5))
    res = pd.read_csv(f'results/{stratum}/pathways.tsv.gz', sep='\t')
    lt = res[(res.SCOPE == 'long_tail') & (res.LAYER != 'combined') & (res.Q_VALUE < 0.1) & ~res.SET.isin(CURATED)]
    for set_id, grp in lt.groupby('SET'):
        best = grp.sort_values(['RATIO', 'Q_VALUE'], ascending=[False, True]).iloc[0]
        genes = set(sets[set_id]) - set(own)
        events.append(cooccurrence.Event(best.LAYER, 'pathway', set_id, set_id, genes, float(best.Q_VALUE)))
    return events


def mine(stratum, all_drivers, sets):
    tumours, layers = load(stratum)
    events = candidate_events(stratum, all_drivers, layers, sets)
    pairs = cooccurrence.test_pairs(events, layers, min_tumours=3)
    pairs['Q_VALUE'], pairs['TESTABLE'] = cooccurrence.tarone_bh(pairs.P_VALUE.values, pairs.P_MIN.values, 0.1)
    pairs['Q_VALUE_EXCLUSIVITY'], pairs['TESTABLE_EXCLUSIVITY'] = cooccurrence.tarone_bh(
        pairs.P_VALUE_EXCLUSIVITY.values, pairs.P_MIN_EXCLUSIVITY.values, 0.1)
    pairs['A'] = pairs.EVENT_1.str.split(':', n=1).str[1]
    pairs['B'] = pairs.EVENT_2.str.split(':', n=1).str[1]
    known = [KNOWN.get(frozenset([a, b])) for a, b in zip(pairs.A, pairs.B)]
    pairs['KNOWN'] = [k[0] if k else '' for k in known]
    pairs['KNOWN_NOTE'] = [k[1] if k else '' for k in known]
    pairs['STRATUM'] = stratum
    pairs['N_TUMOURS'] = len(tumours)
    pairs.to_csv(f'results/{stratum}/network_pairs.tsv.gz', sep='\t', index=False)
    modules = cooccurrence.modules(pairs, 0.1)
    return pairs, modules, len(events), len(tumours)


def site_associations(all_drivers):
    """Driver frequency by location within MSS, and MSS vs MSI"""
    groups = pd.read_csv('tumour_groups.tsv', sep='\t').set_index('PATIENT')
    rows = []
    for stratum in ['MSS', 'MSI']:
        tumours, layers = load(stratum)
        loc = groups.reindex(tumours)['LOCATION'].fillna('NA').values
        for gene in sorted(all_drivers):
            hits, _ = layers['mutation'].event([gene])
            if hits.sum() < 5:
                continue
            if stratum == 'MSS':
                ok = loc != 'NA'
                table = pd.crosstab(pd.Categorical(loc[ok], ['proximal', 'distal', 'rectum']), hits[ok])
                if table.shape[1] == 2:
                    chi2, p, _, _ = sps.chi2_contingency(table.values)
                    freq = (table[True] / table.sum(axis=1)).to_dict()
                    rows.append({'TEST': 'location (MSS)', 'GENE': gene, 'P': p, 'MUTATED': int(hits.sum()),
                                 'FREQ_PROXIMAL': freq['proximal'], 'FREQ_DISTAL': freq['distal'],
                                 'FREQ_RECTUM': freq['rectum']})
    # MSI vs MSS (driver frequency)
    t_mss, l_mss = load('MSS')
    t_msi, l_msi = load('MSI')
    for gene in sorted(all_drivers):
        a = int(l_mss['mutation'].event([gene])[0].sum())
        b = int(l_msi['mutation'].event([gene])[0].sum())
        if a + b < 5:
            continue
        odds, p = sps.fisher_exact([[b, len(t_msi) - b], [a, len(t_mss) - a]])
        rows.append({'TEST': 'MSI vs MSS', 'GENE': gene, 'P': p, 'MUTATED': a + b, 'FREQ_MSS': a / len(t_mss),
                     'FREQ_MSI': b / len(t_msi), 'ODDS': odds})
    df = pd.DataFrame(rows)
    for test in df.TEST.unique():
        idx = df.TEST == test
        df.loc[idx, 'Q'] = fdr_bh(df.loc[idx, 'P'].values)
    df.to_csv('results/site_msi_associations.tsv', sep='\t', index=False)
    return df


if __name__ == '__main__':
    sets_table = read_gene_sets('gene_sets_crc.tsv.gz')
    sets = {s: list(g) for s, g in sets_table.groupby('SET')['SYMBOL']}
    all_drivers = set()
    for s in ['ALL', 'MSS', 'MSI']:
        all_drivers |= set(drivers(s))
    print('drivers in any stratum:', len(all_drivers))
    summary = {}
    for stratum in STRATA:
        pairs, modules, n_events, n_tumours = mine(stratum, all_drivers, sets)
        ex = pairs[pairs.Q_VALUE_EXCLUSIVITY < 0.1].sort_values('P_VALUE_EXCLUSIVITY')
        co = pairs[pairs.Q_VALUE < 0.1].sort_values('P_VALUE')
        print(f'== {stratum}: {n_tumours} tumours, {n_events} events, {len(pairs)} pairs, '
              f'{len(co)} co-occurring, {len(ex)} exclusive')
        for label, sub, col in [('EXCLUSIVE', ex, 'Q_VALUE_EXCLUSIVITY'), ('CO-OCCURRING', co, 'Q_VALUE')]:
            for _, r in sub.iterrows():
                tag = f'known {r.KNOWN}' if r.KNOWN else 'CANDIDATE'
                if r.KNOWN and ((r.KNOWN == 'exclusive') != (label == 'EXCLUSIVE')):
                    tag = f'OPPOSITE to known {r.KNOWN}'
                print(f'   {label:12s} {r.A} | {r.B}: {r.TUMOURS_BOTH} vs {r.EXPECTED_BOTH:.1f} expected '
                      f'(n1={r.TUMOURS_1}, n2={r.TUMOURS_2}) q={r[col]:.2g} [{tag}]')
        summary[stratum] = {'tumours': n_tumours, 'events': n_events, 'pairs': int(len(pairs)),
                            'modules': {m: sorted(e for e, x in modules.items() if x == m) for m in set(modules.values())}}
        print('   modules:', summary[stratum]['modules'])
    assoc = site_associations(all_drivers)
    for test in assoc.TEST.unique():
        sub = assoc[(assoc.TEST == test) & (assoc.Q < 0.1)].sort_values('P')
        print(f'== {test}: {len(sub)} genes with q < 0.1')
        print(sub.drop(columns=['TEST']).round(3).to_string(index=False))
    json.dump(summary, open('results/network_summary.json', 'w'), indent=1)

    # coherence: are the long tails of the curated pathways exclusive with the drivers of the same pathway?
    checks = {'MSS': {'CRC_RTK_RAS_LONG_TAIL': ['KRAS', 'NRAS', 'BRAF'], 'CRC_TGFB_BMP_LONG_TAIL': ['SMAD4', 'SMAD2', 'SMAD3', 'ACVR2A'],
                      'CRC_WNT_LONG_TAIL': ['APC', 'CTNNB1', 'RNF43', 'AMER1', 'TCF7L2', 'SOX9'],
                      'CRC_PI3K_LONG_TAIL': ['PIK3CA', 'PTEN', 'PIK3R1'], 'CRC_TP53_DDR_LONG_TAIL': ['TP53', 'ATM']},
              'MSI': {'CRC_IMMUNE_ESCAPE_LONG_TAIL': ['B2M', 'HLA-B', 'TAP2', 'CD58', 'CASP8'],
                      'CRC_TGFB_BMP_LONG_TAIL': ['ACVR1B', 'BMPR2'], 'CRC_RTK_RAS_LONG_TAIL': ['KRAS', 'BRAF']}}
    coherence = []
    for stratum, d in checks.items():
        p = pd.read_csv(f'results/{stratum}/network_pairs.tsv.gz', sep='\t')
        for tail, genes in d.items():
            sub = p[((p.A == tail) & p.B.isin(genes)) | ((p.B == tail) & p.A.isin(genes))]
            chi = -2 * np.log(np.clip(sub.P_VALUE_EXCLUSIVITY.values, 1e-300, 1)).sum()
            coherence.append({'stratum': stratum, 'long_tail': tail, 'observed': int(sub.TUMOURS_BOTH.sum()),
                              'expected': float(sub.EXPECTED_BOTH.sum()), 'p_exclusivity': float(sps.chi2.sf(chi, 2 * len(sub)))})
    print(pd.DataFrame(coherence).round(3).to_string(index=False))
    pd.DataFrame(coherence).to_csv('results/coherence.tsv', sep='\t', index=False)
