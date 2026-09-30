"""Gene sets called by each approach on the real cohorts"""
import json, sys
import numpy as np, pandas as pd
from scipy import stats as sps
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'core'))
from intogen_core.omics.stats import fdr_bh
from intogen_core.pathways.genesets import GeneSets, read_gene_sets
ARTEFACT = ('NEURONAL', 'SYNAP', 'POTASSIUM', 'GPCR', 'MUSCLE', 'EXTRACELLULAR_MATRIX', 'L1CAM', 'NEUROTRANSMITTER',
            'ION_CHANNEL', 'SENSORY', 'OLFACTORY', 'CARDIAC', 'CHEMICAL_SYNAPSES', 'NMDA', 'GABA', 'CALCIUM')
out = {}
for c in ['KIRC', 'LUAD', 'BRCA', 'COAD']:
    gm = pd.read_csv(f'data/{c}/genemuts.tsv', sep='\t')
    drivers = set(pd.read_csv(f'results/{c}/vet.tsv', sep='\t').SYMBOL)
    universe = set(gm.gene_name)
    sets = GeneSets(read_gene_sets('data/gene_sets.tsv.gz'), universe, 10, 500)
    # pre-fork + ORA: over-representation of the drivers
    n, N = len(drivers & universe), len(universe)
    ora = {s: (sps.hypergeom.sf(len(sets.genes[s] & drivers) - 1, N, len(sets.genes[s]), n)
               if sets.genes[s] & drivers else 1.0) for s in sets.ids}
    ora_q = dict(zip(ora, fdr_bh(np.array(list(ora.values())))))
    # dNdScv genesetdnds on the long tails (one-sided, all substitutions)
    gsd = pd.read_csv(f'results/{c}/gsd.tsv', sep='\t')
    gsd = gsd[gsd.UNIT.str.startswith('LT|')].copy()
    gsd['SET'] = gsd.UNIT.str[3:]
    gsd['P'] = np.where(gsd.wall > 1, gsd.p_wall / 2, 1 - gsd.p_wall / 2)
    gsd['Q'] = fdr_bh(gsd.P.values)
    res = {'ORA of drivers (pre-fork + standard enrichment)': {s for s, q in ora_q.items() if q < 0.1},
           'dNdScv genesetdnds (long tail)': set(gsd[gsd.Q < 0.1].SET)}
    for tag, label in [('fork', 'fork as shipped (long tail)'), ('final', 'fork corrected (long tail)')]:
        df = pd.read_csv(f'results/{c}/{tag}.pathways.tsv.gz', sep='\t')
        lt = df[(df.SCOPE == 'long_tail') & (df.LAYER == 'combined')]
        res[label] = set(lt[lt.Q_VALUE < 0.1].SET)
    out[c] = {}
    for label, called in res.items():
        art = sorted(s for s in called if any(a in s for a in ARTEFACT))
        out[c][label] = {'n': len(called), 'neuronal_muscle_ecm_channel': len(art), 'sets': sorted(called)}
        print(f'{c} {label:52s} {len(called):4d} sets, {len(art):3d} neuronal/muscle/ECM/channel-type')
    top = gsd.sort_values('P').head(6)
    print('   genesetdnds top:', '; '.join(f'{s[:40]} w={w:.2f} q={q:.1e}' for s, w, q in zip(top.SET, top.wall, top.Q)))
    top_ora = sorted(ora_q.items(), key=lambda x: x[1])[:5]
    print('   ORA top:', '; '.join(f'{s[:40]} q={q:.1e}' for s, q in top_ora))
json.dump(out, open('results/sets_compare.json', 'w'), indent=1)
