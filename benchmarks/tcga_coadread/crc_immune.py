"""
Immune escape in hypermutated TCGA COADREAD (a headline of Cornish et al. 2024): fraction of tumours with
protein-altering or truncating mutations in antigen-presentation / interferon-response genes, by group, and
selection on the gene set (fork's long-tail test) in MSI tumours.
"""
import json
import numpy as np
import pandas as pd

IMMUNE = ['B2M', 'HLA-A', 'HLA-B', 'HLA-C', 'TAP1', 'TAP2', 'TAPBP', 'NLRC5', 'JAK1', 'JAK2', 'IRF1', 'STAT1',
          'IFNGR1', 'IFNGR2', 'CIITA', 'CD58', 'CASP8']
CORE = ['B2M', 'HLA-A', 'HLA-B', 'HLA-C', 'JAK1', 'JAK2', 'CASP8']
AFFECTING = {'Missense_Mutation', 'Nonsense_Mutation', 'Splice_Site', 'Frame_Shift_Del', 'Frame_Shift_Ins',
             'In_Frame_Del', 'In_Frame_Ins', 'Nonstop_Mutation', 'Translation_Start_Site'}
TRUNCATING = {'Nonsense_Mutation', 'Splice_Site', 'Frame_Shift_Del', 'Frame_Shift_Ins'}

m = pd.read_csv('coadread_mutations.tsv.gz', sep='\t', dtype=str)
groups = pd.read_csv('tumour_groups.tsv', sep='\t').set_index('PATIENT')
m = m[m['class'].isin(AFFECTING)]
out = {}
for group in ['MSS', 'MSI', 'POLE']:
    patients = groups.index[groups.GROUP == group]
    sub = m[m.PATIENT.isin(patients)]
    row = {'tumours': int(len(patients))}
    for label, genes in [('any_immune', IMMUNE), ('core_immune', CORE)]:
        hit = sub[sub.gene.isin(genes)]
        row[f'{label}_altered'] = float(hit.PATIENT.nunique() / len(patients))
        row[f'{label}_truncated'] = float(hit[hit['class'].isin(TRUNCATING)].PATIENT.nunique() / len(patients))
    per_gene = sub[sub.gene.isin(IMMUNE)].groupby('gene').PATIENT.nunique().sort_values(ascending=False)
    row['per_gene'] = {g: int(n) for g, n in per_gene.items()}
    trunc_gene = sub[sub.gene.isin(IMMUNE) & sub['class'].isin(TRUNCATING)].groupby('gene').PATIENT.nunique()
    row['per_gene_truncating'] = {g: int(n) for g, n in trunc_gene.sort_values(ascending=False).items()}
    out[group] = row
    print(group, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()})

# selection on the immune-escape set in MSI and MSS (fork, long tail and all genes)
for stratum in ['MSI', 'MSS']:
    df = pd.read_csv(f'results/{stratum}/pathways.tsv.gz', sep='\t')
    sel = df[df.SET == 'CRC_IMMUNE_ESCAPE'][['SCOPE', 'LAYER', 'N_GENES', 'OBSERVED', 'EXPECTED', 'RATIO', 'P_VALUE', 'Q_VALUE', 'TOP_GENES']]
    print(f'== selection of CRC_IMMUNE_ESCAPE in {stratum}\n', sel.to_string(index=False))
    out[f'selection_{stratum}'] = sel.to_dict('records')
json.dump(out, open('results/immune_escape.json', 'w'), indent=1, default=float)
