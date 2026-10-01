"""Hypermutation, POLE and MSI groups of the TCGA COADREAD tumours (from their mutations), with their location"""
import numpy as np
import pandas as pd

m = pd.read_csv('coadread_mutations.tsv.gz', sep='\t', dtype=str)
clin = pd.read_csv('clinical.tsv', sep='\t', dtype=str).drop_duplicates('PATIENT').set_index('PATIENT')
g = m.groupby('PATIENT')
t = pd.DataFrame({'CODING': g.size(), 'INDELS': g['type'].apply(lambda s: s.isin(['DEL', 'INS']).sum())})
t['INDEL_FRACTION'] = t.INDELS / t.CODING
# POLE exonuclease domain (codons 268-471) missense mutations
pole = m[(m.gene == 'POLE') & (m['class'] == 'Missense_Mutation')].copy()
pole['codon'] = pole.protein.str.extract(r'p\.[A-Z](\d+)', expand=False).astype(float)
t['POLE_ED'] = t.index.isin(pole[(pole.codon >= 268) & (pole.codon <= 471)].PATIENT)
t = t.join(clin[['MMR_IHC_LOSS', 'SITE', 'PROJECT']])
# the coding mutation burden is bimodal with a gap between 300 and 600 mutations
t['GROUP'] = np.where(t.POLE_ED, 'POLE', np.where(t.CODING > 500, 'MSI', 'MSS'))
site = t.SITE.fillna('NA')
t['LOCATION'] = np.select([site.isin(['Cecum', 'Ascending Colon', 'Hepatic Flexure', 'Transverse Colon']),
                           site.isin(['Splenic Flexure', 'Descending Colon', 'Sigmoid Colon']),
                           site.isin(['Rectosigmoid Junction', 'Rectum'])], ['proximal', 'distal', 'rectum'], 'NA')
t.index.name = 'PATIENT'
t.to_csv('tumour_groups.tsv', sep='\t')
print(t.GROUP.value_counts().to_dict(), t.groupby('GROUP')[['CODING', 'INDEL_FRACTION']].median().to_dict())
