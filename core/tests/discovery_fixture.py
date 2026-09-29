"""
Minimal datasets and inputs to run the driver discovery (drivers-discovery).

Genes:

- DRV1: CGC gene
- DRV2: non-CGC gene with 2 significant methods (expressed according to TCGA)
- DRV3: non-CGC gene with 2 significant methods, not expressed according to TCGA
- DRV4: non-CGC gene, significant for a mutation-based method and methylation
- DRV5: non-CGC gene, significant only for omics-based methods
- NODRV: tier 4 gene
"""

import json
import os

import pandas as pd


GENES = ['DRV1', 'DRV2', 'DRV3', 'DRV4', 'DRV5', 'NODRV']


def make_datasets(folder):
    for sub in ['others', 'postprocess', 'cgc', 'regions']:
        os.makedirs(os.path.join(folder, sub), exist_ok=True)
    p = lambda *parts: os.path.join(folder, *parts)
    with open(p('others', 'non_expressed_genes_tcga.tsv'), 'w') as fd:
        fd.write('PRAD\tDRV3,OTHER\nPANCANCER\tOTHER\n')
    pd.DataFrame({'gene': GENES, 'canonical': True, 'oe_syn': 1.0, 'oe_lof': 0.5, 'oe_mis': 0.9}).to_csv(
        p('postprocess', 'constraint.txt.gz'), sep='\t', index=False, compression='gzip')
    pd.DataFrame({'Symbol': ['OR1A1']}).to_csv(p('others', 'olfactory_receptors.tsv'), sep='\t', index=False)
    with open(p('postprocess', 'artifacts.json'), 'w') as fd:
        json.dump({'suspects': [], 'known': []}, fd)
    for name in ['black_listed.txt', 'white_listed.txt']:
        open(p('postprocess', name), 'w').close()
    pd.DataFrame({'gene_normalized': ['DRV2', 'DRV3', 'DRV4', 'DRV5'], 'pmid': [1, 2, 3, 4]}).to_csv(
        p('postprocess', 'cancermine_sentences.tsv'), sep='\t', index=False)
    pd.DataFrame({'Gene Symbol': ['DRV1'], 'cancer_type': ['PRAD'], 'Tier': [1], 'Role in Cancer': ['TSG']}).to_csv(
        p('cgc', 'cancer_gene_census_parsed.tsv'), sep='\t', index=False)
    rows = [[f'ENSG{i}', g, f'ENSP{i}', '1', 1, 2, 1, 2, 3, 1, f'ENST{i}', 1, 2] for i, g in enumerate(GENES)]
    pd.DataFrame(rows).to_csv(p('regions', 'cds_biomart.tsv'), sep='\t', header=False, index=False)


def make_inputs(folder):
    os.makedirs(folder, exist_ok=True)
    p = lambda name: os.path.join(folder, name)
    files = {}

    bidders = {'DRV1': 'dndscv,oncodrivefml,cbase', 'DRV2': 'dndscv,oncodrivefml', 'DRV3': 'dndscv,cbase',
               'DRV4': 'oncodrivefml,methylation', 'DRV5': 'methylation,expression', 'NODRV': 'dndscv'}
    pd.DataFrame({
        'SYMBOL': GENES, 'TIER': [1, 1, 1, 1, 1, 4],
        'All_Bidders': [bidders[g] for g in GENES], 'Significant_Bidders': [bidders[g] for g in GENES],
        'QVALUE_stouffer_w': [1e-8, 1e-5, 1e-5, 1e-4, 1e-4, 0.3], 'QVALUE_CGC_stouffer_w': [1e-8] + [None] * 5,
        'RANKING': [1, 2, 3, 4, 5, 6], 'MUTS': [10, 6, 6, 5, 5, 3], 'SAMPLES': [10, 6, 6, 5, 5, 3],
        'wmis_cv': [5, 3, 1, 1, 1, 1], 'wnon_cv': [1, 1, 8, 1, 1, 1], 'wspl_cv': 1.0, 'n_mis': 5, 'n_non': 2,
        'ROLE': ['Act', 'Act', 'LoF', 'ambiguous', 'ambiguous', 'ambiguous'],
    }).to_csv(p('combination.05.out.gz'), sep='\t', index=False, compression='gzip')
    files['combination'] = p('combination.05.out.gz')

    rows = []
    for g, n in zip(GENES, [10, 6, 6, 5, 5, 3]):
        for i in range(n):
            rows.append(['1', 1000 + i, 'C', 'T', f'S{i}', f'ID{g}{i}', g, 'ACG', 'A[C>T]G'])
    pd.DataFrame(rows, columns=['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE', 'ID', 'GENE', 'CONTEXT',
                                'MUTATION_TYPE']).to_csv(p('mutations.in.tsv.gz'), sep='\t', index=False,
                                                         compression='gzip')
    files['mutations'] = p('mutations.in.tsv.gz')

    pd.DataFrame({'Sample': [f'S{i}' for i in range(10)], 'Mutation_type': 'A[C>T]G',
                  'Signature.1': 0.8, 'Signature.9': 0.1, 'Signature.10': 0.1}).to_csv(
        p('sig_likelihood'), sep='\t', index=False)
    files['sig_likelihood'] = p('sig_likelihood')

    pd.DataFrame({'HUGO_SYMBOL': ['DRV1'], 'Q_VALUE': [0.01], 'OBSERVED_REGION': [5], 'MEAN_SIMULATED': [1],
                  'REGION': ['ENST0:PF00001:10:50']}).to_csv(p('smregions.tsv.gz'), sep='\t', index=False)
    files['smregions'] = p('smregions.tsv.gz')
    pd.DataFrame({'SYMBOL': ['DRV2'], 'P': [0.01], 'COORDINATES': ['30,20']}).to_csv(
        p('clusters.tsv'), sep='\t', index=False)
    files['clustl_clusters'] = p('clusters.tsv')
    pd.DataFrame({'HUGO Symbol': ['DRV1'], 'CRAVAT Res': [12], 'q-value': [0.01]}).to_csv(
        p('hotmaps.clusters.gz'), sep='\t', index=False, compression='gzip')
    files['hotmaps'] = p('hotmaps.clusters.gz')
    pd.DataFrame({'gene_name': GENES, 'n_mis': 5, 'wmis_cv': 3.0, 'n_non': 2, 'wnon_cv': 1.5, 'n_spl': 0,
                  'wspl_cv': 1.0}).to_csv(p('dndscv.tsv.gz'), sep='\t', index=False, compression='gzip')
    files['dndscv'] = p('dndscv.tsv.gz')
    return files


def make_omics(path):
    """Omics features: DRV2 not expressed in the cohort, DRV3 expressed (unlike TCGA)"""
    pd.DataFrame({
        'SYMBOL': ['DRV1', 'DRV2', 'DRV3', 'DRV4'],
        'EXPRESSED': [True, False, True, True],
        'QVALUE_METHYLATION': [0.5, None, None, 0.001],
        'METHYLATION_FUNCTIONAL': [None, None, None, True],
        'QVALUE_EXPRESSION': [0.01, 0.9, 0.9, 0.9],
        'EXPRESSION_OUTLIER_DIRECTION': ['over', 'under', 'under', 'under'],
        'OMICS_ROLE_SUPPORT': ['Act', None, None, 'LoF'],
    }).to_csv(path, sep='\t', index=False, compression='gzip')
    return path
