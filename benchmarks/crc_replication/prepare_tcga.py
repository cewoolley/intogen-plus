"""
Inputs of the replication from the TCGA COADREAD benchmark (dry run of the pipeline on the discovery cohort).

    python prepare_tcga.py <benchmarks/tcga_coadread run dir> <input dir>
"""
import os
import sys

import pandas as pd


def main(tcga, out):
    os.makedirs(out, exist_ok=True)
    m = pd.read_csv(os.path.join(tcga, 'coadread_mutations.tsv.gz'), sep='\t', dtype=str)
    reads = pd.read_csv(os.path.join(tcga, 'mutations_reads.tsv.gz'), sep='\t', dtype={'chr': str, 'SAMPLE': str})
    reads = reads.drop_duplicates(['SAMPLE', 'chr', 'pos'])[['SAMPLE', 'chr', 'pos', 't_ref', 't_alt']]
    m['pos'] = m.pos.astype(float).astype(int)          # written as 9e+05 by older exports
    reads['pos'] = reads.pos.astype(float).astype(int)
    m = m.merge(reads, on=['SAMPLE', 'chr', 'pos'], how='left')
    mut = pd.DataFrame({'SAMPLE': m.PATIENT, 'CHROM': m.chr, 'POS': m.pos, 'REF': m.ref, 'ALT': m.mut,
                        'T_ALT': m.t_alt, 'T_DEPTH': m.t_alt + m.t_ref}).drop_duplicates(['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT'])
    mut.to_csv(os.path.join(out, 'mutations.tsv.gz'), sep='\t', index=False)

    groups = pd.read_csv(os.path.join(tcga, 'tumour_groups.tsv'), sep='\t').set_index('PATIENT')
    clin = pd.read_csv(os.path.join(tcga, 'clinical.tsv'), sep='\t').drop_duplicates('PATIENT').set_index('PATIENT')
    s = groups[['GROUP', 'LOCATION']].join(clin[['AGE', 'SEX', 'CDR_STAGE', 'OS', 'OS_TIME', 'PFI', 'PFI_TIME']])
    s = s.rename(columns={'GROUP': 'MSI_STATUS', 'CDR_STAGE': 'STAGE'})
    s['LOCATION'] = s.LOCATION.fillna('NA')
    s.index.name = 'SAMPLE'
    s.to_csv(os.path.join(out, 'samples.tsv'), sep='\t')
    print(f'mutations {len(mut)} ({mut.T_DEPTH.notna().mean():.1%} with read counts), tumours {len(s)}')


if __name__ == '__main__':
    main(*sys.argv[1:3])
