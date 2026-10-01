"""
MSS / MSI / POLE groups, coding mutation burden and location of every tumour (README.md, samples.tsv).

    python classify.py --check-status <samples.tsv>          exit 0 when every tumour has an MSI_STATUS
    python classify.py <samples.tsv> <coding.tsv.gz> <sample_vaf.tsv> <groups.tsv> [<annotmuts.tsv.gz>]

Groups come from the cohort's MSI_STATUS when it is given for every tumour. Otherwise they are derived from the
mutations as in TCGA (spec.py): POLE = missense mutation in the POLE exonuclease domain and hypermutated,
MSI = other hypermutated tumours, MSS = the rest.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spec  # noqa: E402

STATUS = {'MSS': 'MSS', 'MSI-L': 'MSS', 'MSI-LOW': 'MSS', 'MSI': 'MSI', 'MSI-H': 'MSI', 'MSI-HIGH': 'MSI', 'POLE': 'POLE'}


def read_samples(path):
    s = pd.read_csv(path, sep='\t', dtype={'SAMPLE': str})
    if s.SAMPLE.duplicated().any():
        raise SystemExit('samples.tsv: one row per tumour sample is expected')
    return s


def has_status(s):
    return 'MSI_STATUS' in s.columns and s.MSI_STATUS.notna().all()


def location(s):
    if 'LOCATION' in s.columns:
        return s.LOCATION.where(s.LOCATION.isin(['proximal', 'distal', 'rectum']), 'NA')
    if 'SITE_ICD10' in s.columns:
        code = s.SITE_ICD10.astype(str).str.upper().str.strip()
        return pd.Series(np.select([code.isin(spec.PROXIMAL), code.isin(spec.DISTAL), code.isin(spec.RECTUM)],
                                   ['proximal', 'distal', 'rectum'], 'NA'), index=s.index)
    return pd.Series('NA', index=s.index)


def main(samples_path, coding_path, vaf_path, out_path, annotation=None):
    s = read_samples(samples_path).set_index('SAMPLE')
    m = pd.read_csv(coding_path, sep='\t', dtype={'SAMPLE': str, 'REF': str, 'ALT': str})
    indel = (m.REF.fillna('').str.replace('-', '', regex=False).str.len()
             != m.ALT.fillna('').str.replace('-', '', regex=False).str.len())
    g = m.assign(INDEL=indel).groupby('SAMPLE')
    s['CODING'] = g.size().reindex(s.index).fillna(0).astype(int)
    s['INDEL_FRACTION'] = (g.INDEL.sum() / g.size()).reindex(s.index)
    vaf = pd.read_csv(vaf_path, sep='\t', dtype={'SAMPLE': str}).set_index('SAMPLE').MEDIAN_VAF
    s['MEDIAN_VAF'] = vaf.reindex(s.index)
    if has_status(s):
        st = s.MSI_STATUS.astype(str).str.upper().str.strip()
        unknown = sorted(set(st) - set(STATUS))
        if unknown:
            raise SystemExit(f'MSI_STATUS values not understood: {unknown}')
        s['GROUP'] = st.map(STATUS)
        s['GROUP_SOURCE'] = 'cohort'
    else:
        a = pd.read_csv(annotation, sep='\t', dtype={'sampleID': str})
        pole = a[(a.gene == 'POLE') & (a.impact == 'Missense')].copy()
        pole['codon'] = pole.aachange.str.extract(r'^[A-Z*](\d+)', expand=False).astype(float)
        lo, hi = spec.POLE_DOMAIN
        ed = set(pole[(pole.codon >= lo) & (pole.codon <= hi)].sampleID)
        hyper = s.CODING > spec.HYPERMUTATED_CODING
        s['GROUP'] = np.where(s.index.isin(ed) & hyper, 'POLE', np.where(hyper, 'MSI', 'MSS'))
        s['GROUP_SOURCE'] = 'mutations'
    s['LOCATION'] = location(s.reset_index()).values
    s.index.name = 'SAMPLE'
    s.to_csv(out_path, sep='\t')
    print('groups:', s.GROUP.value_counts().to_dict(), '| source:', s.GROUP_SOURCE.iloc[0],
          '| median coding mutations:', s.groupby('GROUP').CODING.median().to_dict())


if __name__ == '__main__':
    if sys.argv[1] == '--check-status':
        sys.exit(0 if has_status(read_samples(sys.argv[2])) else 1)
    main(*sys.argv[1:6])
