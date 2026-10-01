"""
Inputs of the replication from an IntOGen-style cohort mutation table (columns CHROMOSOME, POSITION, REF, ALT, SAMPLE,
optionally DONOR and read counts), e.g. the table prepared for IntOGen in the 100,000 Genomes Research Environment.

    python prepare_intogen.py <input dir> <table.tsv.gz> [<table2.tsv.gz> ...]

- One sample per donor when a DONOR column is present (the first in alphabetical order, as IntOGen does).
- Read counts, when present, are kept for allele fractions: T_ALT and T_DEPTH, or t_alt_count with t_depth or
  t_ref_count, or a VAF / AF column.
- VCF-style indels (shared first base) are written MAF-style.
- Use the table before IntOGen's parse-variants step: that step removes hypermutated samples (WGS: more than 10,000
  SNVs and above Q3 + 1.5 IQR of the cohort), which drops MSI and POLE tumours and can drop high-burden MSS tumours.
  A warning is printed when the table looks filtered.
"""
import gzip
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from prepare_vcf import maf_alleles  # noqa: E402

REQUIRED = {'CHROMOSOME': ['CHROMOSOME', 'CHROM', 'CHR'], 'POSITION': ['POSITION', 'POS', 'START'],
            'REF': ['REF'], 'ALT': ['ALT'], 'SAMPLE': ['SAMPLE']}


def find(columns, names):
    upper = {c.upper(): c for c in columns}
    return next((upper[n] for n in names if n in upper), None)


CHUNK = 2_000_000


def columns(path):
    """Column names of a table: required ones, DONOR and read counts"""
    header = list(pd.read_csv(path, sep='\t', nrows=0).columns)
    cols = {k: find(header, v) for k, v in REQUIRED.items()}
    missing = [k for k, v in cols.items() if v is None]
    if missing:
        raise SystemExit(f'{path}: columns {missing} not found (have {header[:12]})')
    cols['DONOR'] = find(header, ['DONOR'])
    cols['T_ALT'] = find(header, ['T_ALT', 'T_ALT_COUNT'])
    cols['T_DEPTH'] = find(header, ['T_DEPTH'])
    cols['T_REF'] = find(header, ['T_REF', 'T_REF_COUNT'])
    cols['VAF'] = find(header, ['VAF', 'AF', 'TUMOUR_VAF', 'TUMOR_VAF'])
    if cols['T_ALT'] and (cols['T_DEPTH'] or cols['T_REF']):
        cols['READS'] = f"{cols['T_ALT']} and {cols['T_DEPTH'] or cols['T_REF']}"
    elif cols['VAF']:
        cols['READS'] = cols['VAF']
    else:
        cols['READS'] = None
    return cols


def chunks(path, cols):
    use = [c for k, c in cols.items() if k != 'READS' and c]
    for t in pd.read_csv(path, sep='\t', dtype=str, usecols=use, chunksize=CHUNK):
        out = pd.DataFrame({'SAMPLE': t[cols['SAMPLE']],
                            'CHROM': t[cols['CHROMOSOME']].str.replace('^chr', '', regex=True),
                            'POS': t[cols['POSITION']].astype(float).astype(int),
                            'REF': t[cols['REF']].fillna('-'), 'ALT': t[cols['ALT']].fillna('-')})
        out['DONOR'] = t[cols['DONOR']] if cols['DONOR'] else out.SAMPLE
        if cols['T_ALT'] and (cols['T_DEPTH'] or cols['T_REF']):
            out['T_ALT'] = pd.to_numeric(t[cols['T_ALT']], errors='coerce')
            out['T_DEPTH'] = (pd.to_numeric(t[cols['T_DEPTH']], errors='coerce') if cols['T_DEPTH']
                              else out.T_ALT + pd.to_numeric(t[cols['T_REF']], errors='coerce'))
        elif cols['VAF']:
            out['VAF'] = pd.to_numeric(t[cols['VAF']], errors='coerce')
        yield out


def is_snv(m):
    return (m.REF.str.len() == 1) & (m.ALT.str.len() == 1) & (m.REF != '-') & (m.ALT != '-')


def main(out_dir, *tables):
    specs = {p: columns(p) for p in tables}
    if len({c['READS'] is None for c in specs.values()}) > 1:
        raise SystemExit('some tables have read counts and others do not')
    for p in tables:
        stats_file = p + '.stats.json'
        if os.path.exists(stats_file):
            removed = json.load(open(stats_file)).get('hypermutators', {}).get('hypermutators', [])
            if removed:
                print(f'WARNING: {p} is the output of parse-variants, which removed {len(removed)} hypermutated samples. '
                      'Use the table before parse-variants.')
    # pass 1: donors and SNVs per sample
    donors, snvs = {}, {}
    for p, c in specs.items():
        for m in chunks(p, c):
            for d, s in m[['DONOR', 'SAMPLE']].drop_duplicates().itertuples(index=False):
                donors.setdefault(d, set()).add(s)
            for s, n in m[is_snv(m)].groupby('SAMPLE').size().items():
                snvs[s] = snvs.get(s, 0) + int(n)
    keep = {min(s) for s in donors.values()}
    dropped = set().union(*donors.values()) - keep
    # pass 2: kept samples, MAF-style indels
    os.makedirs(out_dir, exist_ok=True)
    total, converted, header = 0, 0, True
    with gzip.open(os.path.join(out_dir, 'mutations.tsv.gz'), 'wt') as out:
        for p, c in specs.items():
            for m in chunks(p, c):
                m = m[m.SAMPLE.isin(keep)].drop(columns='DONOR')
                vcf_style = ((m.REF.str.len() != m.ALT.str.len()) & (m.REF.str[0] == m.ALT.str[0])
                             & ~m.REF.str.contains('-', regex=False) & ~m.ALT.str.contains('-', regex=False))
                if vcf_style.any():
                    conv = [maf_alleles(q, r, a) for q, r, a in zip(m.POS[vcf_style], m.REF[vcf_style], m.ALT[vcf_style])]
                    m.loc[vcf_style, ['POS', 'REF', 'ALT']] = pd.DataFrame(conv, index=m.index[vcf_style],
                                                                          columns=['POS', 'REF', 'ALT'])
                m = m.drop_duplicates(['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT'])
                m.to_csv(out, sep='\t', index=False, header=header)
                header = False
                total += len(m)
                converted += int(vcf_style.sum())

    per_sample = pd.Series({s: n for s, n in snvs.items() if s in keep}, dtype=float).reindex(sorted(keep)).fillna(0)
    reads = next(iter(specs.values()))['READS']
    print(f'tumours {len(keep)}, mutations {total}, VCF-style indels converted {converted}, '
          f'samples dropped as second samples of a donor {len(dropped)}')
    print(f'read counts: {reads or "none (allele fractions and the purity proxy are unavailable; give PURITY in samples.tsv)"}')
    q = np.percentile(per_sample, [50, 90, 99]).astype(int)
    print(f'SNVs per tumour: median {q[0]}, 90th percentile {q[1]}, 99th {q[2]}, maximum {int(per_sample.max())}; '
          f'tumours with more than 10,000 SNVs: {int((per_sample > 10000).sum())}')
    if per_sample.median() > 2000 and (per_sample > 10000).sum() == 0:
        print('WARNING: whole genomes but no tumour above 10,000 SNVs: the table looks filtered by parse-variants '
              '(hypermutated tumours removed). Use the table before parse-variants.')


if __name__ == '__main__':
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    main(sys.argv[1], *sys.argv[2:])
