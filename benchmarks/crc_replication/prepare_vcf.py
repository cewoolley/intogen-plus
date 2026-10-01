"""
Inputs of the replication from somatic VCFs (e.g. the 100,000 Genomes Project cancer programme, Strelka2 calls).

    python prepare_vcf.py <manifest.tsv> <input dir> [--regions regions.bed] [--check N]

manifest.tsv: SAMPLE and VCF columns, optionally TUMOUR_COLUMN (name of the tumour sample in the VCF; default TUMOR,
else the last sample). Only PASS records are kept. Read counts: Strelka2 tier-1 counts ({REF}U/{ALT}U for SNVs,
TAR/TIR for indels), else AD. Indels are written MAF-style (shared first base removed, '-' for empty alleles).
With --regions (dnds.R regions) only coding and splice-site mutations are written, and the median allele fraction of
every tumour over all its PASS mutations goes to sample_vaf.tsv (run.sh uses it as the purity proxy).
--check N prints the first N parsed records of the first VCF and exits: compare them with the VCF before a full run.
"""
import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from coding_filter import load_regions, overlaps  # noqa: E402


def counts(fmt, values, ref, alt):
    f = dict(zip(fmt.split(':'), values.split(':')))
    if 'TIR' in f and 'TAR' in f:                                     # Strelka2 indel
        a, r = int(f['TIR'].split(',')[0]), int(f['TAR'].split(',')[0])
        return a, a + r
    if len(ref) == 1 and len(alt) == 1 and f'{alt}U' in f and f'{ref}U' in f:      # Strelka2 SNV
        a, r = int(f[f'{alt}U'].split(',')[0]), int(f[f'{ref}U'].split(',')[0])
        return a, a + r
    if 'AD' in f:
        ad = [int(x) for x in f['AD'].split(',') if x not in ('.', '')]
        if len(ad) >= 2:
            return ad[1], sum(ad)
    return None, None


def maf_alleles(pos, ref, alt):
    if len(ref) != len(alt) or len(ref) > 1:
        while ref and alt and ref[0] == alt[0]:
            ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref or '-', alt or '-'


def parse(sample, path, tumour_column=None, limit=None):
    opener = gzip.open if path.endswith('.gz') else open
    rows, column = [], None
    with opener(path, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                continue
            fields = line.rstrip('\n').split('\t')
            if line.startswith('#'):
                names = fields[9:]
                name = tumour_column if tumour_column else ('TUMOR' if 'TUMOR' in names else names[-1])
                column = 9 + names.index(name)
                continue
            if fields[6] != 'PASS':
                continue
            alt = fields[4].split(',')[0]
            ref = fields[3]
            t_alt, depth = counts(fields[8], fields[column], ref, alt)
            pos, r, a = maf_alleles(int(fields[1]), ref, alt)
            rows.append((sample, fields[0], pos, r, a, t_alt, depth))
            if limit and len(rows) >= limit:
                break
    return pd.DataFrame(rows, columns=['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'T_ALT', 'T_DEPTH']), name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('manifest')
    ap.add_argument('out')
    ap.add_argument('--regions')
    ap.add_argument('--check', type=int, default=0)
    a = ap.parse_args()
    man = pd.read_csv(a.manifest, sep='\t', dtype=str)
    col = man.TUMOUR_COLUMN if 'TUMOUR_COLUMN' in man.columns else pd.Series([None] * len(man))
    if a.check:
        df, name = parse(man.SAMPLE.iloc[0], man.VCF.iloc[0], col.iloc[0], limit=a.check)
        print(f'tumour column: {name}')
        print(df.to_string(index=False))
        return
    os.makedirs(a.out, exist_ok=True)
    regions = load_regions(a.regions) if a.regions else None
    vafs = []
    with gzip.open(os.path.join(a.out, 'mutations.tsv.gz'), 'wt') as out:
        out.write('SAMPLE\tCHROM\tPOS\tREF\tALT\tT_ALT\tT_DEPTH\n')
        for i, (s, v, c) in enumerate(zip(man.SAMPLE, man.VCF, col)):
            df, _ = parse(s, v, c)
            vaf = (df.T_ALT / df.T_DEPTH.where(df.T_DEPTH > 0)).dropna()
            vafs.append((s, float(vaf.median()) if len(vaf) else np.nan))
            if regions is not None:
                chrom = df.CHROM.str.replace('^chr', '', regex=True).values
                ref_len = df.REF.str.replace('-', '', regex=False).str.len().clip(lower=1).values
                df = df[overlaps(regions, chrom, df.POS.values, df.POS.values + ref_len - 1)]
            df.to_csv(out, sep='\t', header=False, index=False)
            if (i + 1) % 100 == 0:
                print(f'{i + 1} of {len(man)} VCFs', flush=True)
    pd.DataFrame(vafs, columns=['SAMPLE', 'MEDIAN_VAF']).to_csv(os.path.join(a.out, 'sample_vaf.tsv'), sep='\t', index=False)


if __name__ == '__main__':
    main()
