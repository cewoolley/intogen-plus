"""
Keep the mutations that touch a coding or splice-site region (dnds.R regions) and record the median allele fraction
of every tumour over all its mutations (purity proxy when the cohort has no purity estimates).

    python coding_filter.py <mutations.tsv.gz> <regions.bed> <coding.tsv.gz> <sample_vaf.tsv>
"""
import sys

import numpy as np
import pandas as pd


def load_regions(path):
    bed = pd.read_csv(path, sep='\t', header=None, names=['chr', 'start', 'end'], dtype={'chr': str})
    regions = {}
    for c, d in bed.groupby('chr'):
        d = d.sort_values('start')
        starts, ends = d.start.values + 1, d.end.values          # 1-based, inclusive
        # merge overlapping intervals so that searchsorted finds the only candidate
        merged_s, merged_e = [starts[0]], [ends[0]]
        for s, e in zip(starts[1:], ends[1:]):
            if s <= merged_e[-1] + 1:
                merged_e[-1] = max(merged_e[-1], e)
            else:
                merged_s.append(s)
                merged_e.append(e)
        regions[c] = (np.array(merged_s), np.array(merged_e))
    return regions


def overlaps(regions, chrom, start, end):
    out = np.zeros(len(chrom), dtype=bool)
    for c in np.unique(chrom):
        if c not in regions:
            continue
        s, e = regions[c]
        idx = chrom == c
        i = np.searchsorted(s, end[idx], side='right') - 1      # last region starting at or before the mutation end
        ok = i >= 0
        hit = np.zeros(idx.sum(), dtype=bool)
        hit[ok] = e[i[ok]] >= start[idx][ok]
        out[idx] = hit
    return out


def main(mutations, regions_bed, coding_out, vaf_out):
    regions = load_regions(regions_bed)
    kept, vafs = [], []
    total = 0
    for chunk in pd.read_csv(mutations, sep='\t', dtype={'SAMPLE': str, 'CHROM': str, 'REF': str, 'ALT': str}, chunksize=2_000_000):
        total += len(chunk)
        chrom = chunk.CHROM.str.replace('^chr', '', regex=True).values
        start = chunk.POS.values.astype(np.int64)
        ref_len = chunk.REF.fillna('').str.replace('-', '', regex=False).str.len().clip(lower=1).values
        end = start + ref_len - 1
        if {'T_ALT', 'T_DEPTH'} <= set(chunk.columns):
            vafs.append(pd.DataFrame({'SAMPLE': chunk.SAMPLE.values, 'VAF': chunk.T_ALT / chunk.T_DEPTH.where(chunk.T_DEPTH > 0)}))
        kept.append(chunk[overlaps(regions, chrom, start, end)])
    coding = pd.concat(kept)
    coding.to_csv(coding_out, sep='\t', index=False)
    if vafs:
        v = pd.concat(vafs).dropna()
        v.groupby('SAMPLE').VAF.median().rename('MEDIAN_VAF').to_csv(vaf_out, sep='\t')
    else:
        pd.DataFrame(columns=['SAMPLE', 'MEDIAN_VAF']).to_csv(vaf_out, sep='\t', index=False)
    print(f'mutations: {total}, in coding or splice-site regions: {len(coding)}, tumours: {coding.SAMPLE.nunique()}')


if __name__ == '__main__':
    main(*sys.argv[1:5])
