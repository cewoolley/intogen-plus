"""
Gene sets (pathways) used by the pathway analyses.

Two formats are accepted:

- TSV with the columns SET, SOURCE, NAME and SYMBOL (one row per gene and set),
  as built in the datasets from Reactome and the MSigDB hallmarks.
- GMT (``.gmt``): one set per line, the name, a description and the genes,
  separated by tabs.
"""

import gzip
import re

import pandas as pd

from intogen_core.omics.io import OmicsError


COLUMNS = ['SET', 'SOURCE', 'NAME', 'SYMBOL']


def _open(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def read_gene_sets(path):
    """
    Returns:
        DataFrame with the columns SET, SOURCE, NAME and SYMBOL
    """
    if re.search(r'\.gmt(\.gz)?$', str(path)):
        rows = []
        with _open(path) as fd:
            for line in fd:
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 3 or fields[0].strip() == '':
                    continue
                name = fields[0].strip()
                rows += [(name, 'custom', name, g.strip()) for g in fields[2:] if g.strip()]
        return pd.DataFrame(rows, columns=COLUMNS)

    df = pd.read_csv(path, sep='\t', dtype=str)
    missing = set(COLUMNS) - set(df.columns)
    if missing:
        raise OmicsError(f'Missing columns in the gene sets file {path}: {sorted(missing)}')
    return df[COLUMNS].dropna()


class GeneSets:
    """
    Gene sets restricted to the analysed genes and filtered by size.

    Args:
        table: gene sets (read_gene_sets)
        universe: genes analysed in the cohort
        min_size, max_size: limits of the number of analysed genes of a set
        symbol_map: outdated to current HUGO symbols
    """

    def __init__(self, table, universe, min_size=10, max_size=500, symbol_map=None):
        universe = set(universe)
        table = table.copy()
        if symbol_map:
            table['SYMBOL'] = [s if s in universe else symbol_map.get(s, s) for s in table['SYMBOL']]
        table = table[table['SYMBOL'].isin(universe)].drop_duplicates(['SET', 'SYMBOL'])
        sizes = table.groupby('SET').size()
        keep = sizes[(sizes >= min_size) & (sizes <= max_size)].index
        table = table[table['SET'].isin(keep)]

        self.genes = table.groupby('SET')['SYMBOL'].apply(frozenset).to_dict()
        info = table.drop_duplicates('SET').set_index('SET')
        self.source = info['SOURCE'].to_dict()
        self.name = info['NAME'].to_dict()
        self.ids = sorted(self.genes)
        self.discarded = int(len(sizes) - len(keep))

    def __len__(self):
        return len(self.ids)
