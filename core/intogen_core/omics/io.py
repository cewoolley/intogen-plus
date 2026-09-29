"""
Input/output helpers shared by the omics steps.

Omics inputs are feature x sample matrices (TSV or CSV, optionally gzipped)
whose first column holds the feature identifiers (probe IDs, gene symbols or
Ensembl gene IDs) and whose remaining columns are samples.

Sample columns can be described with an optional sample sheet (TSV) with
the columns:

- ``ID``: column name in the omics matrix (required)
- ``SAMPLE``: sample identifier used in the mutation data (defaults to ``ID``)
- ``TYPE``: ``tumor`` or ``normal`` (defaults to ``tumor``)
"""

import json
import os
import re

import numpy as np
import pandas as pd

from intogen_core.exceptions import IntogenError


TUMOR = 'tumor'
NORMAL = 'normal'

_TYPES = {'tumor': TUMOR, 'tumour': TUMOR, 't': TUMOR, 'normal': NORMAL, 'n': NORMAL}

SEX_CHROMOSOMES = {'X', 'Y', 'MT', 'M'}

ENSEMBL_GENE = re.compile(r'^(ENSG\d+)(\.\d+)?(_PAR_Y)?$')


class OmicsError(IntogenError):
    pass


def _separator(path):
    return ',' if re.search(r'\.csv(\.gz|\.bz2|\.xz)?$', str(path)) else '\t'


def _to_numeric(df):
    """Convert columns to floats, dropping non-numeric columns (e.g. gene descriptions)"""
    dropped = []
    for c in df.columns:
        if pd.api.types.is_numeric_dtype(df[c]):
            continue
        converted = pd.to_numeric(df[c], errors='coerce')
        if converted.notna().any() and converted.notna().sum() >= 0.5 * df[c].notna().sum():
            df[c] = converted
        else:
            dropped.append(c)
    return df.drop(columns=dropped), dropped


def read_matrix(path, keep=None, chunksize=200000, dtype=np.float32):
    """
    Read a feature x sample matrix.

    Args:
        path: file path
        keep: optional callable receiving the index (feature IDs) of a chunk
            and returning a boolean mask of the rows to keep. Used to avoid
            loading all the probes of large methylation arrays in memory.
        chunksize: rows per chunk
        dtype: numeric type of the values

    Returns:
        (DataFrame, list of dropped non-numeric columns, number of rows read)
    """
    chunks, dropped, total = [], None, 0
    reader = pd.read_csv(path, sep=_separator(path), index_col=0, chunksize=chunksize, low_memory=False)
    for chunk in reader:
        total += len(chunk)
        chunk.index = chunk.index.astype(str).str.strip()
        if keep is not None:
            chunk = chunk[np.asarray(keep(chunk.index), dtype=bool)]
        if dropped is None:
            chunk, dropped = _to_numeric(chunk)
        else:
            chunk = chunk.drop(columns=dropped)
            chunk, _ = _to_numeric(chunk)
        chunks.append(chunk.astype(dtype))

    if dropped is None:
        raise OmicsError(f'Empty matrix: {path}')

    df = pd.concat(chunks) if len(chunks) > 1 else chunks[0]
    df.columns = [str(c).strip() for c in df.columns]
    df.index.name = 'ID'
    if df.shape[1] == 0:
        raise OmicsError(f'No numeric sample columns found in {path}')
    return df, dropped, total


def check_outputs(input_file, *outputs):
    """Refuse to write an output over the input (e.g. through a staged symlink)"""
    for output in outputs:
        if os.path.exists(output) and os.path.realpath(output) == os.path.realpath(input_file):
            raise OmicsError(f'The output {output} would overwrite the input. Rename the input file')


def read_samplesheet(path):
    """
    Read the (optional) sample sheet.

    Returns:
        dict mapping matrix column ID to (SAMPLE, TYPE)
    """
    if path is None:
        return {}
    df = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    df.columns = [c.strip().upper() for c in df.columns]
    if 'ID' not in df.columns:
        raise OmicsError(f'The sample sheet {path} must contain an ID column')

    sheet = {}
    for _, row in df.iterrows():
        id_ = row['ID'].strip()
        if id_ == '':
            continue
        sample = row.get('SAMPLE', '').strip() or id_
        type_ = row.get('TYPE', '').strip().lower() or TUMOR
        if type_ not in _TYPES:
            raise OmicsError(f'Invalid sample type "{type_}" for {id_}. Use tumor or normal')
        if id_ in sheet:
            raise OmicsError(f'Duplicated ID in the sample sheet: {id_}')
        sheet[id_] = (sample, _TYPES[type_])
    return sheet


def assign_samples(columns, sheet):
    """
    Assign each matrix column to a sample and a sample type.

    Only one column per sample and type is kept. As in the processing of the
    mutations, when several columns map to the same sample the first one
    in alphabetical order is selected.

    Returns:
        tumors: dict sample -> column
        normals: dict sample -> column
        discarded: list of columns discarded because of duplicated samples
    """
    tumors, normals, discarded = {}, {}, []
    for col in sorted(columns):
        sample, type_ = sheet.get(col, (col, TUMOR))
        target = tumors if type_ == TUMOR else normals
        if sample in target:
            discarded.append(col)
            continue
        target[sample] = col
    return tumors, normals, discarded


def default_dataset(*parts):
    """Path to a file within the IntOGen datasets, or None if unavailable"""
    folder = os.environ.get('INTOGEN_DATASETS')
    if folder is None:
        return None
    path = os.path.join(folder, *parts)
    return path if os.path.exists(path) else None


def load_gene_annotation(path):
    """
    Load the MANE protein coding genes used by IntOGen (cds_biomart.tsv).

    Returns:
        DataFrame with columns ENSEMBL_GENE, SYMBOL, CHROMOSOME (one row per gene)
    """
    if path is None:
        return None
    df = pd.read_csv(path, sep='\t', header=None, usecols=[0, 1, 3], dtype=str,
                     names=['ENSEMBL_GENE', 'SYMBOL', 'CHROMOSOME'])
    return df.drop_duplicates().dropna(subset=['SYMBOL']).reset_index(drop=True)


def load_symbol_map(path):
    """Mapping of outdated HUGO symbols into current ones"""
    if path is None:
        return {}
    with open(path) as fd:
        return json.load(fd)


class GeneMapper:
    """
    Harmonise gene identifiers into the HUGO symbols used by IntOGen.

    Accepted identifiers are HUGO symbols, Ensembl gene IDs (with or without
    version) and ``SYMBOL|ENTREZ`` pairs (e.g. TCGA legacy matrices).
    When the gene annotation is available, only MANE protein coding genes
    (the genes analysed by IntOGen) are kept.
    """

    def __init__(self, annotation=None, symbol_map=None):
        self.annotation = annotation
        self.symbol_map = symbol_map or {}
        if annotation is not None:
            self.ensembl = dict(zip(annotation['ENSEMBL_GENE'], annotation['SYMBOL']))
            self.symbols = set(annotation['SYMBOL'])
            self.chromosome = dict(zip(annotation['SYMBOL'], annotation['CHROMOSOME']))
        else:
            self.ensembl, self.symbols, self.chromosome = {}, None, {}

    def map(self, identifier):
        identifier = str(identifier).strip()
        if '|' in identifier:
            identifier = identifier.split('|')[0]
        if identifier in ('', '?', '-', 'nan', 'NA'):
            return None
        m = ENSEMBL_GENE.match(identifier)
        if m:
            return self.ensembl.get(m.group(1))
        if self.symbols is not None and identifier not in self.symbols:
            identifier = self.symbol_map.get(identifier, identifier)
        if self.symbols is not None and identifier not in self.symbols:
            return None
        return identifier

    def is_sex_chromosome(self, symbol):
        return self.chromosome.get(symbol, '') in SEX_CHROMOSOMES


def write_matrix(df, path):
    """Write a gene x sample matrix"""
    df.to_csv(path, sep='\t', na_rep='NA', float_format='%.5g', index_label='SYMBOL')


def read_gene_matrix(path):
    """Read a gene x sample matrix written by :func:`write_matrix`"""
    df = pd.read_csv(path, sep='\t', index_col=0, dtype={'SYMBOL': str})
    df.index = df.index.astype(str)
    return df.astype(float)


def write_table(df, path):
    df.to_csv(path, sep='\t', index=False, na_rep='NA', float_format='%.6g')


def _json_default(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return None if np.isnan(o) else float(o)
    if isinstance(o, (np.ndarray, set)):
        return list(o)
    return str(o)


def write_stats(path, stats):
    with open(path, 'w') as fd:
        json.dump(stats, fd, indent=4, default=_json_default)
