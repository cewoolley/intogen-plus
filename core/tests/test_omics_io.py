import numpy as np
import pandas as pd
import pytest

from intogen_core.omics import io


def test_samplesheet_and_sample_assignment(tmp_path):
    sheet = tmp_path / 'samples.tsv'
    sheet.write_text('ID\tSAMPLE\tTYPE\n'
                     'A_rep2\tA\ttumor\n'
                     'A_rep1\tA\t\n'
                     'N1\tA\tNormal\n'
                     'B\t\t\n')
    parsed = io.read_samplesheet(str(sheet))
    assert parsed['N1'] == ('A', io.NORMAL)
    assert parsed['A_rep1'] == ('A', io.TUMOR)
    assert parsed['B'] == ('B', io.TUMOR)

    tumors, normals, discarded = io.assign_samples(['A_rep2', 'A_rep1', 'N1', 'B', 'C'], parsed)
    # one column per sample: the first in alphabetical order
    assert tumors == {'A': 'A_rep1', 'B': 'B', 'C': 'C'}
    assert normals == {'A': 'N1'}
    assert discarded == ['A_rep2']


def test_samplesheet_errors(tmp_path):
    bad = tmp_path / 'bad.tsv'
    bad.write_text('ID\tTYPE\nX\tmetastasis\n')
    with pytest.raises(io.OmicsError):
        io.read_samplesheet(str(bad))
    no_id = tmp_path / 'noid.tsv'
    no_id.write_text('SAMPLE\nX\n')
    with pytest.raises(io.OmicsError):
        io.read_samplesheet(str(no_id))
    assert io.read_samplesheet(None) == {}


def test_gene_mapper():
    annotation = pd.DataFrame({'ENSEMBL_GENE': ['ENSG00000000001', 'ENSG00000000002'],
                               'SYMBOL': ['TP53', 'TSR3'], 'CHROMOSOME': ['17', 'X']})
    mapper = io.GeneMapper(annotation, {'C16orf42': 'TSR3'})
    assert mapper.map('ENSG00000000001.12') == 'TP53'
    assert mapper.map('ENSG00000000001') == 'TP53'
    assert mapper.map('TP53|7157') == 'TP53'
    assert mapper.map('C16orf42') == 'TSR3'
    assert mapper.map('?|100130426') is None
    assert mapper.map('LINC00001') is None       # not a MANE protein coding gene
    assert mapper.map('ENSG99999999999') is None
    assert mapper.is_sex_chromosome('TSR3') and not mapper.is_sex_chromosome('TP53')
    # without annotation symbols are kept
    assert io.GeneMapper().map('ANYGENE') == 'ANYGENE'


def test_read_matrix_drops_text_columns_and_filters_rows(tmp_path):
    path = tmp_path / 'm.csv.gz'
    df = pd.DataFrame({'id': ['a', 'b', 'c'], 'desc': ['x', 'y', 'z'], 'S1': [1, 2, 3], 'S2': ['1', 'NA', '2']})
    df.to_csv(path, index=False)
    m, dropped, total = io.read_matrix(str(path), keep=lambda idx: idx != 'b', chunksize=2)
    assert dropped == ['desc']
    assert total == 3
    assert list(m.index) == ['a', 'c']
    assert list(m.columns) == ['S1', 'S2']
    assert m.loc['c', 'S2'] == 2


def test_check_outputs_refuses_to_overwrite_the_input(tmp_path):
    source = tmp_path / 'input.tsv.gz'
    source.write_text('x')
    link = tmp_path / 'X.expression_matrix.tsv.gz'
    link.symlink_to(source)
    with pytest.raises(io.OmicsError):
        io.check_outputs(str(source), str(link))
    io.check_outputs(str(source), str(tmp_path / 'other.tsv.gz'))


def test_matrix_roundtrip(tmp_path):
    df = pd.DataFrame({'S1': [0.1, np.nan], 'S2': [0.123456789, 1.0]}, index=['G1', 'G2'])
    path = str(tmp_path / 'm.tsv.gz')
    io.write_matrix(df, path)
    back = io.read_gene_matrix(path)
    assert list(back.index) == ['G1', 'G2']
    assert np.isnan(back.loc['G2', 'S1'])
    assert back.loc['G1', 'S2'] == pytest.approx(0.12346)
