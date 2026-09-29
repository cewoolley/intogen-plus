import json

import numpy as np
import pandas as pd
import pytest

from intogen_core.omics import expression, io


def test_detect_units():
    assert expression.detect_units(np.array([[0, 10, 3], [5, 100, 2]], dtype=float)) == 'counts'
    assert expression.detect_units(np.array([[0.5, 10.2], [3.1, 1200.7]])) == 'tpm'
    assert expression.detect_units(np.array([[0.5, 10.2], [3.1, 12.7]])) == 'log2'
    assert expression.detect_units(np.array([[-0.5, 10.2], [3.1, 12.7]])) == 'log2'
    rsem = np.random.default_rng(0).gamma(0.5, 2000, size=(20000, 3)) + 0.37
    assert expression.detect_units(rsem) == 'counts'


def test_normalise():
    counts = np.array([[10, 0], [990, 500], [0, 500.0]])
    log = expression.normalise(counts, 'counts')
    assert log[0, 0] == pytest.approx(np.log2(1e4 + 1))
    assert expression.normalise(np.array([[3.0]]), 'tpm')[0, 0] == pytest.approx(2)
    assert expression.normalise(np.array([[3.0]]), 'log2')[0, 0] == 3
    with pytest.raises(io.OmicsError):
        expression.normalise(np.array([[-1.0, 2.0]]), 'tpm')


def test_parse(cohort_files):
    cohort = cohort_files['cohort']
    m = io.read_gene_matrix(cohort_files['expr_matrix'])
    stats = json.load(open(cohort_files['expr_matrix'] + '.stats.json'))
    assert stats['units'] == 'counts'
    assert list(m.columns) == cohort.tumors
    assert set(m.index) == set(cohort.genes)   # Ensembl IDs mapped to symbols


def test_parse_ignores_count_summary_rows(tmp_path, cohort_files):
    df = pd.DataFrame({'T1': [100, 900, 5000], 'T2': [300, 700, 0]}, index=['SIL1', 'OUT1', '__no_feature'])
    df.loc['N_unmapped'] = [10, 10]
    path = str(tmp_path / 'htseq.tsv')
    df.to_csv(path, sep='\t')
    out = str(tmp_path / 'o.tsv.gz')
    expression.parse(path, out, units='counts', gene_annotation=cohort_files['annotation'])
    m = io.read_gene_matrix(out)
    # library sizes only include genes: 100 / 1000 reads -> 1e5 CPM
    assert m.loc['SIL1', 'T1'] == pytest.approx(np.log2(1e5 + 1), rel=1e-4)
    assert json.load(open(out + '.stats.json'))['summary_rows_removed'] == ['__no_feature', 'N_unmapped']


def test_parse_duplicates_and_normals(tmp_path, cohort_files):
    df = pd.DataFrame({'T1': [1.0, 5.0, 2.0], 'T2': [1.0, 6.0, 3.0], 'N000': [1, 1, 1]},
                      index=['SIL1', 'SIL1', 'OUT1'])
    path = str(tmp_path / 'e.tsv')
    df.to_csv(path, sep='\t')
    out = str(tmp_path / 'o.tsv.gz')
    expression.parse(path, out, samples=cohort_files['samples'], units='log2', gene_annotation=cohort_files['annotation'])
    m = io.read_gene_matrix(out)
    assert list(m.columns) == ['T1', 'T2']
    assert m.loc['SIL1', 'T2'] == 6.0    # most expressed entry
    stats = json.load(open(out + '.stats.json'))
    assert stats['normals_ignored'] == 1 and stats['genes_duplicated'] == 1


def test_outliers_and_expression_status(cohort_files):
    result = pd.read_csv(cohort_files['expr'], sep='\t').set_index('SYMBOL')
    assert result.loc['OUT1', 'DIRECTION'] == 'over' and result.loc['OUT1', 'Q_VALUE'] < 1e-10
    assert result.loc['OUT1', 'OVER_SAMPLES'] == len(cohort_files['cohort'].over)
    assert result.loc['LOW1', 'DIRECTION'] == 'under' and result.loc['LOW1', 'Q_VALUE'] < 1e-10
    assert result.loc['NOEXP', 'EXPRESSED'] in (False, 'False')
    assert result.loc['G1', 'EXPRESSED'] in (True, 'True')
    assert result.loc['YGENE', 'STATUS'] == 'sex_chromosome' and np.isnan(result.loc['YGENE', 'P_VALUE'])
    background = result[result.index.str.match(r'^G\d+$')]
    assert (background['Q_VALUE'] < 0.05).sum() <= 1


def test_non_expressed_threshold_is_inclusive(tmp_path):
    # 8 of 10 samples below 1 TPM -> not expressed; 7 of 10 -> expressed
    values = np.log2(np.array([[0.5] * 8 + [50] * 2, [0.5] * 7 + [50] * 3]) + 1)
    df = pd.DataFrame(values, index=['A', 'B'], columns=[f'S{i}' for i in range(10)])
    path = str(tmp_path / 'm.tsv.gz')
    io.write_matrix(df, path)
    out = str(tmp_path / 'r.tsv.gz')
    expression.analyse(path, out, str(tmp_path / 'e.tsv.gz'))
    result = pd.read_csv(out, sep='\t').set_index('SYMBOL')
    assert result.loc['A', 'EXPRESSED'] in (False, 'False')
    assert result.loc['B', 'EXPRESSED'] in (True, 'True')
