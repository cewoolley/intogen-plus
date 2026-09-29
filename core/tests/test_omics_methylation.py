import json

import numpy as np
import pandas as pd
import pytest

import synthetic_omics as synthetic
from intogen_core.omics import io, methylation


def read(path, **kwargs):
    return pd.read_csv(path, sep='\t', **kwargs)


def test_parse_probe_level(cohort_files):
    cohort = cohort_files['cohort']
    tumors = io.read_gene_matrix(cohort_files['meth_tumor'])
    normals = io.read_gene_matrix(cohort_files['meth_normal'])
    assert list(tumors.columns) == cohort.tumors
    assert list(normals.columns) == cohort.normals
    # genes in sex chromosomes are removed (X inactivation, sex of the patients)
    assert 'XGENE' not in tumors.index and 'YGENE' not in tumors.index
    assert 'SIL1' in tumors.index
    stats = json.load(open(cohort_files['meth_tumor'] + '.stats.json'))
    assert stats['input_level'] == 'probe'
    assert stats['values'] == 'beta'
    assert stats['probes_promoter'] == stats['probes_total'] - 1   # one probe outside promoters
    assert stats['normals'] == len(cohort.normals)


def test_parse_m_values_and_epicv2_ids(tmp_path, cohort_files):
    cohort = synthetic.Cohort(seed=cohort_files['seed'])
    files = cohort.write_methylation(str(tmp_path), 'C', values='m')
    # EPICv2 style identifiers in the matrix
    m = pd.read_csv(files['matrix'], sep='\t', index_col=0)
    m.index = [i + '_TC21' for i in m.index]
    m.to_csv(files['matrix'], sep='\t')
    out, out_n = str(tmp_path / 't.tsv.gz'), str(tmp_path / 'n.tsv.gz')
    methylation.parse(files['matrix'], out, out_n, probes=files['probes'], samples=cohort_files['samples'],
                      gene_annotation=cohort_files['annotation'])
    assert json.load(open(out + '.stats.json'))['values'] == 'm'
    from_m = io.read_gene_matrix(out)
    from_beta = io.read_gene_matrix(cohort_files['meth_tumor'])
    assert np.allclose(from_m.values, from_beta.loc[from_m.index, from_m.columns].values, atol=1e-3, equal_nan=True)


def test_parse_gene_level(tmp_path, cohort_files):
    tumors = io.read_gene_matrix(cohort_files['meth_tumor'])
    matrix = tumors * 100   # percentages
    matrix.index = [f'ENSG{synthetic.genes().index(g):011d}.1' for g in matrix.index]
    path = str(tmp_path / 'gene_level.tsv')
    matrix.to_csv(path, sep='\t')
    out, out_n = str(tmp_path / 't.tsv.gz'), str(tmp_path / 'n.tsv.gz')
    methylation.parse(path, out, out_n, gene_annotation=cohort_files['annotation'])
    stats = json.load(open(out + '.stats.json'))
    assert stats['input_level'] == 'gene' and stats['values'] == 'percent'
    parsed = io.read_gene_matrix(out)
    assert np.allclose(parsed.values, tumors.loc[parsed.index].values, atol=1e-3, equal_nan=True)
    assert io.read_gene_matrix(out_n).shape[1] == 0


def test_probe_level_requires_annotation(tmp_path, cohort_files):
    with pytest.raises(io.OmicsError):
        methylation.parse(cohort_files['matrix'], str(tmp_path / 'a.tsv.gz'), str(tmp_path / 'b.tsv.gz'))


def test_epigenetic_silencing(cohort_files):
    result = read(cohort_files['meth']).set_index('SYMBOL')
    sil1 = result.loc['SIL1']
    assert sil1['STATUS'] == 'tested' and sil1['BASELINE_SOURCE'] == 'normal'
    assert sil1['Q_VALUE'] < 1e-10
    assert sil1['FUNCTIONAL'] is True or sil1['FUNCTIONAL'] == 'True'
    assert sil1['LOG2FC_SILENCING'] < -2
    # hypermethylated but still expressed: not functional, not a candidate
    hyp1 = result.loc['HYP1']
    assert hyp1['P_RECURRENCE'] < 1e-10
    assert hyp1['FUNCTIONAL'] in (False, 'False') and hyp1['P_VALUE'] == 1
    # methylated in normal samples
    assert result.loc['MET1', 'STATUS'] == 'baseline_methylated'
    # background genes
    background = result[result.index.str.match(r'^G\d+$')]
    assert (background['Q_VALUE'] < 0.05).sum() <= 1

    events = read(cohort_files['meth_events'])
    sil_events = events[events['SYMBOL'] == 'SIL1']
    assert len(sil_events) >= len(cohort_files['cohort'].silenced)
    assert sil_events['LOG2_EXPRESSION'].notna().all()
    assert (sil_events['DELTA_BETA'] >= 0.2).all()


def test_epigenetic_silencing_without_normals_or_expression(tmp_path, cohort_files):
    out, events = str(tmp_path / 'm.tsv.gz'), str(tmp_path / 'e.tsv.gz')
    methylation.analyse(cohort_files['meth_tumor'], out, events)
    result = read(out).set_index('SYMBOL')
    assert (result.loc[result['STATUS'] == 'tested', 'BASELINE_SOURCE'] == 'tumor_quantile').all()
    # without expression recurrent hypermethylation cannot be told apart from silencing
    assert result.loc['SIL1', 'Q_VALUE'] < 1e-10 and result.loc['HYP1', 'Q_VALUE'] < 1e-10
    assert result['FUNCTIONAL'].isna().all()
    assert result.loc['MET1', 'STATUS'] == 'baseline_methylated'
    assert read(events)['LOG2_EXPRESSION'].isna().all()
