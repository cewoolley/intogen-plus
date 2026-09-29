import numpy as np
import pandas as pd

from intogen_core.omics import features


def test_features_all_layers(tmp_path, cohort_files):
    out = str(tmp_path / 'COHORT.omics.tsv.gz')
    features.run(cohort_files['vep'], out,
                 expression_matrix=cohort_files['expr_matrix'], expression_results=cohort_files['expr'],
                 methylation_matrix=cohort_files['meth_tumor'], methylation_results=cohort_files['meth'],
                 methylation_events=cohort_files['meth_events'])
    df = pd.read_csv(out, sep='\t').set_index('SYMBOL')
    assert list(df.reset_index().columns) == features.OMICS_COLUMNS
    cohort = cohort_files['cohort']

    # SIL1: mutated in some tumours, silenced in others (alternative hits)
    sil1 = df.loc['SIL1']
    assert sil1['OMICS_ROLE_SUPPORT'] == 'LoF'
    assert sil1['SAMPLES_MUTATION_METHYLATION'] == len(cohort.tumors)
    events = pd.read_csv(cohort_files['meth_events'], sep='\t')
    hyper = set(events.loc[events['SYMBOL'] == 'SIL1', 'SAMPLE'])
    mutated = {cohort.tumors[i] for i in cohort.sil_mutated}
    assert set(cohort.tumors[i] for i in cohort.silenced) <= hyper
    assert sil1['SAMPLES_MUTATED_OR_HYPERMETHYLATED'] == len(hyper | mutated)
    assert sil1['SAMPLES_MUTATED_AND_HYPERMETHYLATED'] == len(hyper & mutated)

    # OUT1: over-expressed, mutated tumours are among the over-expressing ones
    out1 = df.loc['OUT1']
    assert out1['OMICS_ROLE_SUPPORT'] == 'Act'
    assert out1['EXPRESSION_MUTANT_SAMPLES'] == len(cohort.out_mutated)
    # (the wild-type tumours include other over-expressing tumours)
    assert out1['EXPRESSION_LOG2FC_MUTANT'] > 3 and out1['PVALUE_EXPRESSION_MUTANT'] < 0.01

    # TRUNC1: truncating mutations with lower expression
    trunc1 = df.loc['TRUNC1']
    assert trunc1['EXPRESSION_TRUNCATING_SAMPLES'] == len(cohort.trunc_mutated)
    assert trunc1['EXPRESSION_LOG2FC_TRUNCATING'] < -1 and trunc1['PVALUE_EXPRESSION_TRUNCATING'] < 1e-4
    # the loss of expression is recurrent enough to be detected by the outlier test
    assert trunc1['EXPRESSION_OUTLIER_DIRECTION'] == 'under' and trunc1['OMICS_ROLE_SUPPORT'] == 'LoF'

    assert df.loc['LOW1', 'OMICS_ROLE_SUPPORT'] == 'LoF'
    assert df.loc['NOEXP', 'EXPRESSED'] in (False, 'False')


def test_role_support_rules():
    role = features.role_support
    assert role({'QVALUE_METHYLATION': 0.01, 'METHYLATION_FUNCTIONAL': True}) == 'LoF'
    assert role({'QVALUE_METHYLATION': 0.01, 'METHYLATION_FUNCTIONAL': 'False'}) is None
    assert role({'QVALUE_METHYLATION': 0.01, 'METHYLATION_FUNCTIONAL': np.nan}) == 'LoF'   # no expression
    assert role({'QVALUE_EXPRESSION': 0.01, 'EXPRESSION_OUTLIER_DIRECTION': 'over'}) == 'Act'
    assert role({'QVALUE_EXPRESSION_MUTANT': 0.01, 'EXPRESSION_LOG2FC_MUTANT': 1.5}) == 'Act'
    # lower expression of mutants is not LoF support (nonsense-mediated decay of passengers)
    assert role({'QVALUE_EXPRESSION_MUTANT': 0.01, 'EXPRESSION_LOG2FC_MUTANT': -1.5}) is None
    assert role({'QVALUE_METHYLATION': 0.01, 'METHYLATION_FUNCTIONAL': True,
                 'QVALUE_EXPRESSION': 0.01, 'EXPRESSION_OUTLIER_DIRECTION': 'over'}) == 'conflicting'
    assert role({'QVALUE_METHYLATION': 0.5, 'QVALUE_EXPRESSION': np.nan}) is None


def test_features_single_layer(tmp_path, cohort_files):
    out = str(tmp_path / 'X.omics.tsv.gz')
    features.run(cohort_files['vep'], out,
                 methylation_matrix=cohort_files['meth_tumor'], methylation_results=cohort_files['meth'],
                 methylation_events=cohort_files['meth_events'])
    df = pd.read_csv(out, sep='\t').set_index('SYMBOL')
    assert df['EXPRESSED'].isna().all() and df['EXPRESSION_LOG2FC_MUTANT'].isna().all()
    assert df.loc['SIL1', 'QVALUE_METHYLATION'] < 1e-10


def test_features_no_matching_samples(tmp_path, cohort_files):
    vep = pd.read_csv(cohort_files['vep'], sep='\t')
    vep['#Uploaded_variation'] = vep['#Uploaded_variation'].str.replace('__T', '__OTHER')
    path = str(tmp_path / 'vep.tsv.gz')
    vep.to_csv(path, sep='\t', index=False)
    out = str(tmp_path / 'Y.omics.tsv.gz')
    features.run(path, out, expression_matrix=cohort_files['expr_matrix'], expression_results=cohort_files['expr'])
    import json
    stats = json.load(open(out + '.stats.json'))
    assert 'warning_expression_samples' in stats


def test_summary(tmp_path, cohort_files):
    first = str(tmp_path / 'A.omics.tsv.gz')
    features.run(cohort_files['vep'], first, expression_results=cohort_files['expr'])
    output = str(tmp_path / 'omics.tsv')
    features.summary([first], output)
    df = pd.read_csv(output, sep='\t')
    assert list(df.columns[:2]) == ['SYMBOL', 'COHORT']
    assert set(df['COHORT']) == {'A'}
    assert {'OUT1', 'LOW1'} <= set(df['SYMBOL'])
    features.summary([], str(tmp_path / 'empty.tsv'))
