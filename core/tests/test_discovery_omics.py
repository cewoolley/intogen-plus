import pandas as pd
import pytest

import discovery_fixture
from intogen_core.postprocess.drivers import omics


@pytest.fixture()
def discovery(tmp_path, monkeypatch):
    datasets = str(tmp_path / 'datasets')
    discovery_fixture.make_datasets(datasets)
    monkeypatch.setenv('INTOGEN_DATASETS', datasets)
    files = discovery_fixture.make_inputs(str(tmp_path / 'inputs'))
    from intogen_core.postprocess.drivers import discovery as module

    def run(omics_file=None):
        out_drivers, out_vet = str(tmp_path / 'COHORT.drivers.tsv'), str(tmp_path / 'COHORT.vet.tsv')
        module.run(files['combination'], files['mutations'], files['sig_likelihood'], 'COHORT', 'PRAD',
                   files['smregions'], files['clustl_clusters'], files['hotmaps'], files['dndscv'],
                   out_drivers, out_vet, omics=omics_file)
        return pd.read_csv(out_drivers, sep='\t'), pd.read_csv(out_vet, sep='\t').set_index('SYMBOL')

    return run, module, tmp_path


def test_apply_cohort_expression():
    df = pd.DataFrame({'GENE': ['A', 'B', 'C'], 'Warning_Expression': [False, True, False]})
    features = pd.DataFrame({'SYMBOL': ['A', 'B'], 'EXPRESSED': ['False', True]})
    out = omics.apply_cohort_expression(df, features)
    assert list(out['Warning_Expression']) == [True, False, False]
    assert list(out['Warning_Expression_Source']) == ['cohort', 'cohort', 'TCGA']


def test_has_mutational_bidder():
    assert omics.has_mutational_bidder('dndscv,methylation')
    assert not omics.has_mutational_bidder('methylation,expression')
    assert not omics.has_mutational_bidder(float('nan'))


def test_discovery_without_omics(discovery):
    run, module, _ = discovery
    drivers, vet = run()
    assert list(drivers.columns) == module.DRIVERS_COLUMNS
    assert list(vet.reset_index().columns) == module.VET_COLUMNS
    assert vet.loc['DRV2', 'FILTER'] == 'PASS'
    assert vet.loc['DRV3', 'FILTER'] == 'Warning expression'   # TCGA
    assert vet.loc['DRV4', 'FILTER'] == 'PASS'
    # non CGC genes need a mutation-based significant method
    assert vet.loc['DRV5', 'FILTER'] == 'No driver'
    assert vet.loc['NODRV', 'FILTER'] == 'No driver'


def test_summary_with_and_without_omics(discovery, monkeypatch):
    run, module, tmp_path = discovery
    from intogen_core.postprocess.drivers import summary

    outputs = []
    for cohort, omics_file in [('WITH', discovery_fixture.make_omics(str(tmp_path / 'WITH.omics.tsv.gz'))),
                               ('WITHOUT', None)]:
        drivers, vet = run(omics_file)
        d, v = str(tmp_path / f'{cohort}.drivers.tsv'), str(tmp_path / f'{cohort}.vet.tsv')
        drivers.assign(COHORT=cohort).to_csv(d, sep='\t', index=False)
        vet.reset_index().to_csv(v, sep='\t', index=False)
        outputs.append((d, v))
    pd.DataFrame({'TRANSCRIPT': 'ENST0', 'SYMBOL': discovery_fixture.GENES * 2, 'SAMPLES': 3,
                  'COHORT': ['WITH'] * 6 + ['WITHOUT'] * 6}).to_csv(tmp_path / 'mutations.tsv', sep='\t', index=False)
    pd.DataFrame({'COHORT': ['WITH', 'WITHOUT'], 'CANCER_TYPE': 'PRAD', 'PLATFORM': 'WXS', 'MUTATIONS': 50,
                  'SAMPLES': 10}).to_csv(tmp_path / 'cohorts.tsv', sep='\t', index=False)
    monkeypatch.chdir(tmp_path)
    summary.run('mutations.tsv', 'cohorts.tsv', [o[0] for o in outputs], [o[1] for o in outputs])

    drivers = pd.read_csv('drivers.tsv', sep='\t', dtype=str)
    assert list(drivers.columns[-len(omics.OMICS_DRIVERS_COLUMNS):]) == omics.OMICS_DRIVERS_COLUMNS
    with_omics = drivers[drivers['COHORT'] == 'WITH'].set_index('SYMBOL')
    assert with_omics.loc['DRV1', 'EXPRESSED'] == 'True'
    assert drivers.loc[drivers['COHORT'] == 'WITHOUT', 'EXPRESSED'].isna().all()
    unfiltered = pd.read_csv('unfiltered_drivers.tsv', sep='\t')
    assert 'WARNING_EXPRESSION_SOURCE' in unfiltered.columns


def test_discovery_with_omics(discovery):
    run, module, tmp_path = discovery
    omics_file = discovery_fixture.make_omics(str(tmp_path / 'COHORT.omics.tsv.gz'))
    drivers, vet = run(omics_file)
    assert list(drivers.columns) == module.DRIVERS_COLUMNS + omics.OMICS_DRIVERS_COLUMNS
    assert list(vet.reset_index().columns) == module.VET_COLUMNS + omics.OMICS_VET_COLUMNS
    # the expression of the cohort replaces TCGA
    assert vet.loc['DRV2', 'FILTER'] == 'Warning expression'
    assert vet.loc['DRV3', 'FILTER'] == 'PASS'
    assert vet.loc['DRV2', 'WARNING_EXPRESSION_SOURCE'] == 'cohort'
    assert vet.loc['NODRV', 'WARNING_EXPRESSION_SOURCE'] == 'TCGA'
    drivers = drivers.set_index('SYMBOL')
    assert drivers.loc['DRV4', 'QVALUE_METHYLATION'] == pytest.approx(0.001)
    assert drivers.loc['DRV4', 'OMICS_ROLE_SUPPORT'] == 'LoF'
    assert drivers.loc['DRV1', 'OMICS_ROLE_SUPPORT'] == 'Act'
