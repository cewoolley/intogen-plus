import os
import sys
import types

import pytest

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# bgdata is only needed by pipeline steps that are not tested here
# (it is imported but not used by intogen_core.postprocess.drivers.role)
try:
    import bgdata  # noqa: F401
except ImportError:
    sys.modules['bgdata'] = types.ModuleType('bgdata')

import synthetic  # noqa: E402

SEED = 2


@pytest.fixture(scope='session')
def cohort_files(tmp_path_factory):
    """Synthetic cohort with methylation, expression and mutations processed by the omics steps"""
    from intogen_core.omics import expression, methylation

    folder = tmp_path_factory.mktemp('cohort')
    path = lambda name: str(folder / name)
    cohort = synthetic.Cohort(seed=SEED)
    synthetic.write_gene_annotation(path('cds_biomart.tsv'), cohort.genes)
    files = cohort.write_methylation(str(folder), 'COHORT')
    files['annotation'] = path('cds_biomart.tsv')
    files['samples'] = synthetic.write_samplesheet(path('samples.tsv'), cohort)
    files['expression_input'] = cohort.write_expression(str(folder), 'COHORT', ensembl=True)
    files['vep'] = cohort.write_vep(path('COHORT.tsv.gz'))

    files['meth_tumor'] = path('COHORT.promoter_methylation.tsv.gz')
    files['meth_normal'] = path('COHORT.promoter_methylation.normal.tsv.gz')
    methylation.parse(files['matrix'], files['meth_tumor'], files['meth_normal'], probes=files['probes'],
                      samples=files['samples'], gene_annotation=files['annotation'])

    files['expr_matrix'] = path('COHORT.expression_matrix.tsv.gz')
    expression.parse(files['expression_input'], files['expr_matrix'], samples=files['samples'],
                     gene_annotation=files['annotation'])

    files['expr'] = path('COHORT.expression.tsv.gz')
    files['expr_events'] = path('COHORT.expression_events.tsv.gz')
    expression.analyse(files['expr_matrix'], files['expr'], files['expr_events'], gene_annotation=files['annotation'])

    files['meth'] = path('COHORT.methylation.tsv.gz')
    files['meth_events'] = path('COHORT.methylation_events.tsv.gz')
    methylation.analyse(files['meth_tumor'], files['meth'], files['meth_events'],
                        normal_file=files['meth_normal'], expression_file=files['expr_matrix'])
    files['cohort'] = cohort
    files['seed'] = SEED
    files['folder'] = str(folder)
    return files
