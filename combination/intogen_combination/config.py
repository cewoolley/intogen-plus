
import os

from configobj import ConfigObj


CODE_DIR = os.path.dirname(os.path.abspath(__file__))
configs_file = os.path.join(CODE_DIR, 'intogen_qc.cfg')
CONF = ConfigObj(configs_file)


def _flag(method, key):
    return str(CONF[method].get(key, 'False')).strip().lower() in ('true', 'yes', '1')


ALL_METHODS = list(CONF.keys())

# Methods that only take part in the combination when their results are provided
# (omics-based methods: methylation and expression)
OPTIONAL_METHODS = [m for m in ALL_METHODS if _flag(m, 'OPTIONAL')]

# Mutation-based driver identification methods, always part of the combination
METHODS = [m for m in ALL_METHODS if m not in OPTIONAL_METHODS]

# Methods whose results are restricted to the candidate genes of the cohort
# (genes with at least two mutated samples, i.e. with an OncodriveFML q-value)
RESTRICTED_METHODS = [m for m in ALL_METHODS if _flag(m, 'RESTRICT_TO_CANDIDATES')]


def active_methods(files):
    """Methods taking part in the combination given the method results provided"""
    return METHODS + [m for m in OPTIONAL_METHODS if m in files]


REGIONS = os.path.join(os.environ['INTOGEN_DATASETS'], 'regions', 'cds.regions.gz')
