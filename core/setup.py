from os import path
from setuptools import setup, find_packages


VERSION = "0.2"
DESCRIPTION = "Intogen Core functionality"

directory = path.dirname(path.abspath(__file__))

# Get requirements from the requirements.txt file
with open(path.join(directory, 'requirements.txt')) as f:
    required = f.read().splitlines()


# Get the long description from the README file
with open(path.join(directory, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='intogen_core',
    version=VERSION,
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/x-rst',
    url="",
    author="Barcelona Biomedical Genomics Lab",
    author_email="bbglab@irbbarcelona.org",
    license="Apache Software License 2.0",
    packages=find_packages(),
    install_requires=required,
    entry_points={
            'console_scripts': [
                'parse-variants = intogen_core.parsers.variants:cli',
                'format-variants = intogen_core.formatters.main:cli',
                'parse-vep = intogen_core.parsers.vep:cli',
                'parse-nonsynonymous = intogen_core.parsers.nonsynonymous:cli',
                'parse-mnvs = intogen_core.parsers.mnvs:cli',
                'parse-profile = intogen_core.parsers.profile:cli',
                'mutations-summary = intogen_core.postprocess.mutations:cli',
                'drivers-discovery = intogen_core.postprocess.drivers.discovery:cli',
                'drivers-summary = intogen_core.postprocess.drivers.summary:cli',
                'drivers-saturation = intogen_core.postprocess.drivers.saturation:cli',
                'parse-methylation = intogen_core.omics.methylation:parse_cli',
                'methylation-analysis = intogen_core.omics.methylation:analysis_cli',
                'parse-expression = intogen_core.omics.expression:parse_cli',
                'expression-analysis = intogen_core.omics.expression:analysis_cli',
                'omics-features = intogen_core.omics.features:cli',
                'omics-summary = intogen_core.omics.features:summary_cli'
            ]
        },
)