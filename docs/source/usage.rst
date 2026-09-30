
Installation
------------

The IntOGen pipeline requires `Nextflow <https://www.nextflow.io/>`_
and `Singularity <https://sylabs.io/docs/>`_ in order to run.

Beside them, a number of different datasets need to be downloaded,
and other pieces of software installed (as Singularity containers).

For information on how to download and build all these requirements,
check the `README <https://github.com/bbglab/intogen-plus/blob/v2024/build/>`_
file in the build folder.


Usage
^^^^^

Once all the prerequisites are available, running the pipeline
only requires to execute the :file:`intogen.nf` file with the appropiate
parameters. E.g.:

.. code-block:: bash

      nextflow run intogen.nf --input <input>

There are a number of parameters and options that can be added:


-resume  Nextflow feature to allow for resumable executions.

--output <path>   Path of the input. See below for more details.

--container <path>   Path where to store the output. Default: ``intogen_analysis``.

--datasets <path>   Path to the folder containing the singularity images. Default: ``containers``.

--annotations <file>    Path to the default annotations file. Default: ``config/annotations.txt``. See the input section for more details.

--seed <int>   Seed to be used for reproducibility. This applies to 4 methods: smRegions, OncodriveCLUSTL, OncodriveFML, dNdScv.

--debug <bool>    Ask methods for a more verbose output if set to True.

--methylation <files>   Optional DNA methylation matrices of the cohorts. See :doc:`omics`.

--expression <files>   Optional RNA-seq matrices of the cohorts. See :doc:`omics`.

--omics_samples <file>   Optional sample sheet of the methylation and RNA-seq matrices.

--methylation_values <type>   Type of methylation values: ``auto`` (default), ``beta``, ``m`` or ``percent``.

--methylation_probes <file>   Promoter probes of the methylation arrays. Default: :file:`<datasets>/methylation/promoter_probes.tsv.gz`.

--expression_units <units>   Units of the RNA-seq matrices: ``auto`` (default), ``counts``, ``tpm``, ``fpkm`` or ``log2``.

--integrate_omics <bool>   Add the methylation and expression evidence to the combination of driver identification methods. Default: ``false``.

--pathways <bool>   Selection of gene sets and co-occurrence of dysregulation events. Default: ``false``. See :doc:`pathways`.

--gene_sets <file>   Gene sets of the pathway analysis (TSV or GMT). Default: :file:`<datasets>/pathways/gene_sets.tsv.gz`.


Input & output
^^^^^^^^^^^^^^

Input
*****

Although the pipeline does most of its computations at the cohort level,
the pipeline is prepared to work with multiple cohorts at the same time.

Each cohort must contain, at least, the chromosome, position, ref, alt
and sample. Files are expected to be TSV files with a header line.

.. important:: All mutations should be mapped to the positive strand.
   The strand value is ignored.

In addition, each cohort must be associated with:

- cohort ID (``DATASET``): a unique identifier for each cohort.
- a cancer type (``CANCER``): although any acronym can be used here, we
  recommend to restrict to the acronyms that can be found
  in :file:`extra/data/dictionary_long_name.json`.
- a sequencing platform (``PLATFORM``): ``WXS`` for whole exome sequencing
  and ``WGS`` for whole genome sequencing
- a reference genome (``GENOMEREF``): only ``HG38`` and ``HG19`` are supported

Cohort file names, as well as the fields mentioned above
must not contain dots.

The way to provide those values is through `OpenVariant <https://github.com/bbglab/openvariant>`__ , 
a comprehensive Python package that provides different functionalities to read, parse and operate 
different multiple input file formats (e. g. tsv, csv, vcf, maf, bed). 
Whether you are planning to run single or multiple cohorts, you would 
need to provide an annotation file in yaml format to specify the above mentioned structure required by IntOGen. 
Instructions on how to build an annotation file are documented here: `OpenVariant annotation file <https://openvariant.readthedocs.io/en/latest/user_guide/annotation_structure.html>`__ .

Optionally, DNA methylation and RNA-seq data of the cohorts can be provided
(``--methylation`` and ``--expression``). Each file must be named after the
cohort ID followed by a dot (e.g. :file:`TCGA_WXS_BRCA.beta.tsv.gz`).
See :doc:`omics` for the accepted formats.


Output
******

By default this pipeline outputs 4 files:

- :file:`cohorts.tsv`: summary of the cohorts that have been analyzed
- :file:`drivers.tsv`: summary of the results of the driver discovery by cohort
- :file:`mutations.tsv`: summary of all the mutations analyzed by cohort
- :file:`unique_drivers.tsv`: information on the genes reported as drivers (in any cohort)
- :file:`unfiltered_drivers.tsv`: information on the filters applied to the post-processing step: from the output of the combination to the final set of driver genes.
- :file:`omics.tsv`: only when methylation or RNA-seq data are provided. Genes with
  significant epigenetic silencing or expression outliers in any cohort.
  In addition, omics features are added to :file:`drivers.tsv` and
  :file:`unfiltered_drivers.tsv` (see :doc:`omics`).
- :file:`pathways.tsv`, :file:`pathway_cooccurrence.tsv` and :file:`pathway_modules.tsv`:
  only with ``--pathways true``. Gene sets under selection (with and without their
  individually significant genes), co-occurring or mutually exclusive
  dysregulation events and the modules of co-occurring events. In addition, the significant gene sets and modules of each driver
  are added to :file:`drivers.tsv` (see :doc:`pathways`).

Those files can be found in the path indicated with the
``--output`` options.

Moreover, the ``--debug true`` options will generate a
:file:`debug` folder under the output folder, in which
all the input and output files of the different methods are
linked.
