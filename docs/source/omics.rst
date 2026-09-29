DNA methylation and transcriptomics
-----------------------------------

The driver discovery methods of IntOGen detect genes under positive
selection from the somatic point mutations of a cohort. DNA methylation
(e.g. Illumina array beta values) and transcriptomic (RNA-seq) data of the
same tumours can be provided as optional inputs. They are not evidence of
positive selection by themselves, but they inform on:

- whether a gene is **expressed** in the tumours of the cohort: a mutated gene
  that is not expressed is unlikely to be a driver.
- **alternative mechanisms of alteration** of cancer genes: tumour suppressors
  are frequently silenced by promoter hypermethylation (e.g. *MLH1*,
  *CDKN2A*, *MGMT*) and oncogenes can be activated by over-expression
  (e.g. amplifications, fusions).
- the **functional consequences of the mutations** and the **mode of action**
  of the genes (e.g. loss of expression of mutated alleles, higher expression
  of mutated oncogenes).

The omics data are integrated at three levels. The first two are applied
whenever omics data are provided; the third one is optional:

1. **Expression filter**: the expression of the cohort replaces the TCGA-based
   proxy used to flag non-expressed candidate drivers
   (see :doc:`postprocessing`).
2. **Omics features**: every driver is annotated with its epigenetic silencing,
   expression outliers and the expression of mutated vs wild-type tumours.
   Genes altered by these non-mutational mechanisms are reported for all
   cohorts in :file:`omics.tsv`.
3. **Integrative mode** (``--integrate_omics true``): the epigenetic silencing
   and expression outlier tests are added as additional voters to the
   combination of driver identification methods
   (see :doc:`drivers_combination`).

|image1|


Input
^^^^^

Omics data are provided as matrices (features x samples) in TSV or CSV files
(optionally gzipped). The first column contains the feature identifiers and
the remaining columns the samples. Non-numeric columns (e.g. a gene
description) are ignored.

Each file is associated to a cohort by its name: the cohort ID (``DATASET``)
followed by a dot, e.g. :file:`TCGA_WXS_BRCA.beta.tsv.gz` or
:file:`TCGA_WXS_BRCA.counts.tsv.gz`. Omics files of cohorts without mutation
data are processed but not integrated.

.. code-block:: bash

      nextflow run intogen.nf --input <input> \
            --methylation "omics/*.beta.tsv.gz" \
            --expression "omics/*.counts.tsv.gz" \
            --omics_samples omics/samples.tsv


--methylation <files>   DNA methylation matrices.

--expression <files>   RNA-seq matrices.

--omics_samples <file>   Optional sample sheet (see below).

--methylation_values <type>   ``auto`` (default), ``beta``, ``m`` (M-values) or ``percent``.
   When ``auto``, values in [0, 1] are considered beta values, values in [0, 100]
   percentages and any other range M-values.

--methylation_probes <file>   Promoter probes of the methylation arrays.
   Default: :file:`<datasets>/methylation/promoter_probes.tsv.gz`.

--expression_units <units>   ``auto`` (default), ``counts``, ``tpm``, ``fpkm`` or ``log2``.
   When ``auto``, integer values (or values of samples with more than 5 million
   reads, e.g. RSEM expected counts) are considered counts, matrices with
   negative values or values below 30 log-transformed values, and TPM/FPKM otherwise.

--integrate_omics <bool>   Add the omics evidence to the combination of methods. Default: ``false``.


DNA methylation
***************

Two types of matrices are accepted:

- **Probe level**: Illumina HumanMethylation450, EPIC or EPICv2 probes
  (e.g. ``cg00000029`` or ``cg00000029_TC21``). The methylation of the
  promoter of each gene is the mean beta value of the probes whose CpG is
  located within 1,500 bp upstream and 500 bp downstream of the transcription
  start site of its MANE transcript (the transcripts analysed by IntOGen).
  Probes flagged by the general quality mask of Zhou et al. [1]_ (e.g. probes
  overlapping common SNPs or mapping to multiple locations) are discarded.
- **Gene level**: promoter methylation per gene (HUGO symbols or Ensembl
  gene IDs).

Genes in the sex chromosomes are discarded, since their promoter methylation
reflects the X inactivation in female patients.

RNA-seq
*******

Genes are identified by HUGO symbols, Ensembl gene IDs (with or without
version) or ``SYMBOL|ENTREZ`` pairs (as in TCGA legacy matrices). Outdated
symbols are updated and only MANE protein coding genes are kept. When a gene
appears more than once, the entry with the highest mean expression is used.

Counts are normalised to counts per million (CPM). Expression values are
analysed as log2(x + 1), where x is the CPM, TPM or FPKM.

Samples
*******

By default, every column of a matrix is a tumour sample whose ID is the one
used in the mutation data (``SAMPLE``). A sample sheet (``--omics_samples``)
can be used to map the column IDs to the samples of the mutation data and to
identify normal samples:

.. code-block:: text

      ID                      SAMPLE              TYPE
      TCGA-A1-A0SB-01A-11D    TCGA-A1-A0SB-01A    tumor
      TCGA-A1-A0SB-11A-41D                        normal

- ``ID``: column of the omics matrix (required)
- ``SAMPLE``: sample in the mutation data (default: ``ID``)
- ``TYPE``: ``tumor`` (default) or ``normal``

As for the mutations, only one column per sample is used (the first one in
alphabetical order). Normal samples are used as the reference of the promoter
methylation and ignored in the RNA-seq data.

.. important:: The omics layers are analysed using all their tumour samples,
   but the features relating omics and mutations (e.g. expression of
   mutated tumours) only use the tumours with both types of data
   (i.e. matching sample IDs, after the filters applied to the mutations).


Epigenetic silencing
^^^^^^^^^^^^^^^^^^^^

For each gene:

1. **Reference state**: the promoter methylation in normal tissue is the mean
   beta value of the normal samples (at least 3). When not enough normal samples
   are available, it is estimated as the lower quartile of the beta values of the
   tumours, so that genes hypermethylated in up to ~75% of the tumours (e.g.
   *GSTP1* in prostate cancer) can be detected. Only genes with an unmethylated
   reference promoter (beta < 0.2) and at least 10 tumours with data are tested.
2. **Hypermethylation**: a tumour is hypermethylated at a gene if its promoter
   beta value is at least 0.3 and at least 0.2 above the reference. These
   thresholds follow the criteria commonly used to define epigenetic silencing in
   TCGA studies.
3. **Recurrence**: the number of hypermethylated tumours is compared with the
   number expected given the hypermethylation rate of each tumour
   (the fraction of the tested genes that are hypermethylated in it). Under
   the null hypothesis, the number of hypermethylated tumours follows a
   Poisson-binomial distribution. Similarly to how hypermutated samples are
   handled in the mutational analysis, tumours with a genome-wide
   hypermethylation (e.g. CpG island methylator phenotype) contribute little
   evidence.
4. **Functional silencing**: when expression data are available, the
   expression of the hypermethylated tumours is compared with that of the other
   tumours (one-sided Mann-Whitney test; at least 3 tumours in each group).
   As in MethylMix [2]_, only genes whose methylation is associated with lower
   expression are considered functionally silenced: genes without evidence of
   silencing (p-value >= 0.05) get a p-value of 1.

P-values are corrected for multiple testing (Benjamini-Hochberg).


Expression
^^^^^^^^^^

**Expressed genes**: a gene is not expressed in the cohort when at least 80%
of the tumours have an expression below 1 (CPM/TPM/FPKM). This is the same
criterion applied to the TCGA data (see :doc:`postprocessing`).

**Expression outliers**: inspired by the Cancer Outlier Profile Analysis [3]_,
which identified recurrent gene fusions from outlier expression, an
expression outlier is a tumour with a modified z-score [4]_ above 3.5 (median
and median absolute deviation of the gene) and at least a 2-fold change with
respect to the median. Over-expression outliers must be expressed; under-expression
outliers can only be detected for expressed genes. As for the epigenetic
silencing, the number of tumours with an over- (or under-) expression outlier
is compared with the expectation given the outlier rate of each tumour
(Poisson-binomial test). The p-value of the gene is that of the most
significant direction (Bonferroni corrected for the two directions).
Genes in chromosome Y are not tested (their expression reflects the sex of
the patient).


Omics features
^^^^^^^^^^^^^^

For each gene of a cohort:

- **Expression of mutated tumours**: difference between the median
  expression of the tumours with non-synonymous mutations in the gene and the
  median expression of the tumours without them (two-sided Mann-Whitney test;
  at least 2 tumours in each group). The same comparison is done with the
  truncating mutations (nonsense, frameshift and splice site), whose lower
  expression usually reflects nonsense-mediated decay.
- **Mutations and hypermethylation**: number of tumours (with both mutation and
  methylation data) with a non-synonymous mutation and promoter
  hypermethylation of the gene, and with any of them. Tumour suppressors can be
  inactivated by either mechanism.
- **Mode of action supported by the omics data** (``OMICS_ROLE_SUPPORT``):
  ``LoF`` if the gene is functionally silenced by promoter hypermethylation or
  shows recurrent under-expression outliers; ``Act`` if it shows recurrent
  over-expression outliers or higher expression in mutated tumours
  (q-value < 0.1); ``conflicting`` if both. Lower expression of mutated tumours
  is not used as support for LoF, since nonsense-mediated decay also affects
  passenger truncating mutations. This annotation complements the mode of
  action inferred from the mutations (``ROLE``), which is not modified.


Integration in the driver discovery
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Expression filter
*****************

When the cohort has expression data, the genes flagged as not expressed in the
cohort are discarded as drivers (``Warning expression``), instead of those
flagged using TCGA data. Genes not measured in the cohort keep the flag
derived from TCGA. The source of the flag is reported
(``WARNING_EXPRESSION_SOURCE``: ``cohort`` or ``TCGA``).

Integrative mode
****************

With ``--integrate_omics true`` the epigenetic silencing and the expression
outliers tests are used as two additional driver identification methods in the
combination: their rankings take part in the Schulze voting and their p-values
in the weighted Stouffer combination (see :doc:`drivers_combination`). Their
weights are estimated for each cohort, as for the rest of the methods, based on
their ability to rank known cancer genes first. Thus, cohorts where the omics
evidence does not help to prioritise cancer genes give it low weights.

Since the aim of IntOGen is the discovery of driver genes with positive
selection of mutations:

- the omics methods only rank the candidate genes of the mutational analysis
  (genes with at least two mutated samples). Genes altered only by
  non-mutational mechanisms are not part of the combination (they are reported
  in :file:`omics.tsv`). The multiple testing correction of the combination
  is not affected.
- non-CGC genes still need at least two significant methods, one of which must
  be a mutation-based method.

.. note:: In integrative mode, the combined q-values reflect both the
   evidence of positive selection of mutations and the omics evidence.


Output
^^^^^^

:file:`omics.tsv` contains the genes with significant (q-value < 0.1)
epigenetic silencing or expression outliers in any cohort, with the columns
described below.

When omics data are available, the following columns are added to
:file:`drivers.tsv`:

.. list-table::
   :header-rows: 1

   * - Column
     - Description
   * - EXPRESSED
     - Whether the gene is expressed in the cohort
   * - EXPRESSION_OUTLIER_DIRECTION
     - Direction (over/under) of the expression outliers
   * - EXPRESSION_OUTLIER_SAMPLES
     - Tumours with an expression outlier in that direction
   * - QVALUE_EXPRESSION
     - Q-value of the recurrence of expression outliers
   * - EXPRESSION_LOG2FC_MUTANT
     - Log2 fold change of the expression of mutated vs wild-type tumours
   * - QVALUE_EXPRESSION_MUTANT
     - Q-value of the comparison of mutated vs wild-type tumours
   * - EXPRESSION_LOG2FC_TRUNCATING
     - Log2 fold change of the tumours with truncating mutations
   * - PVALUE_EXPRESSION_TRUNCATING
     - P-value of the comparison of tumours with truncating mutations
   * - METHYLATION_HYPER_SAMPLES
     - Tumours with promoter hypermethylation
   * - METHYLATION_FUNCTIONAL
     - Whether hypermethylation is associated with lower expression (empty without expression)
   * - QVALUE_METHYLATION
     - Q-value of the epigenetic silencing test
   * - SAMPLES_MUTATED_AND_HYPERMETHYLATED
     - Tumours with both a mutation and promoter hypermethylation of the gene
   * - SAMPLES_MUTATED_OR_HYPERMETHYLATED
     - Tumours with a mutation or promoter hypermethylation of the gene
   * - OMICS_ROLE_SUPPORT
     - Mode of action supported by the omics data (LoF, Act or conflicting)

:file:`unfiltered_drivers.tsv` includes ``WARNING_EXPRESSION_SOURCE``,
``QVALUE_METHYLATION`` and ``QVALUE_EXPRESSION``.

The complete results of each cohort are in the :file:`steps` folder:
:file:`methylation` (promoter methylation matrices, epigenetic silencing tests
and hypermethylation events), :file:`expression` (normalised matrices,
expression tests and outlier events) and :file:`omics` (omics features of all
genes).


Considerations
^^^^^^^^^^^^^^

- Tumour purity and the infiltration of immune and stromal cells affect both
  methylation and expression. Expression outliers of genes specific of
  non-tumour cells are expected in some cohorts.
- The recurrence tests account for the differences between tumours but not
  for gene-specific variability (e.g. genes that are frequently hypermethylated
  in cancer without consequences, like many Polycomb targets). Expression data
  and normal samples improve the specificity of the epigenetic silencing
  analysis: without expression, recurrent hypermethylation cannot be told apart
  from silencing.
- Batch effects between samples should be corrected before running the
  pipeline.
- Omics data from all the tumours of a cohort are used, independently of the
  filters applied to the mutations (e.g. hypermutated samples). Multiple samples
  of the same patient should be avoided (e.g. using the sample sheet).


.. [1] Zhou W, Laird PW, Shen H. Comprehensive characterization, annotation and innovative use of Infinium DNA methylation BeadChip probes. Nucleic Acids Res. 2017;45(4):e22. doi:10.1093/nar/gkw967

.. [2] Gevaert O. MethylMix: an R package for identifying DNA methylation-driven genes. Bioinformatics. 2015;31(11):1839-1841. doi:10.1093/bioinformatics/btv020

.. [3] Tomlins SA, et al. Recurrent fusion of TMPRSS2 and ETS transcription factor genes in prostate cancer. Science. 2005;310(5748):644-648. doi:10.1126/science.1117679

.. [4] Iglewicz B, Hoaglin DC. How to detect and handle outliers. ASQC Quality Press; 1993.

.. |image1| image:: /_static/omics_integration.svg
   :width: 7in
   :align: middle
