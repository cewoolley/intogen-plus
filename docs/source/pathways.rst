Pathways and co-occurrence of dysregulation
-------------------------------------------

The driver discovery methods of IntOGen test each gene separately. Many
pathways, however, are altered in cancer through many different genes, each
one altered in few tumours: a pathway can be under positive selection even if
none of its rarely altered genes (its *long tail*) reaches significance by
itself. Moreover, tumours accumulate alterations in several pathways and
through several mechanisms (mutations, epigenetic silencing, expression
changes) that may be selected together.

The optional pathway analysis (``--pathways true``) addresses both questions
for each cohort:

1. **Selection of gene sets**: whether the genes of a pathway or hallmark
   accumulate more mutations, epigenetic silencing or expression outliers than
   expected, with and without their individually significant genes.
2. **Co-occurrence of dysregulation events**: whether the drivers, the genes
   with significant epigenetic silencing or expression outliers and the
   selected long tails of the gene sets are altered in the same tumours more
   often than expected by chance. Events that co-occur are grouped into
   modules: clusters of dysregulation that happen together.

The results are reported separately from the driver genes, which are not
modified. Only a reference to the significant gene sets and modules of each
driver is added to :file:`drivers.tsv`.

.. code-block:: bash

      nextflow run intogen.nf --input <input> --pathways true \
            --methylation "omics/*.beta.tsv.gz" \
            --expression "omics/*.counts.tsv.gz"

The mutation layers are always analysed. The epigenetic silencing and
expression layers are added when methylation and RNA-seq data of the cohort
are provided (see :doc:`omics`).


Gene sets
^^^^^^^^^

By default, the gene sets are the pathways of
`Reactome <https://reactome.org>`_ [1]_ and the hallmark gene sets of
`MSigDB <https://www.gsea-msigdb.org>`_ [2]_, built with the rest of the
datasets (:file:`<datasets>/pathways/gene_sets.tsv.gz`). Gene symbols are
harmonised with the genes analysed by IntOGen (outdated symbols are updated
and other genes are discarded).

--pathways <bool>   Run the pathway analysis. Default: ``false``.

--gene_sets <file>   Custom gene sets, either a TSV file with the columns
   ``SET``, ``SOURCE``, ``NAME`` and ``SYMBOL`` (one row per gene) or a
   GMT file. Default: :file:`<datasets>/pathways/gene_sets.tsv.gz`.

Only the gene sets with 10 to 500 genes among the genes analysed in the
cohort are tested.


Layers of alteration
^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1

   * - Layer
     - Alterations
     - Data
   * - mutation
     - Non-synonymous substitutions (missense, nonsense and essential splice site)
     - dNdScv
   * - mutation_missense
     - Missense substitutions
     - dNdScv
   * - mutation_truncating
     - Nonsense and essential splice site substitutions
     - dNdScv
   * - silencing
     - Promoter hypermethylation, excluding genes whose hypermethylation is not
       associated with lower expression
     - Methylation (see :doc:`omics`)
   * - expression_over / expression_under
     - Expression outliers
     - RNA-seq (see :doc:`omics`)


Selection of gene sets
^^^^^^^^^^^^^^^^^^^^^^

**Mutations.** The test extends the model of dNdScv [3]_ from genes to sets
of genes. dNdScv estimates, for each gene, the number of synonymous and
non-synonymous substitutions expected under neutrality given the mutational
profile of the cohort, the sequence of the gene and its local mutation rate,
which is predicted from genomic covariates with gamma-distributed variation
(overdispersion :math:`\theta`, estimated from the synonymous mutations of
all genes). Updated with the synonymous mutations observed in each gene, the
number of non-synonymous mutations of a gene under neutrality follows a
negative binomial distribution. The number of non-synonymous mutations of the
set is compared with the exact distribution of the sum of these negative
binomials (one-sided test). ``RATIO`` is the observed/expected number of
non-synonymous mutations, i.e. the dN/dS of the set.

Without the covariates of dNdScv, a conditional test is used instead:
given the number of mutations of a gene, its number of non-synonymous
mutations is binomial under neutrality, with probability
:math:`e_{ns} / (e_{ns} + e_{s})`, so that the unknown local mutation rate
cancels out. This test is robust but has less power.

**Epigenetic silencing and expression outliers.** The number of events (e.g.
hypermethylated tumours) of the genes of the set is compared with the number
expected given the event rate of each tumour (the sum of the expectations of
the gene-level tests described in :doc:`omics`), using a Poisson distribution.

**Scopes.** Each set is tested with all its genes (``all``) and without the
genes that are individually significant in the layer (``long_tail``). For the
mutation layers, the genes with tier 1 to 3 in the combination of driver
discovery methods (see :doc:`drivers_combination`) are excluded; for the
omics layers, the genes with q-value < 0.1. A significant long tail means that
the set is under selection beyond its known drivers, through genes that are
rarely altered individually.

**Combination.** The mutation, silencing and expression (the most significant
direction, Bonferroni corrected) p-values of a set are combined with Fisher's
method (``combined`` layer).

P-values are corrected for multiple testing with the Benjamini-Hochberg
method, separately for each scope and layer.


Co-occurrence of dysregulation events
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Events.** An event is the alteration of a gene or of a gene set in a tumour,
in one of the layers ``mutation``, ``silencing``, ``expression_over`` and
``expression_under``. The candidate events of a cohort are:

- the drivers (mutations)
- the genes with significant (q-value < 0.1) epigenetic silencing or expression outliers
- the long tails with significant selection (q-value < 0.1) of the gene sets,
  in the layer where they are most significant. A gene set event is the
  alteration of any of the genes of its long tail.

The 100 most significant candidates altered in at least 3 tumours are tested,
in pairs, in the tumours with data of both layers.

**Background model.** Tumours with many alterations (e.g. hypermutated
tumours or tumours with a CpG island methylator phenotype) have many events,
and any pair of events tends to co-occur in them. As in DISCOVER [4]_, the
probability :math:`p_{gs}` that gene :math:`g` is altered in tumour :math:`s`
is estimated, for each layer, with a model that keeps both the number of
tumours altered in each gene and the number of genes altered in each tumour:

.. math::

   p_{gs} = \frac{1}{1 + e^{-(a_g + b_s)}}

The probability that an event made of several genes happens in a tumour is
:math:`q_s = 1 - \prod_g (1 - p_{gs})`. Under independence, the number of
tumours with both events follows a Poisson-binomial distribution with
probabilities :math:`q_{1s} q_{2s}`, which gives the p-value of the
co-occurrence (one-sided test).

**Overlapping events.** The genes shared by two events of the same layer (e.g.
a driver and a pathway it belongs to, or two overlapping pathways) are removed
from both events before testing them; pairs left without genes are not tested.
The same applies to expression events, since the expression of a gene is often
a direct consequence of its other alterations. Mutations and silencing of the
same gene (two hits) are tested.

**Multiple testing.** Rare events cannot co-occur significantly with any other
event: the smallest p-value that a pair can achieve depends on the number of
tumours altered by each event. Tests that cannot reach significance are
discarded before the Benjamini-Hochberg correction (Tarone's procedure [5]_
[6]_), which increases the power to detect the co-occurrence of rare events.
Pairs discarded are reported with ``TESTABLE`` false and q-value 1.

**Modules.** The events connected by significant co-occurrences
(q-value < 0.1) are grouped into modules (connected components).


Output
^^^^^^

:file:`pathways.tsv` contains the gene sets with significant selection
(q-value < 0.1) in any scope or layer, one row per cohort and set:

.. list-table::
   :header-rows: 1

   * - Column
     - Description
   * - COHORT, SET, SOURCE, NAME
     - Cohort and gene set (identifier, source and name)
   * - N_GENES
     - Genes of the set analysed in the cohort
   * - <LAYER>_Q, <LAYER>_RATIO
     - Q-value and observed/expected ratio of the set in the layer (all genes)
   * - <LAYER>_LONG_TAIL_Q, <LAYER>_LONG_TAIL_RATIO
     - The same, without the individually significant genes
   * - MUTATION_LONG_TAIL_TOP_GENES
     - Genes of the long tail with most non-synonymous mutations (gene:mutations)

The layers are ``COMBINED``, ``MUTATION``, ``MUTATION_MISSENSE``,
``MUTATION_TRUNCATING`` and, with omics data, ``SILENCING``,
``EXPRESSION_OVER`` and ``EXPRESSION_UNDER``.

:file:`pathway_cooccurrence.tsv` contains the pairs of events that co-occur
significantly (q-value < 0.1):

.. list-table::
   :header-rows: 1

   * - Column
     - Description
   * - COHORT
     - Cohort
   * - EVENT_1, NAME_1, EVENT_2, NAME_2
     - Events (``<layer>:<gene or gene set>``) and their names
   * - TUMOURS
     - Tumours with data of both layers
   * - TUMOURS_1, TUMOURS_2, TUMOURS_BOTH
     - Tumours with each event and with both
   * - EXPECTED_BOTH, RATIO
     - Tumours expected with both events and observed/expected ratio
   * - P_VALUE, Q_VALUE
     - Co-occurrence p-value and q-value
   * - SHARED_GENES_REMOVED
     - Number of genes shared by the events and removed before testing them
   * - MODULE
     - Module of the events
   * - TESTABLE
     - Whether the pair could reach significance (see above)

:file:`pathway_modules.tsv` lists the events of each module (``COHORT``,
``MODULE``, ``N_EVENTS``, ``N_PAIRS``, ``EVENTS``).

The following columns are added to :file:`drivers.tsv`:

.. list-table::
   :header-rows: 1

   * - Column
     - Description
   * - PATHWAYS
     - Gene sets of the driver whose long tail is significant (combined layer)
   * - MODULES
     - Co-occurrence modules with an event of the driver

The complete results of each cohort are in the :file:`steps/pathways` folder:
the tests of all the gene sets, scopes and layers
(:file:`<COHORT>.pathways.tsv.gz`, including the observed and expected
alterations, the number of altered tumours and the most altered genes),
all the pairs of events tested (:file:`<COHORT>.pathway_cooccurrence.tsv.gz`),
the significant gene sets and modules of every gene
(:file:`<COHORT>.pathway_genes.tsv.gz`) and statistics of the analysis
(:file:`<COHORT>.pathways.tsv.gz.stats.json`: method and overdispersion of the
mutation test, gene sets tested and number of significant results).


Considerations
^^^^^^^^^^^^^^

- Gene sets overlap (e.g. the Reactome hierarchy), so that the results of
  related sets are correlated. Large sets accumulate the signal of many genes
  and can be significant with a small excess of alterations: the ratios, the
  number of altered genes and the most altered genes help to interpret them.
- The long tail excludes the genes that are significant in the combination of
  driver discovery methods (tiers 1 to 3), not only the final drivers, so
  that it is not driven by genes discarded by the filters of the
  post-processing.
- The mutation tests use the substitutions analysed by dNdScv (indels are not
  included). When dNdScv has no results for a cohort, only the omics layers are
  analysed.
- The co-occurrence test is conservative: the background model absorbs part of
  the dependence between events and rare events have little power. Absence of
  co-occurrence is not evidence of independence. Mutual exclusivity is not
  tested.
- Co-occurrence does not imply cooperation: events can co-occur because they
  are associated with a common factor not captured by the tumour burden (e.g.
  subtypes within a cohort, tumour purity for the omics layers).
- The candidate events are limited to the most significant ones to keep the
  number of tests manageable. The analysis can be run for a cohort with other
  parameters with ``pathway-analysis`` (``--min-size``, ``--max-size``,
  ``--threshold``, ``--min-tumours`` and ``--max-events``).


.. [1] Milacic M, et al. The Reactome Pathway Knowledgebase 2024. Nucleic Acids Res. 2024;52(D1):D672-D678. doi:10.1093/nar/gkad1025

.. [2] Liberzon A, et al. The Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Syst. 2015;1(6):417-425. doi:10.1016/j.cels.2015.12.004

.. [3] Martincorena I, et al. Universal patterns of selection in cancer and somatic tissues. Cell. 2017;171(5):1029-1041. doi:10.1016/j.cell.2017.09.042

.. [4] Canisius S, Martens JWM, Wessels LFA. A novel independence test for somatic alterations in cancer shows that biology drives mutual exclusivity but chance explains most co-occurrence. Genome Biol. 2016;17:261. doi:10.1186/s13059-016-1114-x

.. [5] Tarone RE. A modified Bonferroni method for discrete data. Biometrics. 1990;46(2):515-522. doi:10.2307/2531456

.. [6] Gilbert PB. A modified false discovery rate multiple-comparisons procedure for discrete data, applied to human immunodeficiency virus genetics. J R Stat Soc Ser C Appl Stat. 2005;54(1):143-158. doi:10.1111/j.1467-9876.2005.00475.x
