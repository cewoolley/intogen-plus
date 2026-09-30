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
2. **Co-occurrence and mutual exclusivity of dysregulation events**: whether
   the drivers, the genes with significant epigenetic silencing or expression
   outliers and the selected long tails of the gene sets are altered in the
   same tumours more (co-occurrence) or less (mutual exclusivity, e.g.
   alternative alterations of a pathway) often than expected by chance. Events
   that co-occur are grouped into modules: clusters of dysregulation that
   happen together.

The results are reported separately from the driver genes, which are not
modified. Only a reference to the significant gene sets and modules of each
driver is added to :file:`drivers.tsv`.

.. code-block:: bash

      nextflow run intogen.nf --input <input> --pathways true

The analysis only requires the somatic mutations: from sequencing data alone,
it finds gene sets under selection through missense or truncating mutations
and the networks of drivers and gene sets mutated together in the same
tumours. Methylation and RNA-seq data of the cohorts, when provided, add the
epigenetic silencing and expression layers (see :doc:`omics`):

.. code-block:: bash

      nextflow run intogen.nf --input <input> --pathways true \
            --methylation "omics/*.beta.tsv.gz" \
            --expression "omics/*.counts.tsv.gz"


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
     - Mutations that alter the protein
     - Mutations (selection: dNdScv)
   * - mutation_missense
     - Missense mutations
     - Mutations (selection: dNdScv)
   * - mutation_truncating
     - Nonsense, frameshift and essential splice site mutations
     - Mutations (selection: dNdScv)
   * - silencing
     - Promoter hypermethylation, excluding genes whose hypermethylation is not
       associated with lower expression
     - Methylation (see :doc:`omics`)
   * - expression_over / expression_under
     - Expression outliers
     - RNA-seq (see :doc:`omics`)

The selection of the mutation layers is tested with the substitutions analysed
by dNdScv (missense, nonsense and essential splice site). The co-occurrence
analyses use all the mutations that alter the protein (VEP consequences of high
and moderate impact, including indels).


Selection of gene sets
^^^^^^^^^^^^^^^^^^^^^^

**Mutations.** The test extends the model of dNdScv [3]_ from genes to sets
of genes. dNdScv estimates, for each gene, the number of synonymous and
non-synonymous substitutions expected under neutrality given the mutational
profile of the cohort, the sequence of the gene and its local mutation rate,
which is predicted from genomic covariates with gamma-distributed variation
(overdispersion :math:`\theta`, estimated from the synonymous mutations of
all genes). Since sets aggregate many genes, their tests are sensitive to
errors in :math:`\theta`: the lower bound of its 95% profile likelihood
confidence interval is used, which makes the tests conservative when few genes
are analysed and is close to the estimate for a whole exome. Updated with the synonymous mutations observed in each gene, the
number of non-synonymous mutations of a gene under neutrality follows a
negative binomial distribution. The number of non-synonymous mutations of the
set is compared with the exact distribution of the sum of these negative
binomials (one-sided test). ``RATIO`` is the observed/expected number of
non-synonymous mutations, i.e. the dN/dS of the set.

**Background dN/dS.** In real exomes, genes that are not under positive
selection do not have a dN/dS of exactly 1 under the global mutational model of
the cohort, and the deviation depends on their genomic context: genes that are
not expressed in the tissue (e.g. neuronal genes, ion channels, muscle and
extracellular matrix genes) escape purifying selection and
transcription-coupled repair and have dN/dS around 1.1-1.2. Summed over
hundreds of genes, these small deviations make large gene sets significant.
In TCGA exomes, testing against the neutral model (dN/dS = 1) reports tens of
such gene sets per cohort, with random gene sets rejected up to 3 times more
often than expected (see :ref:`case-study`). Gene sets are therefore compared
with genes of similar context (competitive null [7]_, see also [8]_ on the
heterogeneity of mutation rates): the expected non-synonymous
mutations of each gene are multiplied by the background dN/dS of its context,
predicted with a Poisson regression of the non-synonymous mutations of the
genes that are not individually significant (with their neutral expectation as
offset) on the epigenomic covariates used by dNdScv (20 principal components of
chromatin marks across tissues; :file:`<datasets>/pathways/gene_covariates.tsv.gz`).
Missense and truncating mutations have their own background. Without the
covariates file, the covariate-predicted mutation rate and the size of the
genes are used, which is calibrated but removes only part of the confounding.

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

**Combination.** The mutation p-values of a set (``mutation``,
``mutation_missense`` and ``mutation_truncating``) are combined with Simes'
method, so that sets selected through a single type of mutation (e.g. truncating
mutations of tumour suppressors) are not diluted by the passenger mutations of
the other types. Then, the mutation, silencing and expression (the most
significant direction, Bonferroni corrected) p-values are combined with
Fisher's method (``combined`` layer). Without omics data, the combination only
reflects the mutations.

P-values are corrected for multiple testing with the Benjamini-Hochberg
method, separately for each scope and layer.


Co-occurrence of dysregulation events
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Events.** An event is the alteration of a gene or of a gene set in a tumour,
in one of the layers above. The candidate events of a cohort are:

- the drivers (``mutation``)
- the genes with significant (q-value < 0.1) epigenetic silencing or expression outliers
- the long tails with significant selection (q-value < 0.1) of the gene sets,
  in the layer where they are most enriched (highest observed/expected ratio)
  among those where they are significant (one mutation layer and one per omics
  layer). A gene set event is the alteration of any of the genes of its long
  tail in that layer, e.g. the truncating mutations of the gene sets selected
  through truncating mutations. Using the alterations under selection reduces
  the dilution by passenger mutations, which is important to detect networks
  from mutations alone.

With sequencing data only, the events are the drivers and the gene sets
selected through their mutations, and the modules are networks of drivers and
pathways mutated together.

The 100 most significant candidates altered in at least 3 tumours are tested,
in pairs, in the tumours with data of both layers.

**Null model.** Tumours with many alterations (e.g. hypermutated or smoking
tumours, or tumours with a CpG island methylator phenotype) have many events,
so any pair of passenger-rich events co-occurs in them: in TCGA cohorts,
pairwise Fisher's exact tests on the most mutated genes (as commonly done)
report hundreds of co-occurring pairs, most of them between passenger genes.
However, the dependence on the burden differs between events: the selected
mutations of drivers barely depend on it, or even decrease with it (e.g. EGFR
mutations in lung adenocarcinomas of non-smokers). Models that impose the same
dependence on every gene, such as the additive model of DISCOVER [4]_, remove
the passenger artefacts but also the power to detect co-occurrence of drivers.
Each event therefore has its own burden elasticity: the probability that it
happens in tumour :math:`s` is estimated with a logistic regression on the
alteration burden of the tumour in the layer of the event (number of altered
genes, without the genes of the event):

.. math::

   q_s = \frac{1}{1 + e^{-(a + b \log(1 + \text{burden}_s))}}

Under independence given the burden, the number of tumours with both events
follows a Poisson-binomial distribution with probabilities
:math:`q_{1s} q_{2s}`. Its upper tail gives the co-occurrence p-value and its
lower tail the mutual exclusivity p-value.

**Overlapping events.** The genes shared by two events of the same layer (e.g.
a driver and a pathway it belongs to, or two overlapping pathways; the mutation
layers count as one) are removed from both events before testing them; pairs
left without genes are not tested.
The same applies to expression events, since the expression of a gene is often
a direct consequence of its other alterations. Mutations and silencing of the
same gene (two hits) are tested.

**Multiple testing.** Rare events cannot co-occur significantly with any other
event: the smallest p-value that a pair can achieve depends on the number of
tumours altered by each event. Tests that cannot reach significance are
discarded before the Benjamini-Hochberg correction (Tarone's procedure [5]_
[6]_), which increases the power to detect the co-occurrence of rare events.
Pairs discarded are reported with ``TESTABLE`` false and q-value 1. Co-occurrence
and mutual exclusivity are corrected separately.

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

:file:`pathway_cooccurrence.tsv` contains the pairs of events that co-occur or are
mutually exclusive significantly (q-value < 0.1):

.. list-table::
   :header-rows: 1

   * - Column
     - Description
   * - COHORT
     - Cohort
   * - EVENT_1, NAME_1, EVENT_2, NAME_2
     - Events (``<layer>:<gene or gene set>``, e.g. ``mutation:TP53`` or
       ``mutation_truncating:R-HSA-1234``) and their names
   * - TUMOURS
     - Tumours with data of both layers
   * - TUMOURS_1, TUMOURS_2, TUMOURS_BOTH
     - Tumours with each event and with both
   * - EXPECTED_BOTH, RATIO
     - Tumours expected with both events and observed/expected ratio
   * - P_VALUE, Q_VALUE
     - Co-occurrence p-value and q-value
   * - P_VALUE_EXCLUSIVITY, Q_VALUE_EXCLUSIVITY
     - Mutual exclusivity p-value and q-value
   * - SHARED_GENES_REMOVED
     - Number of genes shared by the events and removed before testing them
   * - MODULE
     - Module of the events
   * - TESTABLE, TESTABLE_EXCLUSIVITY
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
mutation test (``mutation_theta``, the value used, and ``mutation_theta_mle``,
the maximum likelihood estimate), background dN/dS (``mutation_background``:
covariates used and quantiles of the background dN/dS), gene sets tested and
number of significant results).


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
- The selection of the mutation layers is tested with the substitutions
  analysed by dNdScv (indels are not included). When dNdScv has no results for
  a cohort, only the omics layers are analysed.
- Single cohorts have limited power: in TCGA cohorts (370-1000 tumours),
  known co-occurrences such as STK11-KEAP1 in lung adenocarcinoma do not survive
  the correction for all the pairs tested, while mutual exclusivity (e.g.
  KRAS-EGFR, BRAF-KRAS, TP53-CDH1) is a stronger signal. Absence of
  co-occurrence is not evidence of independence.
- The background dN/dS reduces, but may not remove, the confounding by the
  genomic context of the genes: large sets of genes not expressed in the tissue
  (e.g. neuronal genes) with modest ratios should be interpreted with care.
  Genes not expressed in the cohort (with RNA-seq data) are a direct way to
  check it.
- Co-occurrence does not imply cooperation: events can co-occur because they
  are associated with a common factor not captured by the tumour burden (e.g.
  subtypes within a cohort, tumour purity for the omics layers).
- The candidate events are limited to the most significant ones to keep the
  number of tests manageable. The analysis can be run for a cohort with other
  parameters with ``pathway-analysis`` (``--min-size``, ``--max-size``,
  ``--threshold``, ``--min-tumours`` and ``--max-events``).


.. _case-study:

Evaluation on TCGA exomes
^^^^^^^^^^^^^^^^^^^^^^^^^

The analysis was evaluated with the TCGA MC3 somatic mutations of four cohorts
(KIRC, 369 tumours; LUAD, 515; BRCA, 1016; COAD, 391), dNdScv 0.0.1.0 (commit 43c5e2f) and the
Reactome and hallmark gene sets of MSigDB 7.5.1. The scripts are in
:file:`benchmarks/tcga_case_study`.

**Calibration.** Fraction of 400 random sets of genes that are not
individually significant with p-value < 0.05 (expected 0.05):

.. list-table::
   :header-rows: 1

   * - Test
     - KIRC
     - LUAD
     - BRCA
     - COAD
   * - Neutral expectation, Poisson
     - 0.065
     - 0.235
     - 0.138
     - 0.220
   * - dNdScv ``genesetdnds`` (one-sided)
     - 0.052
     - 0.075
     - 0.113
     - 0.125
   * - Negative binomial, neutral model
     - 0.028
     - 0.048
     - 0.045
     - 0.100
   * - Negative binomial, background dN/dS (default)
     - 0.022
     - 0.030
     - 0.025
     - 0.040

**Power.** Selection spread over the long tail of a set (40 random genes that
are not individually significant, 25 replicates), detected with q-value < 0.1
among all the gene sets. The genes are rarely significant individually, so the
over-representation of individually significant genes (gene-level driver
discovery followed by enrichment analysis) detects none of them:

.. list-table::
   :header-rows: 1

   * - Selection
     - Gene level + enrichment
     - Background dN/dS (KIRC / LUAD / BRCA / COAD)
   * - Missense, dN/dS 1.6
     - 0 / 0 / 0 / 0
     - 0 / 1.00 / 0.68 / 0.92
   * - Missense, dN/dS 2
     - 0 / 0 / 0 / 0.04
     - 0.60 / 1.00 / 1.00 / 1.00
   * - Truncating, dN/dS 3
     - 0 / 0 / 0 / 0.04
     - 0.08 / 1.00 / 0.84 / 0.96

**Co-occurrence.** Pairwise Fisher's exact tests on the most mutated genes
and drivers report 0, 400, 87 and 714 co-occurring pairs (KIRC, LUAD, BRCA,
COAD; q-value < 0.1), most of them with passenger genes. The burden-aware test
reports none of those, and recovers known mutual exclusivity (KRAS-EGFR in LUAD;
TP53 with CDH1, GATA3, PIK3CA, FOXA1 and MAP3K1, and PIK3CA-AKT1 in BRCA; BRAF-KRAS
and APC-BRAF in COAD) and the co-occurrence of CDH1 and ERBB2 mutations in
breast cancer. With co-occurrence planted between two genes of real tumours,
the burden-aware test keeps its false positive rate below 0.03 whatever the
dependence of the genes on the burden (Fisher's test: up to 1.0), with power
close to Fisher's test for driver-like genes (the additive DISCOVER model has
almost none in COAD).


.. [1] Milacic M, et al. The Reactome Pathway Knowledgebase 2024. Nucleic Acids Res. 2024;52(D1):D672-D678. doi:10.1093/nar/gkad1025

.. [2] Liberzon A, et al. The Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Syst. 2015;1(6):417-425. doi:10.1016/j.cels.2015.12.004

.. [3] Martincorena I, et al. Universal patterns of selection in cancer and somatic tissues. Cell. 2017;171(5):1029-1041. doi:10.1016/j.cell.2017.09.042

.. [4] Canisius S, Martens JWM, Wessels LFA. A novel independence test for somatic alterations in cancer shows that biology drives mutual exclusivity but chance explains most co-occurrence. Genome Biol. 2016;17:261. doi:10.1186/s13059-016-1114-x

.. [5] Tarone RE. A modified Bonferroni method for discrete data. Biometrics. 1990;46(2):515-522. doi:10.2307/2531456

.. [7] Goeman JJ, Bühlmann P. Analyzing gene expression data in terms of gene sets: methodological issues. Bioinformatics. 2007;23(8):980-987. doi:10.1093/bioinformatics/btm051

.. [8] Lawrence MS, et al. Mutational heterogeneity in cancer and the search for new cancer-associated genes. Nature. 2013;499(7457):214-218. doi:10.1038/nature12213

.. [6] Gilbert PB. A modified false discovery rate multiple-comparisons procedure for discrete data, applied to human immunodeficiency virus genetics. J R Stat Soc Ser C Appl Stat. 2005;54(1):143-158. doi:10.1111/j.1467-9876.2005.00475.x
