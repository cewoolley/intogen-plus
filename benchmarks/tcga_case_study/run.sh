#!/bin/bash
# Case study of the pathway analysis on TCGA exomes (see README.md)
set -euo pipefail
cd "$(dirname "$0")"
COHORTS="KIRC LUAD BRCA COAD"
PYTHON=${PYTHON:-python}

Rscript prepare.R data $COHORTS
$PYTHON inputs.py
for c in $COHORTS; do
    Rscript gsd.R data/$c/mats.rds results/$c/gsd_units.tsv results/$c/gsd.tsv > results/$c/gsd.log 2>&1 &
done
wait
$PYTHON corrected.py
$PYTHON calibration.py
$PYTHON spikein.py $COHORTS
$PYTHON compare_sets.py
$PYTHON cooccurrence.py
$PYTHON consolidate.py
