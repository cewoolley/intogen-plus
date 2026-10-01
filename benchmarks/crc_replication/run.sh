#!/bin/bash
# Pre-registered replication of the TCGA co-regulator findings (README.md, PREREGISTRATION.md)
#   [FEASIBILITY=1] ./run.sh <input dir with mutations.tsv.gz and samples.tsv> <hg19|hg38> <output dir> [Cornish et al. driver table]
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
IN=$1; BUILD=$2; OUT=$3; CORNISH=${4:-}
PYTHON=${PYTHON:-python}
WORK=$OUT/work
mkdir -p "$WORK"

Rscript -e 'cat(as.character(packageVersion("dndscv")))' > "$WORK/dndscv_version.txt"
Rscript "$HERE/dnds.R" regions "$BUILD" "$WORK/regions.bed"
$PYTHON "$HERE/coding_filter.py" "$IN/mutations.tsv.gz" "$WORK/regions.bed" "$WORK/coding.tsv.gz" "$WORK/sample_vaf.tsv"
# median allele fractions over all mutations from prepare_vcf.py --regions, when the input is already coding only
if [ -f "$IN/sample_vaf.tsv" ]; then cp "$IN/sample_vaf.tsv" "$WORK/sample_vaf.tsv"; fi
if $PYTHON "$HERE/classify.py" --check-status "$IN/samples.tsv"; then
  $PYTHON "$HERE/classify.py" "$IN/samples.tsv" "$WORK/coding.tsv.gz" "$WORK/sample_vaf.tsv" "$WORK/groups.tsv"
else
  Rscript "$HERE/dnds.R" annotate "$BUILD" "$WORK/coding.tsv.gz" "$WORK/annotate"
  $PYTHON "$HERE/classify.py" "$IN/samples.tsv" "$WORK/coding.tsv.gz" "$WORK/sample_vaf.tsv" "$WORK/groups.tsv" \
    "$WORK/annotate/annotmuts.tsv.gz"
fi
Rscript "$HERE/dnds.R" strata "$BUILD" "$WORK/coding.tsv.gz" "$WORK/groups.tsv" "$WORK/dnds"
if [ "${FEASIBILITY:-0}" = 1 ]; then
  # blinded: cohort size, neutral expectation and stage mix only
  $PYTHON "$HERE/replicate.py" --work "$WORK" --out "$OUT" --feasibility
else
  $PYTHON "$HERE/replicate.py" --work "$WORK" --out "$OUT" ${CORNISH:+--cornish "$CORNISH"}
fi
