#!/bin/bash
# Driver networks in TCGA COADREAD (see README.md)
set -euo pipefail
cd "$(dirname "$0")"
PYTHON=${PYTHON:-python}
mkdir -p data results

Rscript mutations.R
$PYTHON classify.py
Rscript dnds.R
$PYTHON gene_sets.py
$PYTHON crc_pathways.py ALL MSS MSI
$PYTHON crc_networks.py
$PYTHON crc_robust.py
$PYTHON crc_immune.py
$PYTHON crc_power.py
$PYTHON crc_summary.py
