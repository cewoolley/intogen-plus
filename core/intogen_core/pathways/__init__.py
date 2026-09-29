"""
Pathway-level analyses: selection of gene sets through mutations, epigenetic
silencing and expression outliers, and co-occurrence of dysregulation events
in the same tumours.

Gene-level driver discovery misses genes that are rarely altered. When many
such genes belong to the same pathway, the pathway as a whole can show a
significant signal. Each gene set is tested using all its genes and using
only its "long tail": the genes that are not individually significant.
"""
