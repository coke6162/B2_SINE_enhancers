#!/bin/bash

# Get names for top 750 IRGs (by ascending log2FC)
awk '(NR>1){if ($3 < 0 && $7 <= 0.05) print $0}' picc_BMDM_IFNG_4h_vs_UT.txt \
| sort -k3,3rn \
| awk '(NR<=750){print $1}' \
> picc_BMDM_IFNG_4h_vs_UT_IRGs_top750_names.txt

# Get coordindates for top 750 IRGs
awk 'FNR==NR {a[$5] = $0; next} $1 in a {print a[$1]}' picc_BMDM_IFNG_4h_vs_UT.bed picc_BMDM_IFNG_4h_vs_UT_IRGs_top750_names.txt \
| bedtools sort -i - \
> picc_BMDM_IFNG_4h_vs_UT_IRGs_top750_names.bed