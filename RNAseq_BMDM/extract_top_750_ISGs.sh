#!/bin/bash

## Piccolo et al. 2017 IFNG 4h vs UT

# Get names for top 750 ISGs (by descending log2FC)
awk '(NR>1){if ($3 > 0 && $7 <= 0.05) print $0}' picc_BMDM_IFNG_4h_vs_UT.txt \
| sort -k3,3rn \
| awk '(NR<=750){print $1}' \
> picc_BMDM_IFNG_4h_vs_UT_ISGs_top750_names.txt

# Get coordindates for top 750 ISGs
awk 'FNR==NR {a[$5] = $0; next} $1 in a {print a[$1]}' picc_BMDM_IFNG_4h_vs_UT.bed picc_BMDM_IFNG_4h_vs_UT_ISGs_top750_names.txt \
| bedtools sort -i - \
> picc_BMDM_IFNG_4h_vs_UT_ISGs_top750.bed

## Piccolo et al. 2017 IFNG 2h vs UT

# Get names for top 750 ISGs (by descending log2FC)
awk '(NR>1){if ($3 > 0 && $7 <= 0.05) print $0}' picc_BMDM_IFNG_2h_vs_UT.txt \
| sort -k3,3rn \
| awk '(NR<=750){print $1}' \
> picc_BMDM_IFNG_2h_vs_UT_ISGs_top750_names.txt

# Get coordindates for top 750 ISGs
awk 'FNR==NR {a[$5] = $0; next} $1 in a {print a[$1]}' picc_BMDM_IFNG_2h_vs_UT.bed picc_BMDM_IFNG_2h_vs_UT_ISGs_top750_names.txt \
| bedtools sort -i - \
> picc_BMDM_IFNG_2h_vs_UT_ISGs_top750.bed

## Platanitis et al. 2019 IFNG 2h vs UT

# Get names for top 750 ISGs (by descending log2FC)
awk '(NR>1){if ($3 > 0 && $7 <= 0.05) print $0}' plat_BMDM_IFNG_2h_vs_UT.txt \
| sort -k3,3rn \
| awk '(NR<=750){print $1}' \
> plat_BMDM_IFNG_2h_vs_UT_ISGs_top750_names.txt

# Get coordindates for top 750 ISGs
awk 'FNR==NR {a[$5] = $0; next} $1 in a {print a[$1]}' plat_BMDM_IFNG_2h_vs_UT.bed plat_BMDM_IFNG_2h_vs_UT_ISGs_top750_names.txt \
| bedtools sort -i - \
> plat_BMDM_IFNG_2h_vs_UT_ISGs_top750.bed