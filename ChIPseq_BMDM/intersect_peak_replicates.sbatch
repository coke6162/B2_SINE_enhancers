#!/bin/bash

# Script for getting peaks that intersect across replicates

# Example usage:
# sbatch intersect_peak_replicates.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J intersect
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Load bedtools
module load bedtools/2.28.0

# Not enough Piccolo et al. UT H3K27ac replicates

# Not enough Piccolo et al. IFNG 2h H3K27ac replicates

# Not enough Piccolo et al. UT STAT1 replicates

# Piccolo et al. IFNG 2h STAT1
cat picc_BMDM_WT_IFNG_2h_STAT1_R1_peaks.narrowPeak picc_BMDM_WT_IFNG_2h_STAT1_R2_peaks.narrowPeak \
| bedtools sort -i - \
| bedtools merge -i - -c 5,4 -o max,collapse \
| awk '{if ($5 ~ "_R1" && $5 ~ "_R2") print $1 "\t" $2 "\t" $3 "\t" "." "\t" $4}' \
> picc_BMDM_WT_IFNG_2h_STAT1_intersected.bed

# Platanitis et al. UT STAT1
cat plat_BMDM_WT_UT_STAT1_R1_peaks.narrowPeak plat_BMDM_WT_UT_STAT1_R2_peaks.narrowPeak \
| bedtools sort -i - \
| bedtools merge -i - -c 5,4 -o max,collapse \
| awk '{if ($5 ~ "_R1" && $5 ~ "_R2") print $1 "\t" $2 "\t" $3 "\t" "." "\t" $4}' \
> plat_BMDM_WT_UT_STAT1_intersected.bed

# Platanitis et al. IFNG 1.5h STAT1
cat plat_BMDM_WT_IFNG_1.5h_STAT1_R1_peaks.narrowPeak plat_BMDM_WT_IFNG_1.5h_STAT1_R2_peaks.narrowPeak \
| bedtools sort -i - \
| bedtools merge -i - -c 5,4 -o max,collapse \
| awk '{if ($5 ~ "_R1" && $5 ~ "_R2") print $1 "\t" $2 "\t" $3 "\t" "." "\t" $4}' \
> plat_BMDM_WT_IFNG_1.5h_STAT1_intersected.bed

# Cuartero et al. 2018 UT RAD21
cat cuar_BMDM_WT_UT_RAD21_R1_peaks.narrowPeak cuar_BMDM_WT_UT_RAD21_R2_peaks.narrowPeak \
| bedtools sort -i - \
| bedtools merge -i - -c 5,4 -o max,collapse \
| awk '{if ($5 ~ "_R1" && $5 ~ "_R2") print $1 "\t" $2 "\t" $3 "\t" "." "\t" $4}' \
> cuar_BMDM_WT_UT_RAD21_intersected.bed

# Gualdrini et al. 2022 UT CTCF
cat gual_BMDM_shNT_UT_CTCF_R1_peaks.narrowPeak gual_BMDM_shNT_UT_CTCF_R2_peaks.narrowPeak \
| bedtools sort -i - \
| bedtools merge -i - -c 5,4 -o max,collapse \
| awk '{if ($5 ~ "_R1" && $5 ~ "_R2") print $1 "\t" $2 "\t" $3 "\t" "." "\t" $4}' \
> gual_BMDM_shNT_UT_CTCF_intersected.bed