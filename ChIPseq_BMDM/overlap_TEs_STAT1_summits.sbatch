#!/bin/bash

# Script for identifying TEs that overlap STAT1 summits
# "mm10_dfam.bed" is available in gzip'd format in the "RNAseq_BMDM" repos
# To generate "GAS_fimo_out.bed" and "ISRE_fimo_out.bed" see "motif_analysis" repos

# Example usage:
# sbatch overlap_TEs_STAT1_summits.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J overlapping_TEs
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Load bedtools
module load bedtools/2.28.0

# Extract B2 coordinates
grep "B2_Mm2" mm10_dfam.bed > B2_Mm2.bed

# Convert intersected Piccolo et al. IFNG 2h STAT1 narrowPeaks to summits

# Get peak names
cat picc_BMDM_WT_IFNG_2h_STAT1_R1_peaks.narrowPeak picc_BMDM_WT_IFNG_2h_STAT1_R2_peaks.narrowPeak \
| bedtools sort -i - \
| bedtools merge -i - -c 5,4 -o max,collapse \
| awk '{if ($5 ~ "_R1" && $5 ~ "_R2") print $1 "\t" $2 "\t" $3 "\t" "." "\t" $4 "\t" $5}' \
| tr "," "\t" \
| cut -f6 \
> picc_BMDM_WT_IFNG_2h_STAT1_intersected_peak_names_first_occurrence.txt

# Concatenate summit files
cat picc_BMDM_WT_IFNG_2h_STAT1_R1_summits.bed picc_BMDM_WT_IFNG_2h_STAT1_R2_summits.bed \
| bedtools sort -i - \
> picc_BMDM_WT_IFNG_2h_STAT1_concatenated_summits.bed

# Extract summits by matching peak name
awk 'FNR==NR {a[$4]=$0; next} $1 in a {print a[$0]}' picc_BMDM_WT_IFNG_2h_STAT1_concatenated_summits.bed picc_BMDM_WT_IFNG_2h_STAT1_intersected_peak_names_first_occurrence.txt \
> picc_BMDM_WT_IFNG_2h_STAT1_intersected_summits_first_occurrence.bed

# Identify TEs that overlap STAT1 summits
bedtools intersect -a picc_BMDM_WT_IFNG_2h_STAT1_intersected_summits_first_occurrence.bed -b mm10_dfam_filtered.bed -wa -wb \
cut -f6- \
bedtools sort -i - \
> picc_BMDM_WT_IFNG_2h_STAT1_intersected_summits_first_occurrence_repeats.bed

# Intersect with predicted GAS and ISRE motifs
cat GAS_fimo_out.bed ISRE_fimo_out.bed \
| bedtools sort -i - \
| bedtools intersect -a picc_BMDM_WT_IFNG_2h_STAT1_intersected_summits_first_occurrence_repeats.bed -b - -wa -wb \
> picc_BMDM_WT_IFNG_2h_STAT1_intersected_summits_first_occurrence_repeats_GAS_ISRE.bed