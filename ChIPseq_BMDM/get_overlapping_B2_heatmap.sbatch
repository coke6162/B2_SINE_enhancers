#!/bin/bash

# Script for identifying B2_Mm2 bound by STAT1 only, CTCF only, both STAT1 and CTCF, or neither
# For convenience all four bed files are available on the "ChIPseq_BMDM" repos
# "mm10_dfam.bed" is available in gzip'd format in the "RNAseq_BMDM" repos
# These files are meant to be used with "B2_Mm2_heatmap.sbatch"

# Example usage:
# sbatch get_overlapping_B2_heatmap.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J overlapping_B2
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Load bedtools
module load bedtools/2.28.0

# Extract B2 coordinates
grep "B2_Mm2" mm10_dfam.bed > B2_Mm2.bed

# Set shuffle seed
# Seeding adopted from https://stackoverflow.com/a/41962458/7820599
# Option parsing adopted from https://stackoverflow.com/a/14203146
get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}

seed=1;

REST=""
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
	-s)
	    seed="$2"
	    shift
	    shift
	    ;;
	*)    # unknown option
	    REST="$REST $1"
	    shift # past argument
	    ;;
    esac
done

# Make subdirectory
mkdir -p ChIP_heatmap/

# Identify IFNG STAT1-bound B2_Mm2
bedtools intersect -a B2_Mm2.bed -b picc_BMDM_WT_IFNG_2h_STAT1_intersected.bed -u > ChIP_heatmap/B2_Mm2_picc_BMDM_WT_IFNG_2h_STAT1_intersected.bed

# Identify UT CTCF-bound B2_Mm2
bedtools intersect -a B2_Mm2.bed -b gual_BMDM_shNT_UT_CTCF_intersected.bed -u > ChIP_heatmap/B2_Mm2_gual_BMDM_shNT_UT_CTCF_intersected.bed

# Randomly subset for 1000 control B2_Mm2 elements that do not overlap CTCF or STAT1
bedtools intersect -a B2_Mm2.bed -b ChIP_heatmap/B2_Mm2*.bed -v | awk '{if ($1 != "chrY") print}' | shuf -n 1000 --random-source=<(get_seeded_random ${seed}) ${REST} | bedtools sort -i - > ChIP_heatmap/B2_Mm2_control_no_STAT1_or_CTCF_shuf1000.bed

# Identify B2_Mm2 elements bound ONLY by STAT1 (i.e. without CTCF)
bedtools intersect -a ChIP_heatmap/B2_Mm2_picc_BMDM_WT_IFNG_2h_STAT1_intersected.bed -b ChIP_heatmap/B2_Mm2_gual_BMDM_shNT_UT_CTCF_intersected.bed -v > ChIP_heatmap/B2_Mm2_picc_BMDM_WT_IFNG_2h_STAT1_intersected_only.bed

# Identify B2_Mm2 elements bound ONLY by CTCF (i.e. without STAT1)
bedtools intersect -a ChIP_heatmap/B2_Mm2_gual_BMDM_shNT_UT_CTCF_intersected.bed -b ChIP_heatmap/B2_Mm2_picc_BMDM_WT_IFNG_2h_STAT1_intersected.bed -v > ChIP_heatmap/B2_Mm2_gual_BMDM_shNT_UT_CTCF_intersected_only.bed

# Identify B2_Mm2 elements bound by both STAT1 and CTCF
bedtools intersect -a ChIP_heatmap/B2_Mm2_picc_BMDM_WT_IFNG_2h_STAT1_intersected.bed -b ChIP_heatmap/B2_Mm2_gual_BMDM_shNT_UT_CTCF_intersected.bed -u > ChIP_heatmap/B2_Mm2_picc_BMDM_WT_IFNG_2h_STAT1_intersected_and_CTCF.bed