#!/bin/bash

# Script for identifying STAT1-bound B2_Mm2, shuffled B2_Mm2, and unbound B2_Mm2, B2_Mm1a, and B2_Mm1t
# "mm10_dfam.bed" is available in gzip'd format in the "RNAseq_BMDM" repos
# These files are meant to be used with "bedtools_closest.sbatch" and "nearest_neighbor_histograms.R"

# Example usage:
# sbatch get_overlapping_B2_nearest_neighbor.sbatch

## General settings
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
grep "B2_Mm1a" mm10_dfam.bed > B2_Mm1a.bed
grep "B2_Mm1t" mm10_dfam.bed > B2_Mm1t.bed

# Generate mm10.chrom.sizes
samtools faidx mm10.fa
awk '{print $1 "\t" $2}' mm10.fa.fai > mm10.chrom.sizes

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
mkdir -p nearest_neighbor/

# Identify IFNG STAT1-bound B2_Mm2
bedtools intersect -a B2_Mm2.bed -b picc_BMDM_WT_IFNG_2h_STAT1_intersected.bed -u > nearest_neighbor/B2_Mm2_picc_BMDM_WT_IFNG_2h_STAT1_intersected.bed
bedtools intersect -a B2_Mm2.bed -b plat_BMDM_WT_IFNG_1.5h_STAT1_intersected.bed -u > nearest_neighbor/B2_Mm2_plat_BMDM_WT_IFNG_1.5h_STAT1_intersected.bed

# Shuffle IFNG STAT1-bound B2_Mm2
for i in B2_Mm2_picc_BMDM_WT_IFNG_2h_STAT1_intersected.bed B2_Mm2_plat_BMDM_WT_IFNG_1.5h_STAT1_intersected.bed
do
	bedtools shuffle -i ${i} -g mm10.chrom.sizes -chrom -sed 1 -noOverlapping | bedtools sort -i - > nearest_neighbor/${i%.bed}_shuffled.bed
done

# Randomly subset for 2122 (Piccolo et al. IFNG 2h STAT1) and 838 (Platanitis et al. IFNG 1.5h STAT1) unbound B2_Mm2 elements
bedtools intersect -a B2_Mm2.bed -b picc_BMDM_WT_IFNG_2h_STAT1_intersected.bed -v | shuf -n 2122 --random-source=<(get_seeded_random ${seed}) ${REST} | bedtools sort -i - > nearest_neighbor/B2_Mm2_null_picc_BMDM_WT_IFNG_2h_STAT1_intersected_random2122.bed
bedtools intersect -a B2_Mm2.bed -b plat_BMDM_WT_IFNG_1.5h_STAT1_intersected.bed -v | shuf -n 838 --random-source=<(get_seeded_random ${seed}) ${REST} | bedtools sort -i - > nearest_neighbor/B2_Mm2_null_plat_BMDM_WT_IFNG_1.5h_STAT1_intersected_random838.bed

# Randomly subset for 2122 (Piccolo et al. IFNG 2h STAT1) and 838 (Platanitis et al. IFNG 1.5h STAT1) unbound B2_Mm1a elements
bedtools intersect -a B2_Mm1a.bed -b picc_BMDM_WT_IFNG_2h_STAT1_intersected.bed -v | shuf -n 2122 --random-source=<(get_seeded_random ${seed}) ${REST} | bedtools sort -i - > nearest_neighbor/B2_Mm1a_null_picc_BMDM_WT_IFNG_2h_STAT1_intersected_random2122.bed
bedtools intersect -a B2_Mm1a.bed -b plat_BMDM_WT_IFNG_1.5h_STAT1_intersected.bed -v | shuf -n 838 --random-source=<(get_seeded_random ${seed}) ${REST} | bedtools sort -i - > nearest_neighbor/B2_Mm1a_null_plat_BMDM_WT_IFNG_1.5h_STAT1_intersected_random838.bed

# Randomly subset for 2122 (Piccolo et al. IFNG 2h STAT1) and 838 (Platanitis et al. IFNG 1.5h STAT1) unbound B2_Mm1t elements
bedtools intersect -a B2_Mm1t.bed -b picc_BMDM_WT_IFNG_2h_STAT1_intersected.bed -v | shuf -n 2122 --random-source=<(get_seeded_random ${seed}) ${REST} | bedtools sort -i - > nearest_neighbor/B2_Mm1t_null_picc_BMDM_WT_IFNG_2h_STAT1_intersected_random2122.bed
bedtools intersect -a B2_Mm1t.bed -b plat_BMDM_WT_IFNG_1.5h_STAT1_intersected.bed -v | shuf -n 838 --random-source=<(get_seeded_random ${seed}) ${REST} | bedtools sort -i - > nearest_neighbor/B2_Mm1t_null_plat_BMDM_WT_IFNG_1.5h_STAT1_intersected_random838.bed