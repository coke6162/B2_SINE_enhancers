#!/bin/bash

# Script for subsetting enhancers by ABC score and intersecting with H3K27ac peaks

# Example usage:
# sbatch subset_and_intersect_enhancers.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=6:00:00
#SBATCH --mem=8G

# Job name and output
#SBATCH -J subset_and_intersect
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Load modules
module load bedtools/2.28.0

# Decompress file
pwd; hostname; date

echo "Decompressing enhancer txt file..."

gzip -d EnhancerPredictionsAllPutative.txt.gz

# Convert predicted enhancer list to bed
echo $(date +"[%b %d %H:%M:%S] Converting predicted enhancers to bed format...")

awk '(NR>1){print $1 "\t" $2 "\t" $3 "\t" $7 "\t" $21 "\t" "."}' EnhancerPredictionsAllPutative.txt \
| bedtools sort -i - \
> EnhancerPredictionsAllPutative.bed

# Merge overlapping peaks & assign each peak an arbitrary name
echo $(date +"[%b %d %H:%M:%S] Assigning each peak with an arbitrary name...")

bedtools merge -i picc_BMDM_WT_UT_H3K27ac_peaks.narrowPeak \
| awk '{print $1 "\t" $2 "\t" $3 "\t" "region"NR "\t" "."}' \
> ${enhancerDir}/picc_BMDM_WT_UT_H3K27ac_peaks_regions.bed

# Intersect complete enhancer list with IFNG-inducible peaks
echo $(date +"[%b %d %H:%M:%S] Intersecting enhancers with peaks...")

bedtools intersect -a EnhancerPredictionsAllPutative.bed -b picc_BMDM_WT_UT_H3K27ac_peaks_regions.bed -loj \
| grep "region" \
| awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "." "\t" $10}' \
| uniq \
> EnhancerPredictionsAllPutative_intersect_IFNG_H3K27ac.bed

# Filter for enhancer-gene pairs with ABC score > 0.01
echo $(date +"[%b %d %H:%M:%S] Intersecting subset enhancers with peaks...")

cat EnhancerPredictionsAllPutative_intersect_IFNG_H3K27ac.bed \
| awk '{if ($5!="NaN") print}' \
| awk '{if ($5>0.01) print}' \
| sort -nrk5 \
> EnhancerPredictionsAllPutative_intersect_IFNG_H3K27ac_scoreOver0.01.bed

echo $(date +"[%b %d %H:%M:%S] Done")