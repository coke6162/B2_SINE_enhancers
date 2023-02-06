#!/bin/bash

# Script for calling candidate regions for ABC model
# 'mm10.chrom.sizes' supplied in repos
# 'mm10-blocklist.v2.bed' can be obtained through link in README

# Example usage:
# sbatch call_candidate_regions.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=6:00:00
#SBATCH --mem=8G

# Job name and output
#SBATCH -J candidate_regions
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Purge modules
module purge

# Load modules
module load gcc/7.1.0 bedtools/2.28.0 python/3.7.4

# Specify names of bam and peak files
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Setting variable to bam and peak files...")

bamFile=$(ls *IFNG*ATAC_R1*.sorted.dedup.shifted.bam)
peakFiles=$(ls *IFNG*ATAC*_peaks.narrowPeak)

# Call candidate regions
echo $(date +"[%b %d %H:%M:%S] Calling candidate regions with makeCandidateRegions.py...")

python ABC-Enhancer-Gene-Prediction-0.2.2/src/makeCandidateRegions.py \
--narrowPeak plat_BMDM_IFNG_2h_ATAC_R1_peaks.narrowPeak \
--bam plat_BMDM_IFNG_2h_ATAC_R1.sorted.dedup.bam \
--chrom_sizes mm10.chrom.sizes \
--regions_blacklist mm10-blocklist.v2.bed \
--outDir ./ \
--peakExtendFromSummit 250 \
--nStrongestPeaks 150000

echo $(date +"[%b %d %H:%M:%S] Done")