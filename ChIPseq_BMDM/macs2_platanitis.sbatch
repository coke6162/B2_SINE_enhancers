#!/bin/bash

# Script for calling STAT1 peaks using data from Platanitis et al. 2019
# STAT1 peaks are called relative to provided input using p-value cutoff of 0.01

# Example usage:
# sbatch macs2_platanitis.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=12:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J macs2
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Load MACS2
module load python/2.7.14/MACS/2.1.1

# Call peaks for UT STAT1 R1
echo $(date +"[%b %d %H:%M:%S] Calling UT STAT1 R1 peaks...")

macs2 callpeak \
-t plat_BMDM_WT_UT_STAT1_R1.sorted.dedup.bam \
-c plat_BMDM_input.sorted.dedup.bam \
--name plat_BMDM_WT_UT_STAT1_R1 \
--format BAM --gsize mm --pvalue 0.01 \
--bdg --SPMR --call-summits

# Call peaks for UT STAT1 R2
echo $(date +"[%b %d %H:%M:%S] Calling UT STAT1 R2 peaks...")

macs2 callpeak \
-t plat_BMDM_WT_UT_STAT1_R2.sorted.dedup.bam \
-c plat_BMDM_input.sorted.dedup.bam \
--name plat_BMDM_WT_UT_STAT1_R2 \
--format BAM --gsize mm --pvalue 0.01 \
--bdg --SPMR --call-summits

# Call peaks for IFNG 1.5h STAT1 R1
echo $(date +"[%b %d %H:%M:%S] Calling IFNG 1.5h STAT1 R1 peaks...")

macs2 callpeak \
-t plat_BMDM_WT_IFNG_1.5h_STAT1_R1.sorted.dedup.bam \
-c plat_BMDM_input.sorted.dedup.bam \
--name plat_BMDM_WT_IFNG_1.5h_STAT1_R1 \
--format BAM --gsize mm --pvalue 0.01 \
--bdg --SPMR --call-summits

# Call peaks for IFNG 1.5h STAT1 R2
echo $(date +"[%b %d %H:%M:%S] Calling IFNG 1.5h STAT1 R2 peaks...")

macs2 callpeak \
-t plat_BMDM_WT_IFNG_1.5h_STAT1_R2.sorted.dedup.bam \
-c plat_BMDM_input.sorted.dedup.bam \
--name plat_BMDM_WT_IFNG_1.5h_STAT1_R2 \
--format BAM --gsize mm --pvalue 0.01 \
--bdg --SPMR --call-summits

echo $(date +"[%b %d %H:%M:%S] Done")