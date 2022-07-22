#!/bin/bash

# Script for calling H3K27ac and STAT1 peaks using data from Piccolo et al. 2017
# H3K27ac peaks are called relative to provided input using q-value cutoff of 0.01
# STAT1 peaks are called relative to provided input using p-value cutoff of 0.01
# Relaxed cutoff used because only peaks intersecting both replicates will be retained (see intersect_peak_replicates.sbatch)

# Example usage:
# sbatch macs2_piccolo.sbatch

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

# Call peaks for UT H3K27ac
echo $(date +"[%b %d %H:%M:%S] Calling UT H3K27ac peaks...")

macs2 callpeak \
-t picc_BMDM_WT_UT_H3K27ac.sorted.dedup.bam \
-c picc_BMDM_WT_UT_input.sorted.dedup.bam \
--name picc_BMDM_WT_UT_H3K27ac \
--format BAM --gsize mm --qvalue 0.01 \
--bdg --SPMR --call-summits

# Call peaks for IFNG 2h H3K27ac
echo $(date +"[%b %d %H:%M:%S] Calling IFNG 2h H3K27ac peaks...")

macs2 callpeak \
-t picc_BMDM_WT_IFNG_2h_H3K27ac.sorted.dedup.bam \
-c picc_BMDM_WT_UT_input.sorted.dedup.bam \
--name picc_BMDM_WT_IFNG_2h_H3K27ac \
--format BAM --gsize mm --qvalue 0.01 \
--bdg --SPMR --call-summits

# Call peaks for UT STAT1
echo $(date +"[%b %d %H:%M:%S] Calling UT STAT1 peaks...")

macs2 callpeak \
-t picc_BMDM_WT_UT_STAT1.sorted.dedup.bam \
-c picc_BMDM_WT_UT_input.sorted.dedup.bam \
--name picc_BMDM_WT_UT_STAT1 \
--format BAM --gsize mm --pvalue 0.01 \
--bdg --SPMR --call-summits

# Call peaks for IFNG 2h STAT1 R1
echo $(date +"[%b %d %H:%M:%S] Calling IFNG 2h STAT1 R1 peaks...")

macs2 callpeak \
-t picc_BMDM_WT_IFNG_2h_STAT1_R1.sorted.dedup.bam \
-c picc_BMDM_WT_UT_input.sorted.dedup.bam \
--name picc_BMDM_WT_IFNG_2h_STAT1_R1 \
--format BAM --gsize mm --pvalue 0.01 \
--bdg --SPMR --call-summits

# Call peaks for IFNG 2h STAT1 R2
echo $(date +"[%b %d %H:%M:%S] Calling IFNG 2h STAT1 R2 peaks...")

macs2 callpeak \
-t picc_BMDM_WT_IFNG_2h_STAT1_R2.sorted.dedup.bam \
-c picc_BMDM_WT_UT_input.sorted.dedup.bam \
--name picc_BMDM_WT_IFNG_2h_STAT1_R2 \
--format BAM --gsize mm --pvalue 0.01 \
--bdg --SPMR --call-summits

echo $(date +"[%b %d %H:%M:%S] Done")