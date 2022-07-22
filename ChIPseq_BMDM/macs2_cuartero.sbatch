#!/bin/bash

# Script for calling RAD21 peaks using data from Cuartero et al. 2018
# Peaks are called relative to provided input using p-value cutoff < 0.01
# Relaxed cutoff used because only peaks intersecting both replicates will be retained (see intersect_peak_replicates.sbatch)

# Example usage:
# sbatch macs2_cuartero.sbatch

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

# Call peaks for UT RAD21 R1
echo $(date +"[%b %d %H:%M:%S] Calling UT RAD21 R1 peaks...")

macs2 callpeak \
-t cuar_BMDM_WT_UT_RAD21_R1.sorted.dedup.bam \
-c cuar_BMDM_WT_UT_RAD21_input_R1.sorted.dedup.bam \
--name cuar_BMDM_WT_UT_RAD21_R1 \
--format BAM --gsize mm --pvalue 0.01 \
--bdg --SPMR --call-summits

# Call peaks for UT RAD21 R2
echo $(date +"[%b %d %H:%M:%S] Calling UT RAD21 R2 peaks...")

macs2 callpeak \
-t cuar_BMDM_WT_UT_RAD21_R2.sorted.dedup.bam \
-c cuar_BMDM_WT_UT_RAD21_input_R2.sorted.dedup.bam \
--name cuar_BMDM_WT_UT_RAD21_R2 \
--format BAM --gsize mm --pvalue 0.01 \
--bdg --SPMR --call-summits

echo $(date +"[%b %d %H:%M:%S] Done")