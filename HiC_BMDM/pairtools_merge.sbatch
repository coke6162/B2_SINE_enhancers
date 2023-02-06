#!/bin/bash

# Script for merging sorted pairsam files with pairtools merge
# 'sample_names.txt' supplied through github repos

# Example usage:
# sbatch --array 0-2 pairtools_merge.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=12:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Set constant variables
numThreads=8

# Define query
queries=($(awk '(NR>1){print $2}' sample_names.txt | sed 's/_R[0-9]//g' | sort | uniq))

# Define bam files
bams=$(ls ${queries[$SLURM_ARRAY_TASK_ID]}*.sorted.pairs.gz)

# Define temporary directory
tmpDir=pairtools_tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Make temporary directory
pwd; hostname; date

echo "Making temporary directories..."

mkdir -p ${tmpDir}

# Sort pairsam file
echo $(date +"[%b %d %H:%M:%S] Sorting pairsam file...")

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

pairtools merge \
--nproc ${numThreads} \
--memory 32G \
--tmpdir ${tmpDir} \
--output ${queries[$SLURM_ARRAY_TASK_ID]}.merged.sorted.pairs.gz \
${bams}

# Remove temporary directory
echo $(date +"[%b %d %H:%M:%S] Removing temporary directory...")

rm -r ${tmpDir}

echo $(date +"[%b %d %H:%M:%S] Done")