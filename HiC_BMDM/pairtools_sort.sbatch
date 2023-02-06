#!/bin/bash

# Script for sorting pairsam with pairtools sort

# Example usage:
# sbatch --array 0-5 pairtools_sort.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=12:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J sort
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Set constant variables
numThreads=8

# Define query files
queries=($(ls *.pairs.gz | xargs -n 1 basename | sed 's/.pairs.gz//g' | uniq))

# Define temporary directory
tmpDir=pairtools_tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Make temporary directory
pwd; hostname; date

echo "Making temporary directories..."

mkdir -p ${tmpDir}

# Sort pairsam file
echo $(date +"[%b %d %H:%M:%S] Sorting pairsam file...")

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

pairtools sort \
--nproc ${numThreads} \
--memory 32G \
--tmpdir ${tmpDir} \
--output ${queries[$SLURM_ARRAY_TASK_ID]}.sorted.pairs.gz \
${queries[$SLURM_ARRAY_TASK_ID]}.pairs.gz

# Remove temporary directory
echo $(date +"[%b %d %H:%M:%S] Removing temporary directory...")

rm -r ${tmpDir}

echo $(date +"[%b %d %H:%M:%S] Done")