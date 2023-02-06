#!/bin/bash

# Script for marking duplicates with pairtools dedup

# Example usage:
# sbatch --array 0-2 pairtools_dedup.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=24:00:00
#SBATCH --mem=16GB

# Job name and output
#SBATCH -J markasdup
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Set constant variables
numThreads=4

# Define query files
queries=($(ls *.merged.sorted.pairs.gz | xargs -n 1 basename | sed 's/.merged.sorted.pairs.gz//g' | uniq))

# Mark duplicates
pwd; hostname; date

echo "Marking duplicates in pairsam file..."
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

pairtools dedup \
--mark-dups \
--output-dups ${queries[$SLURM_ARRAY_TASK_ID]}.dups.merged.sorted.pairs.gz \
--output-unmapped ${queries[$SLURM_ARRAY_TASK_ID]}.unmapped.merged.sorted.pairs.gz \
--output ${queries[$SLURM_ARRAY_TASK_ID]}.marked.merged.sorted.pairs.gz \
${queries[$SLURM_ARRAY_TASK_ID]}.merged.sorted.pairs.gz

# Index marked file
echo $(date +"[%b %d %H:%M:%S] Indexing marked file with pairix (sanity check)...")

pairix ${queries[$SLURM_ARRAY_TASK_ID]}.marked.merged.sorted.pairs.gz

echo $(date +"[%b %d %H:%M:%S] Done")