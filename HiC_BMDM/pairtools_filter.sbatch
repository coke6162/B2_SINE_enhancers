#!/bin/bash

# Script for filtering pair files using pairtools
# 'mm10.chrom.sizes' supplied through github repos

# Example usage:
# sbatch --array 0-2 pairtools_filter.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=24:00:00
#SBATCH --mem=4GB

# Job name and output
#SBATCH -J filter
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Set constant variables
numThreads=4

# Define query files
queries=($(ls *.marked.merged.sorted.pairs.gz | xargs -n 1 basename | sed 's/.marked.merged.sorted.pairs.gz//g' | uniq))

# Retain only unique alignments with classification UU and UC
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Selecting alignments...")

pairtools select \
'(pair_type=="UU") or (pair_type=="UC")' \
--output ${queries[$SLURM_ARRAY_TASK_ID]}.tmp1.gz \
${queries[$SLURM_ARRAY_TASK_ID]}.marked.merged.sorted.pairs.gz

# Run pairtools split
echo $(date +"[%b %d %H:%M:%S] Running pairtools split...")

pairtools split \
--output-pairs ${queries[$SLURM_ARRAY_TASK_ID]}.tmp2.gz \
${queries[$SLURM_ARRAY_TASK_ID]}.tmp1.gz

# Run pairtools select
echo $(date +"[%b %d %H:%M:%S] Running pairtools select...")

pairtools select \
'True' \
--chrom-subset mm10.chrom.sizes \
--output ${queries[$SLURM_ARRAY_TASK_ID]}.filtered.pairs.gz \
${queries[$SLURM_ARRAY_TASK_ID]}.tmp2.gz

# Index filtered file
echo $(date +"[%b %d %H:%M:%S] Indexing filtered file with pairix (sanity check)...")

pairix ${queries[$SLURM_ARRAY_TASK_ID]}.filtered.pairs.gz

echo $(date +"[%b %d %H:%M:%S] Done")