#!/bin/bash

# Parse sam file for ligation pairs & write pairs file with pairtools parse
# 'mm10.chrom.sizes' supplied through github repos

# Example usage:
# sbatch --array 0-5 pairtools_parse.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=24:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J parse
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Set constant variables
numThreads=2

# Load module
module load pairtools/1.0.2

# Define query files
queries=($(ls *.sam | xargs -n 1 basename | sed 's/.sam//g' | uniq))

# Convert bam to pairsam
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo "Converting bam to pairsam..."

pairtools parse \
--chroms-path mm10.chrom.sizes \
--assembly mm10 \
--add-columns mapq \
--output ${queries[$SLURM_ARRAY_TASK_ID]}.pairs.gz \
${queries[$SLURM_ARRAY_TASK_ID]}.sam

echo $(date +"[%b %d %H:%M:%S] Done")