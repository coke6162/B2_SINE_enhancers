#!/bin/bash

# Script for running salmon quant

# Example usage:
# sbatch --array 0-29 salmon_quant.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=4:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J salmon_quant
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Set constant variables
numThreads=4

# Load salmon
module load salmon/1.9.0

# Set query files
# note: remove R1/R2 to generate a unique identifier for each pair of files
queries=($(ls *.fastq.gz | xargs -n 1 basename | sed 's/_R1_trimmed.fastq.gz//g' | sed 's/_R2_trimmed.fastq.gz//g' | uniq))

# Run salmon
pwd; hostname; date

echo "Processing sample: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Starting salmon pseudo-alignment...")

salmon quant \
-i index/salmon \
-p ${numThreads} \
--geneMap all_samples_merged_transcripts_decoy.fa \
-1 <(gunzip -c ${queries[$SLURM_ARRAY_TASK_ID]}_R1_trimmed.fastq.gz) \
-2 <(gunzip -c ${queries[$SLURM_ARRAY_TASK_ID]}_R2_trimmed.fastq.gz) \
-o ${queries[$SLURM_ARRAY_TASK_ID]} \
--validateMappings --rangeFactorizationBins 4 --gcBias --libType A

echo $(date +"[%b %d %H:%M:%S] Done")