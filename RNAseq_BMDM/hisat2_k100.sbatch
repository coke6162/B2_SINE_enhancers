#!/bin/bash

# Script for aligning trimmed reads to the reference genome using hisat2
# hisat2 index located in a subdirectory named "index" (not supplied)
# all hisat2 index files share the same "mm10" prefix
# IMPORTANT: '-k 100' flag should be included for TEtranscripts analysis

# Example usage:
# sbatch --array 0-14 hisat2_k100.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=4:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J hisat2_k100
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Set constant variables
numThreads=4

# Load hisat2 and samtools
module load hisat2/2.1.0 samtools/1.10

# Define query files
# Removes suffix to generate a unique identifier for each pair of files
queries=($(ls *_trimmed.fastq.gz | sed 's/_trimmed.fastq.gz//g' | uniq))

# Align trimmed reads
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Aligning trimmed fastq with hisat2...")

hisat2 \
-k 100 -p ${numThreads} --no-softclip \
-x index/mm10 \
-U ${queries[$SLURM_ARRAY_TASK_ID]}_trimmed.fastq.gz \
| samtools view -Sb - \
| samtools sort -@ ${numThreads} - -o ${queries[$SLURM_ARRAY_TASK_ID]}_k100.sorted.bam

# Index filtered bam
samtools index ${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam

echo $(date +"[%b %d %H:%M:%S] Done")