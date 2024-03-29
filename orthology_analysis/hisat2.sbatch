#!/bin/bash

# Script for aligning trimmed reads to the reference genome using hisat2
# hisat2 index located in a subdirectory named "index/hisat2" (not supplied)
# all hisat2 index files share the same "hg38" prefix

# Example usage:
# sbatch --array 0-3 hisat2.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=4:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J hisat2
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Set constant variables
numThreads=4

# Load hisat2 and samtools
module load hisat2/2.1.0 samtools/1.14

# Define query files
# Removes suffix to generate a unique identifier for each pair of files
queries=($(ls *_trimmed.fastq.gz | sed 's/_trimmed.fastq.gz//g' | uniq))

# Align trimmed reads
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Aligning trimmed fastq with hisat2 and filtering poorly mapped/multimapped fragments...")

hisat2 \
-p ${numThreads} \
--no-softclip \
-x index/hisat2/hg38 \
-U /${queries[$SLURM_ARRAY_TASK_ID]}_trimmed.fastq.gz \
| samtools view -q 10 -Sb - \
| samtools sort -@ ${numThreads} - -o ${queries[$SLURM_ARRAY_TASK_ID]}.sorted.uniq.bam

# Index filtered bam
samtools index ${queries[$SLURM_ARRAY_TASK_ID]}.sorted.uniq.bam

echo $(date +"[%b %d %H:%M:%S] Done")