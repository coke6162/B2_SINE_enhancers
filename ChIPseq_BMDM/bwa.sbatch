#!/bin/bash

# Script for aligning trimmed, single-end reads to the reference genome using bwa mem
# bwa index located in a subdirectory named "index" (not supplied)
# all bwa index files share the same "mm10" prefix

# Example usage:
# sbatch --array 0-11 bwa.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=4:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J bwa
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Set constant variables
numThreads=4

# Load modules
module load bwa/0.7.15 samtools/1.8

# Define query files
queries=($(ls *_trimmed.fastq.gz | xargs -n 1 basename | sed 's/_trimmed.fastq.gz//'))

# Define tmpDir var
tmpDir=tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Make temporary directory
pwd; hostname; date

echo "Making temporary directories..."

mkdir -p ${tmpDir}

# Align trimmed reads to the reference genome
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Aligning to assembly...")

bwa mem \
-t ${numThreads} \
index/mm10 \
${queries[$SLURM_ARRAY_TASK_ID]}_trimmed.fastq.gz \
| samtools view -@ ${numThreads} -Sb -q 10 -F 4 \
| samtools sort -@ ${numThreads} -T ${tmpDir} \
> ${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam

# Index sorted bam files
echo $(date +"[%b %d %H:%M:%S] Indexing sorted bam...")

samtools index \
${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam \
${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam.bai

# Remove temporary directory
echo $(date +"[%b %d %H:%M:%S] Removing temporary directory...")

rm -r ${tmpDir}/

echo $(date +"[%b %d %H:%M:%S] Done")