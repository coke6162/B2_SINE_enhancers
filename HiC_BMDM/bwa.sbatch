#!/bin/bash

# Script for aligning paired-end reads to the reference genome using bwa mem
# bwa index located in a subdirectory named "index/bwa" (not supplied)
# all bwa index files share the same "mm10" prefix

# Example usage:
# sbatch --array 0-5 bwa.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=24:00:00
#SBATCH --mem=16GB

# Job name and output
#SBATCH -J bwa
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Set constant variables
numThreads=8

# Purge modules
module purge

# Load modules
module load samtools/1.14 bwa/0.7.17

# Define query files
queries=($(ls *.fastq.gz | xargs -n 1 basename | sed 's/_1.fastq.gz//g' | sed 's/_2.fastq.gz//g' | uniq))

# Define temporary directory
tmpDir=samtools_tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Make temporary directory
pwd; hostname; date

echo "Making temporary directories..."

mkdir -p ${tmpDir}

# Align reads to the reference genome
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Aligning to assembly...")

bwa mem -SP -t ${numThreads} \
bwa/index/mm10 \
${queries[$SLURM_ARRAY_TASK_ID]}_1.fastq.gz \
${queries[$SLURM_ARRAY_TASK_ID]}_2.fastq.gz \
> ${queries[$SLURM_ARRAY_TASK_ID]}.sam

# Generate flagstat files from unfiltered sam files
echo $(date +"[%b %d %H:%M:%S] Making flagstat file from unfiltered sam...")

samtools flagstat -@ ${numThreads} ${queries[$SLURM_ARRAY_TASK_ID]}.sam \
> ${queries[$SLURM_ARRAY_TASK_ID]}_flagstat.txt

# Convert sam files to bam files
echo $(date +"[%b %d %H:%M:%S] Converting sam to unfiltered bam & sorting")

samtools view -@ ${numThreads} -Sb ${queries[$SLURM_ARRAY_TASK_ID]}.sam \
| samtools sort -@ ${numThreads} -T ${tmpDir} - \
> ${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam

# Index unfiltered sorted bam files
echo $(date +"[%b %d %H:%M:%S] Indexing unfiltered, sorted bams...")

samtools index \
${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam \
> ${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam.bai

# Remove intermediate files
echo $(date +"[%b %d %H:%M:%S] Removing intermediate files...")

rm -r ${tmpDir}

echo $(date +"[%b %d %H:%M:%S] Done")