#!/bin/bash

## Script for running deeptools to shift fragments +4 and -5 for positive and negative strand, respectively

# Example usage:
# sbatch --array 0-8 shift_fragments

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=4:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J shift_fragments
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Load modules
module load samtools/1.14

# Define query files
queries=($(ls *.sorted.dedup.bam | xargs -n 1 basename))

# Define variables
tmpDir=tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
numCPU=8

# Make temporary directories
pwd; hostname; date

echo "Making temporary directories..."

mkdir -p ${tmpDir}/samtools
mkdir -p ${tmpDir}/deeptools

# Define temporary directory for deeptools
TMPDIR="${tmpDir}/deeptools"
export TMPDIR

# Run "alignmentSieve" to shift reads
echo $(date +"[%b %d %H:%M:%S] Starting alignmentSieve...")

echo "Processing sample: "${queries[$SLURM_ARRAY_TASK_ID]}

alignmentSieve \
--bam ${queries[$SLURM_ARRAY_TASK_ID]} \
--outFile ${queries[$SLURM_ARRAY_TASK_ID]%.sorted.dedup.bam}.sorted.dedup.shifted.bam.tmp \
--ATACshift \
--numberOfProcessors ${numCPU}

# Sort shifted bams
echo $(date +"[%b %d %H:%M:%S] Sorting shifted bams...")

samtools sort \
-@ ${numCPU} \
${queries[$SLURM_ARRAY_TASK_ID]%.sorted.dedup.bam}.sorted.dedup.shifted.bam.tmp \
> ${queries[$SLURM_ARRAY_TASK_ID]%.sorted.dedup.bam}.sorted.dedup.shifted.bam

# Index sorted bam files
echo $(date +"[%b %d %H:%M:%S] Indexing sorted shifted bams...")

samtools index \
${queries[$SLURM_ARRAY_TASK_ID]%.sorted.dedup.bam}.sorted.dedup.shifted.bam \
> ${queries[$SLURM_ARRAY_TASK_ID]%.sorted.dedup.bam}.sorted.dedup.shifted.bam.bai

# Remove intermediate bam files
echo $(date +"[%b %d %H:%M:%S] Removing intermediate files...")

rm ${queries[$SLURM_ARRAY_TASK_ID]%.sorted.dedup.bam}.sorted.dedup.shifted.bam.tmp

echo $(date +"[%b %d %H:%M:%S] Done")