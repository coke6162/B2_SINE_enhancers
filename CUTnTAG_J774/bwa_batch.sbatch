#!/bin/bash

## Script for aligning to bwa

## Example usage:
## Note: array should equal the number of .bam files want at end
## inDir=path outDir=path sbatch --array 0-21 bwa_batch.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=6:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J bwa_batch.sbatch
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# load modules
module load bwa
module load samtools

# define genome index
genomeIndex="/Shares/CL_Shared/db/genomes/mm10/index/bwa/mm10.primary"

# define parameters
numCPUs=60
numThreads=4

# define queries
# breaks files down into base identifier ("sed" command strips R1 or R2 so can add it later to pair files together)
queries=($(ls ${inDir}/*fastq.gz | xargs -n 1 basename | sed 's/_R1_trimmed.fastq.gz//g' | sed 's/_R2_trimmed.fastq.gz//g' | uniq))

# print time and date
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Starting bwa alignment...")

# Run bwa and save .bam file in outDir
bwa mem \
-t ${numCPUs} \
${genomeIndex} \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R1_trimmed.fastq.gz \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R2_trimmed.fastq.gz \
| samtools view -Sb -q 10 \
| samtools sort -@ ${numThreads} - \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.fastq.gz}.sorted.bam

echo $(date +"[%b %d %H:%M:%S] Done!")
