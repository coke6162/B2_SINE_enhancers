#!/bin/bash

## Script for running hisat2
## Example usage:
## inDir=path outDir=path hisatIdxDir=/Shares/CL_Shared/db/genomes/mm10/index/hisat2 hisatIdx=mm10.main rnaStrandness=RF sbatch --array 0-15 hisat2_PE.sbatch

# important note: 
# for single-ended data, rnaStrandness is R (for reverse) or F (for forward)
# for paired-ended data, rnaStrandness is RF (for reverse) or FR (for forward)

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=4:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J hisat2_PE
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Set constant variables
numThreads=8

# Set query files
# note: remove R1/R2 to generate a unique identifier for each pair of files
queries=($(ls $inDir/*_trimmed.fastq.gz | xargs -n 1 basename | sed 's/_R1_trimmed.fastq.gz//g' | sed 's/_R2_trimmed.fastq.gz//g' | uniq))

# Load modules
module load hisat2 samtools

# Run hisat2
pwd; hostname; date

echo "Starting HISAT2 alignment..."
echo "HISAT2 version: "$(hisat2 --version)
echo "samtools version: "$(samtools --version)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

hisat2 --no-softclip \
--rna-strandness ${rnaStrandness} \
-p ${numThreads} \
-x ${hisatIdxDir}/${hisatIdx} \
-1 $inDir/${queries[$SLURM_ARRAY_TASK_ID]}_R1_trimmed.fastq.gz \
-2 $inDir/${queries[$SLURM_ARRAY_TASK_ID]}_R2_trimmed.fastq.gz \
| samtools view -q 10 -Sb - \
| samtools sort -@ ${numThreads} - -o $outDir/${queries[$SLURM_ARRAY_TASK_ID]%.fastq.gz}.uniq.bam

# Create a bai index file for upload to UCSC browser
samtools index $outDir/${queries[$SLURM_ARRAY_TASK_ID]%.fastq.gz}.uniq.bam

echo $(date +"[%b %d %H:%M:%S] Done!")
