#!/bin/bash

## Script for running fastqc using an array
## Example usage:
## inDir=path \
## outDir=path \
## sbatch --array 0-21 fastqcreport.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=0:30:00
#SBATCH --mem=4GB

# Job name and output
#SBATCH -J fastqcreport
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Load FastQC
module load fastqc

# Define constant variables
numThreads=8

# Define query files
queries=($(ls ${inDir}/*.fastq.gz | xargs -n 1 basename))

# Run FastQC
pwd; hostname; date

echo "Running FastQC..."
echo "FastQC version: "$(fastqc --version)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

fastqc \
-o ${outDir} \
-f fastq \
-t ${numThreads} \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# '-f <file_type>' - input file type; options: fastq, bam, sam, bam_mapped, and sam_mapped
# '-t <int>' - number of threads
