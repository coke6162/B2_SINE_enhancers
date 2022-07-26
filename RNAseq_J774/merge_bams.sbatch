#!/bin/bash

## Script for running picard MergeSamFiles
## Example usage:
## inDir=path \
## outDir=path \
## sbatch --array 0-2 mergeBams.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=24:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J mergeBams
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Set constant variables
numThreads=8

# Define path variables
tmpDir=${inDir}/samtools_tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Define query samples
# note: samples.txt is the text file where you define which .bam files are to be merged
queries=($(cat ${inDir}/samples.txt))

# Load modules
module load samtools

# Make temporary directory
pwd; hostname; date

echo "samtools version: "$(samtools --version)

echo "Processing sample: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}

echo "Making temporary directory..."

mkdir -p ${tmpDir}

# Create a text file specifying all bam files to merge
# This works if all your samples share a common suffix & end in ".bam"
echo $(date +"[%b %d %H:%M:%S] Creating a txt file giving names of all files to merge...")

ls ${inDir} | grep "^${queries[$SLURM_ARRAY_TASK_ID]}.*bam$" - | sed 's|^|'${inDir}'/|' > ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_merge_list.txt

# Merge bams with samtools
echo $(date +"[%b %d %H:%M:%S] Merging bam files with samtools...")

samtools merge \
-@ ${numThreads} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_merged.bam \
$(cat ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_merge_list.txt)

# Sort merged bam
echo $(date +"[%b %d %H:%M:%S] Sorting merged bam...")

samtools sort \
-T ${tmpDir} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_merged.bam \
-o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_merged.uniq.bam

# Index merged bam
echo $(date +"[%b %d %H:%M:%S] Indexing merged bam...")

samtools index \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_merged.uniq.bam \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_merged.uniq.bam.bai

echo $(date +"[%b %d %H:%M:%S] Done")
