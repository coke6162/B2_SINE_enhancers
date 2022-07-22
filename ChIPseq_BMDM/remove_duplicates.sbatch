#!/bin/bash

# Script for removing duplicates using Picard MarkDuplicates

# Example usage:
# sbatch --array 0-11 remove_duplicates.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=4:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J remove_duplicates
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Set constant variables
numThreads=1

# Define path variables
tmpDir=java_tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Define query files
queries=($(ls *.sorted.bam | xargs -n 1 basename | sed 's/.sorted.bam//g' | sort | uniq))

# Make temporary directory
pwd; hostname; date

echo "Making temporary directory..."

mkdir -p ${tmpDir}

# Sort bam file
echo $(date +"[%b %d %H:%M:%S] Sorting bam files with picard...")

echo "Processing sample: "${queries[$SLURM_ARRAY_TASK_ID]}

java -Xmx4G -Djava.io.tmpdir=${tmpDir} -XX:ParallelGCThreads=${numThreads} \
-jar /opt/picard/2.6.0/picard-2.6.0.jar SortSam \
INPUT=${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam \
OUTPUT=${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam.tmp \
SORT_ORDER=coordinate \
VALIDATION_STRINGENCY=STRICT

# Remove duplicates
echo $(date +"[%b %d %H:%M:%S] Removing duplicates with picard...")

java -Xmx4G -Djava.io.tmpdir=${tmpDir} -XX:ParallelGCThreads=${numThreads} \
-jar /opt/picard/2.6.0/picard-2.6.0.jar MarkDuplicates \
INPUT=${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam.tmp \
OUTPUT=${queries[$SLURM_ARRAY_TASK_ID]}.sorted.dedup.bam.tmp \
METRICS_FILE=${queries[$SLURM_ARRAY_TASK_ID]}_dedup_metrics.txt \
VALIDATION_STRINGENCY=STRICT \
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
REMOVE_DUPLICATES=true

# Sort dedup'd bams
echo $(date +"[%b %d %H:%M:%S] Sorting dedup'd bams")

samtools sort \
-@ ${numThreads} \
${queries[$SLURM_ARRAY_TASK_ID]}.sorted.dedup.bam.tmp \
-o ${queries[$SLURM_ARRAY_TASK_ID]}.sorted.dedup.bam

# Index sorted dedup'd bams
echo $(date +"[%b %d %H:%M:%S] Indexing sorted dedup'd bams")

samtools index \
-@ {numThreads} \
${queries[$SLURM_ARRAY_TASK_ID]}.sorted.dedup.bam

# Remove temporary files and directories
echo $(date +"[%b %d %H:%M:%S] Removing temporary files & directories")

rm ${queries[$SLURM_ARRAY_TASK_ID]}*.bam.tmp
rm -r ${tmpDir}/

echo $(date +"[%b %d %H:%M:%S] Done")