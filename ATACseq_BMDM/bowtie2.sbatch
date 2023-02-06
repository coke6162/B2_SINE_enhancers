#!/bin/bash

# Script for aligning trimmed, paired-end reads to the reference genome using bowtie2
# bowtie2 index located in a subdirectory named "index/bowtie2" (not supplied)
# all bowtie2 index files share the same "mm10" prefix

# sbatch --array=0-8 bowtie2.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=24:00:00
#SBATCH --mem=8G

# Job name and output
#SBATCH -J bowtie2
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Set constant variables
numThreads=8

# Module load bowtie2 and samtools
module load bowtie/2.2.9 samtools/1.14

# Define query files
queries=($(ls *fastq.gz | xargs -n 1 basename | sed 's/_trimmed_R1.fastq.gz//g' | sed 's/_trimmed_R2.fastq.gz//g' | uniq))

# Define tmpDir var
tmpDir=tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Make temporary directory
pwd; hostname; date

echo "Making temporary directories..."

mkdir -p ${tmpDir}

# Align reads to reference
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Aligning reads to the genome...")

bowtie2 \
--end-to-end --very-sensitive -X 1000 --fr --threads ${numThreads} \
-x index/bowtie2/mm10 \
-1 ${queries[$SLURM_ARRAY_TASK_ID]}_trimmed_R1.fastq.gz \
-2 ${queries[$SLURM_ARRAY_TASK_ID]}_trimmed_R2.fastq.gz \
| samtools view -@ ${numThreads} -Sb -q 10 -F 4 \
| samtools sort -@ ${numThreads} -T ${tmpDir} \
> ${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam

# Index sorted bam files
echo $(date +"[%b %d %H:%M:%S] Indexing sorted bam...")

samtools index \
${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam \
> ${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam.bai

# Determine percentage of total reads that mapped to mitochondrial genes
echo $(date +"[%b %d %H:%M:%S] Quantifying chrM proportion...")

chrMreads=`samtools view -c ${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam chrM`
totalReads=`samtools view -c ${queries[$SLURM_ARRAY_TASK_ID]}.sorted.bam`
fractionMreads=`echo "100 * ${chrMreads} / ${totalReads}" | bc -l`
touch ${queries[$SLURM_ARRAY_TASK_ID]}_chrMreadsFraction.txt
echo ${queries[$SLURM_ARRAY_TASK_ID]} >> ${queries[$SLURM_ARRAY_TASK_ID]}_chrMreadsFraction.txt
echo ${totalReads} 'total mapped reads' >> ${queries[$SLURM_ARRAY_TASK_ID]}_chrMreadsFraction.txt
echo ${chrMreads} 'mitochondrial reads' >> ${queries[$SLURM_ARRAY_TASK_ID]}_chrMreadsFraction.txt
echo ${fractionMreads} 'percentage of mitochondrial reads from total mapped reads' >> ${queries[$SLURM_ARRAY_TASK_ID]}_chrMreadsFraction.txt

# Remove tmpDir
echo $(date +"[%b %d %H:%M:%S] Removing tmpDir...")

rm -r ${tmpDir}

echo $(date +"[%b %d %H:%M:%S] Done")