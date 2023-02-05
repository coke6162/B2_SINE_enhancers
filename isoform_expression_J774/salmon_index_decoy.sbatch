#!/bin/bash

# Script for generating a decoy-aware index using the merged GTF (from Stringtie)
# Assumes mm10 fasta file "mm10.fa" is in current directory

# Example usage:
# sbatch salmon_index_decoy.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=2:00:00
#SBATCH --mem=64GB

# Job name and output
#SBATCH -J salmon_index
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Set constant variables
numThreads=8
kmerLen=31

# Load modules
module load cufflinks/2.2.1 samtools/1.14 salmon/1.9.0

# Convert merged GTF to fasta
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Converting merged GTF to fasta...")

gffread all_samples_merged.gtf -g mm10.fa -w all_samples_merged_transcripts.fa

# Append mm10 FASTA to transcriptome FASTA
echo $(date +"[%b %d %H:%M:%S] Appending mm10 fasta to transcriptome fasta...")

cat all_samples_merged_transcripts.fa mm10.fa > all_samples_merged_transcripts_decoy.fa

# Index mm10 fasta
echo $(date +"[%b %d %H:%M:%S] Indexing mm10 fasta...")

samtools faidx mm10.fa

# Generate list of mm10 chromosome names
echo $(date +"[%b %d %H:%M:%S] Generating list of mm10 chromosome names...")

awk '{print $1}' mm10.fa.fai > decoys.txt

# Run salmon index
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Generating decoy-aware index...")

salmon index \
-p ${numThreads} \
-k ${kmerLen} \
--gencode \
-t all_samples_merged_transcripts_decoy.fa \
-d decoys.txt \
-i index/

echo $(date +"[%b %d %H:%M:%S] Done")