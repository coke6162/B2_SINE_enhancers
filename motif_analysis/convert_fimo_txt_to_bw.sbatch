#!/bin/bash

# Script for converting FIMO txt output to bigWig format

# Example usage:
# sbatch convert_fimo_txt_to_bw.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=1:00:00
#SBATCH --mem=2GB

# Job name and output
#SBATCH -J fimo_txt_to_bw
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Load modules
module purge
module load python/2.7.14/numpy/1.14.1
module load bedtools/2.28.0
module load samtools/1.10

# Convert FIMO txt output to bed
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Convert fimo output to bed format...")

awk '(NR>1){print $3 "\t" $4 "\t" $5 "\t" $1 "\t" $8 "\t" $6}' GAS_fimo_out.txt | bedtools sort -i - > GAS_fimo_out.bed
awk '(NR>1){print $3 "\t" $4 "\t" $5 "\t" $1 "\t" $8 "\t" $6}' ISRE_fimo_out.txt | bedtools sort -i - > ISRE_fimo_out.bed
awk '(NR>1){print $3 "\t" $4 "\t" $5 "\t" $1 "\t" $8 "\t" $6}' CTCF_fimo_out.txt | bedtools sort -i - > CTCF_fimo_out.bed

# Create a bedgraph from the bed output
echo $(date +"[%b %d %H:%M:%S] Create a bedgraph from bed output...")

bedtools sort -i GAS_fimo_out.bed | bedtools merge -i - -o max -c 5 > GAS_fimo_out.bdg
bedtools sort -i ISRE_fimo_out.bed | bedtools merge -i - -o max -c 5 > ISRE_fimo_out.bdg
bedtools sort -i CTCF_fimo_out.bed | bedtools merge -i - -o max -c 5 > CTCF_fimo_out.bdg

# Generate mm10.chrom.sizes
echo $(date +"[%b %d %H:%M:%S] Generating mm10.chrom.sizes...")

samtools faidx mm10
awk '{print $1 "\t" $2}' mm10.fa.fai > mm10.chrom.sizes

# Create a bigwig from the bedgraph output
echo $(date +"[%b %d %H:%M:%S] Create a bigwig from bedgraph output...")

bedGraphToBigWig GAS_fimo_out.bdg mm10.chrom.sizes GAS_fimo_out.bw
bedGraphToBigWig ISRE_fimo_out.bdg mm10.chrom.sizes ISRE_fimo_out.bw
bedGraphToBigWig CTCF_fimo_out.bdg mm10.chrom.sizes CTCF_fimo_out.bw

echo $(date +"[%b %d %H:%M:%S] Done")