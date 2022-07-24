#!/bin/bash

# Script for identifying enriched motifs using XSTREME v5.4.1
# Motif db file available through MEME Suite v5.4.1 isntallation

# Example usage:
# sbatch xstreme.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J xstreme
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Load bedtools
module load bedtools/2.28.0

# Convert bed to fasta
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Converting bed to fasta...")

bedtools getfasta -fi mm10.fa -bed picc_BMDM_WT_IFNG_2h_STAT1_intersected.bed -fo picc_BMDM_WT_IFNG_2h_STAT1_intersected.fa
bedtools getfasta -fi mm10.fa -bed plat_BMDM_WT_IFNG_1.5h_STAT1_intersected.bed -fo plat_BMDM_WT_IFNG_1.5h_STAT1_intersected.fa

# Run XSTREME
echo $(date +"[%b %d %H:%M:%S] Running XSTREME...")

for i in picc_BMDM_WT_IFNG_2h_STAT1_intersected.fa plat_BMDM_WT_IFNG_1.5h_STAT1_intersected.fa
do
	echo "Sample: "${i}
	xstreme --verbosity 2 --oc ${i%.fa} --dna --p ${i} --m JASPAR2018_CORE_vertebrates_non-redundant.meme --minw 6 --maxw 20 --streme-nmotifs 20 --align center
done

echo $(date +"[%b %d %H:%M:%S] Done")