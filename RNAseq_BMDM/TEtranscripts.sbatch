#!/bin/bash

# Script for running TEtranscripts

# Example usage:
# sbatch TEtranscripts.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=24-0:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J TEtranscripts
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Load TEtranscripts
module load python/2.7.14/tetranscripts/2.1.4

# Specify R library paths
export R_LIBS='/Users/path/to/R/x86_64-redhat-linux-gnu-library/3.6/'

# Run TEtranscripts
pwd; hostname; date

echo "TEtranscripts version: "$(TEtranscripts --version)

# (1) piccolo IFNG 2h vs UT
echo $(date +"[%b %d %H:%M:%S] Running TEtranscripts (Piccolo IFNG 2h & UT)...")

TEtranscripts \
--sortByPos --mode multi \
--iteration 100 --stranded no \
--GTF gencode.vM18.annotation.gtf \
--TE mm10_dfam_B2_Mm2_picc_IFNG_2h_STAT1_modified.gtf \
--project piccolo_IFNG_2h_vs_UT \
--treatment picc_BMDM_IFNG_2h_R1_k100.sorted.bam picc_BMDM_IFNG_2h_R2_k100.sorted.bam picc_BMDM_IFNG_2h_R3_k100.sorted.bam \
--control picc_BMDM_UT_R1_k100.sorted.bam picc_BMDM_UT_R2_k100.sorted.bam picc_BMDM_UT_R3_k100.sorted.bam

# (2) piccolo IFNG 4h vs UT
echo $(date +"[%b %d %H:%M:%S] Running TEtranscripts (Piccolo IFNG 4h & UT)...")

TEtranscripts \
--sortByPos --mode multi \
--iteration 100 --stranded no \
--GTF gencode.vM18.annotation.gtf \
--TE mm10_dfam_B2_Mm2_picc_IFNG_2h_STAT1_modified.gtf \
--project piccolo_IFNG_4h_vs_UT \
--treatment picc_BMDM_IFNG_4h_R1_k100.sorted.bam picc_BMDM_IFNG_4h_R2_k100.sorted.bam picc_BMDM_IFNG_4h_R3_k100.sorted.bam \
--control picc_BMDM_UT_R1_k100.sorted.bam picc_BMDM_UT_R2_k100.sorted.bam picc_BMDM_UT_R3_k100.sorted.bam

# (3) platanitis IFNG 2h vs UT
echo $(date +"[%b %d %H:%M:%S] Running TEtranscripts (Platanitis IFNG 2h & UT)...")

TEtranscripts \
--sortByPos --mode multi \
--iteration 100 --stranded no \
--GTF gencode.vM18.annotation.gtf \
--TE mm10_dfam_B2_Mm2_plat_IFNG_1.5h_STAT1_modified.gtf \
--project platanitis_IFNG_2h_vs_UT \
--treatment plat_BMDM_IFNG_2h_R1_k100.sorted.uniq.bam plat_BMDM_IFNG_2h_R2_k100.sorted.uniq.bam plat_BMDM_IFNG_2h_R3_k100.sorted.uniq.bam \
--control plat_BMDM_UT_R1_k100.sorted.uniq.bam plat_BMDM_UT_R2_k100.sorted.uniq.bam plat_BMDM_UT_R3_k100.sorted.uniq.bam

# Filter out genes and remove class info
echo $(date +"[%b %d %H:%M:%S] Filtering out genes and removing TE class information...")

grep -v "ENSMUS" piccolo_IFNG_2h_vs_UT.cntTable | tr ':' '\t' | cut -f1,4- > piccolo_IFNG_2h_vs_UT_TEsonly.cntTable
grep -v "ENSMUS" piccolo_IFNG_4h_vs_UT.cntTable | tr ':' '\t' | cut -f1,4- > piccolo_IFNG_4h_vs_UT_TEsonly.cntTable
grep -v "ENSMUS" platanitis_IFNG_2h_vs_UT.cntTable | tr ':' '\t' | cut -f1,4- > platanitis_IFNG_2h_vs_UT_TEsonly.cntTable

echo $(date +"[%b %d %H:%M:%S] Done")