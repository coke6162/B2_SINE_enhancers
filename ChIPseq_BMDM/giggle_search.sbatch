#!/bin/bash

# Script to search genomic intervals using GIGGLE v0.6.3

# sbatch giggle_search.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time=1:00:00
#SBATCH --mem=100MB
#SBATCH --propagate=NONE

# Job name and output
#SBATCH -J giggle
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# load modules
module load samtools/1.14

# Make a temporary directory
pwd; hostname; date

echo "Making a temporary directory..."

mkdir -p giggle_tmp/

# Copy files to temporary directory and add ".bed" suffix (if necessary; required for GIGGLE)
echo $(date +"[%b %d %H:%M:%S] Adding '.bed' suffix...")

cp picc_BMDM_WT_UT_H3K27ac_peaks.narrowPeak giggle_tmp/picc_BMDM_WT_UT_H3K27ac_peaks.narrowPeak.bed
cp picc_BMDM_WT_IFNG_2h_H3K27ac_peaks.narrowPeak giggle_tmp/picc_BMDM_WT_IFNG_2h_H3K27ac_peaks.narrowPeak.bed
cp picc_BMDM_WT_UT_STAT1_peaks.narrowPeak giggle_tmp/picc_BMDM_WT_UT_STAT1_peaks.narrowPeak.bed
cp picc_BMDM_WT_IFNG_2h_STAT1_intersected.bed giggle_tmp/picc_BMDM_WT_IFNG_2h_STAT1_intersected.bed
cp plat_BMDM_WT_UT_STAT1_intersected.bed giggle_tmp/plat_BMDM_WT_UT_STAT1_intersected.bed
cp plat_BMDM_WT_IFNG_1.5h_STAT1_intersected.bed giggle_tmp/plat_BMDM_WT_IFNG_1.5h_STAT1_intersected.bed

# Compress each file
echo $(date +"[%b %d %H:%M:%S] Compressing files ending in '.bed'...")

bgzip giggle_tmp/*.bed

# Run GIGGLE search
echo $(date +"[%b %d %H:%M:%S] Starting GIGGLE search...")

for i in giggle_tmp/*.bed.gz
do
	giggle search -i giggle_index/indexed -q ${i} -s -g 2730871774 > ${i%.bed.gz}_giggle.txt.tmp
done

# Remove TE paths from output and sort by descending GIGGLE score
echo $(date +"[%b %d %H:%M:%S] Sorting by descending combo score...")

for i in giggle_tmp/*_giggle.txt.tmp
do
	awk '(NR>1){print $0}' ${i} | sed 's/.bed.gz//g' | sed "s|giggle_index/sorted/||g" | sort -k8,8rn ${i} > ${i%.tmp}
	cp ${i%.tmp} .
done

# Remove temporary directory
echo $(date +"[%b %d %H:%M:%S] Removing temporary directory...")

rm -r giggle_tmp/

echo $(date +"[%b %d %H:%M:%S] Done!")