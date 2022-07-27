#!/bin/bash

# Script for getting top unique GAS motif (by p-value) overlapping B2 elements
# This script assumes you have genome-wide motif tracks where FIMO was run with a p-value threshold of 1 (see fimo_pval_1.sbatch)

# Example usage:
# sbatch get_fimo_pval_box_whisker.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J fimo_pval
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Load bedtools
module load bedtools/2.28.0

# Identify top unique GAS motif (by p-value) overlapping B2 elements
for i in box_and_whisker/B2_Mm2_picc_BMDM_WT_IFNG_2h_STAT1_intersected.bed box_and_whisker/B2_Mm2_control_no_STAT1_or_CTCF_shuf5000.bed box_and_whisker/B2_Mm1a_control_no_STAT1_or_CTCF_shuf5000.bed box_and_whisker/B2_Mm1t_control_no_STAT1_or_CTCF_shuf5000.bed box_and_whisker/B2_Mm2_gual_BMDM_shNT_UT_CTCF_intersected_no_STAT1.bed
do
	bedtools intersect -a ${i} -b box_and_whisker/GAS_FIMO_out_pval_1.bed -wao | sort -k4,4 -k11,11g - | sort -k4,4 -u | bedtools sort -i - > ${i%.bed}_GAS.bed
done

# Generate summary file for STAT1-bound B2_Mm2
awk -v '{if ($13 > 0) print "B2_Mm2_STAT1" "\t" $11; else print "B2_Mm2_STAT1" "\t" 1}' box_and_whisker/B2_Mm2_picc_BMDM_WT_IFNG_2h_STAT1_intersected_GAS.bed > box_and_whisker/B2_Mm2_picc_BMDM_WT_IFNG_2h_STAT1_intersected_GAS_summarized.txt

# Generate summary file for unbound B2_Mm2
awk -v '{if ($13 > 0) print "B2_Mm2_no_STAT1_CTCF" "\t" $11; else print "B2_Mm2_no_STAT1_CTCF" "\t" 1}' box_and_whisker/B2_Mm2_control_no_STAT1_or_CTCF_shuf5000_GAS.bed > box_and_whisker/B2_Mm2_control_no_STAT1_or_CTCF_shuf5000_GAS_summarized.txt

# Generate summary file for unbound B2_Mm1a
awk -v '{if ($13 > 0) print "B2_Mm1a_no_STAT1_CTCF" "\t" $11; else print "B2_Mm1a_no_STAT1_CTCF" "\t" 1}' box_and_whisker/B2_Mm1a_control_no_STAT1_or_CTCF_shuf5000_GAS.bed > box_and_whisker/B2_Mm1a_control_no_STAT1_or_CTCF_GAS_summarized.txt

# Generate summary file for unbound B2_Mm1t
awk -v '{if ($13 > 0) print "B2_Mm1t_no_STAT1_CTCF" "\t" $11; else print "B2_Mm1t_no_STAT1_CTCF" "\t" 1}' box_and_whisker/B2_Mm1t_control_no_STAT1_or_CTCF_shuf5000_GAS.bed > box_and_whisker/B2_Mm1t_control_no_STAT1_or_CTCF_shuf5000_GAS_summarized.txt

# Generate summary file for B2_Mm2 bound by CTCF but not STAT1
awk -v '{if ($13 > 0) print "B2_Mm2_CTCF_no_STAT1" "\t" $11; else print "B2_Mm2_CTCF_no_STAT1" "\t" 1}' box_and_whisker/B2_Mm2_gual_BMDM_shNT_UT_CTCF_intersected_no_STAT1_GAS.bed > box_and_whisker/B2_Mm2_picc_BMDM_WT_IFNG_2h_STAT1_intersected_GAS_summarized.txt

# Concatenate summary files
cat box_and_whisker/*_summarized.txt > box_and_whisker/B2_GAS_pval_all_summarized.txt