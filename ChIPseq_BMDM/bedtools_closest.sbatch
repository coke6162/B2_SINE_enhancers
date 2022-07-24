#!/bin/bash

# Script for getting absolute distances from STAT1-bound/unbound B2 SINEs to nearest ISG, IRG, or nonresponsive gene
# Note that only the top 750 ISGs and IRGs & random 750 nonresponsive genes are included for this analysis
# To generate these gene lists (including Piccolo et al. 2017 IFNG 2h vs UT & Platanitis et al. 2019 IFNG 2h vs UT), see "RNAseq_BMDM" repos

# Example usage:
# sbatch get_overlapping_B2_nearest_neighbor.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J B2_Mm2
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# Load bedtools
module load bedtools/2.28.0

# Make subdirectory
mkdir -p nearest_neighbor/piccolo_RNAseq_IFNG_4h_vs_UT
mkdir -p nearest_neighbor/piccolo_RNAseq_IFNG_2h_vs_UT
mkdir -p nearest_neighbor/platanitis_RNAseq_IFNG_2h_vs_UT

# Navigate to "nearest_neighbor" subdirectory
cd nearest_neighbor/

## Compare distance from B2 (defined using Piccolo et al. 2017 STAT1 ChIP-seq) to nearest Piccolo et al. 2017 IFNG 4h vs UT gene (Fig 1B)

# Get distance from B2 (defined using Piccolo et al. 2017 STAT1 ChIP-seq) to nearest Piccolo et al. 2017 IFNG 4h vs UT ISG
for i in B2*picc_BMDM_WT_IFNG_2h_STAT1_intersected*bed
do
	bedtools closest -d -t first -a ${i} -b picc_BMDM_IFNG_4h_vs_UT_ISGs_top750.bed | sort -k10,10rg > piccolo_RNAseq_IFNG_4h_vs_UT/${i%.bed}_nearest_4h_vs_UT_ISG_top750.bed
done

# Get distance from B2 (defined using Piccolo et al. 2017 STAT1 ChIP-seq) to nearest Piccolo et al. 2017 IFNG 4h vs UT IRG
for i in B2*picc_BMDM_WT_IFNG_2h_STAT1_intersected*bed
do
	bedtools closest -d -t first -a ${i} -b picc_BMDM_IFNG_4h_vs_UT_IRGs_top750.bed | sort -k10,10rg > piccolo_RNAseq_IFNG_4h_vs_UT/${i%.bed}_nearest_4h_vs_UT_IRG_top750.bed
done

# Get distance from B2 (defined using Piccolo et al. 2017 STAT1 ChIP-seq) to nearest Piccolo et al. 2017 IFNG 4h vs UT nonresponsive gene
for i in B2*picc_BMDM_WT_IFNG_2h_STAT1_intersected*bed
do
	bedtools closest -d -t first -a ${i} -b picc_BMDM_IFNG_4h_vs_UT_nonresponsive_padj0.90_baseMean100_log2FC0.10_random750.bed | sort -k10,10rg > piccolo_RNAseq_IFNG_4h_vs_UT/${i%.bed}_nearest_4h_vs_UT_nonresponsive_random750.bed
done

## Compare distance from B2 (defined using Piccolo et al. 2017 STAT1 ChIP-seq) to nearest Piccolo et al. 2017 IFNG 2h vs UT gene (Fig S1A)

# Get distance from B2 (defined using Piccolo et al. 2017 STAT1 ChIP-seq) to nearest Piccolo et al. 2017 IFNG 2h vs UT ISG
for i in B2*picc_BMDM_WT_IFNG_2h_STAT1_intersected*bed
do
	bedtools closest -d -t first -a ${i} -b picc_BMDM_IFNG_2h_vs_UT_ISGs_top750.bed | sort -k10,10rg > piccolo_RNAseq_IFNG_2h_vs_UT/${i%.bed}_nearest_2h_vs_UT_ISG_top750.bed
done

# Get distance from B2 (defined using Piccolo et al. 2017 STAT1 ChIP-seq) to nearest Piccolo et al. 2017 IFNG 2h vs UT IRG
for i in B2*picc_BMDM_WT_IFNG_2h_STAT1_intersected*bed
do
	bedtools closest -d -t first -a ${i} -b picc_BMDM_IFNG_2h_vs_UT_IRGs_top750.bed | sort -k10,10rg > piccolo_RNAseq_IFNG_2h_vs_UT/${i%.bed}_nearest_2h_vs_UT_IRG_top750.bed
done

# Get distance from B2 (defined using Piccolo et al. 2017 STAT1 ChIP-seq) to nearest Piccolo et al. 2017 IFNG 2h vs UT nonresponsive gene
for i in B2*picc_BMDM_WT_IFNG_2h_STAT1_intersected*bed
do
	bedtools closest -d -t first -a ${i} -b picc_BMDM_IFNG_2h_vs_UT_nonresponsive_padj0.90_baseMean100_log2FC0.10_random750.bed | sort -k10,10rg > piccolo_RNAseq_IFNG_2h_vs_UT/${i%.bed}_nearest_2h_vs_UT_nonresponsive_random750.bed
done

## Compare distance from B2 (defined using Platanitis et al. 2019 STAT1 ChIP-seq) to nearest Platanitis et al. 2019 IFNG 2h vs UT gene (Fig S1B)

# Get distance from B2 (defined using Platanitis et al. 2019 STAT1 ChIP-seq) to nearest Platanitis et al. 2019 IFNG 2h vs UT ISG
for i in B2*plat_BMDM_WT_IFNG_1.5h_STAT1_intersected*bed
do
	bedtools closest -d -t first -a ${i} -b plat_BMDM_IFNG_2h_vs_UT_ISGs_top750.bed | sort -k10,10rg > platanitis_RNAseq_IFNG_2h_vs_UT/${i%.bed}_nearest_2h_vs_UT_ISG_top750.bed
done

# Get distance from B2 (defined using Platanitis et al. 2019 STAT1 ChIP-seq) to nearest Platanitis et al. 2019 IFNG 2h vs UT IRG
for i in B2*plat_BMDM_WT_IFNG_1.5h_STAT1_intersected*bed
do
	bedtools closest -d -t first -a ${i} -b plat_BMDM_IFNG_2h_vs_UT_IRGs_top750.bed | sort -k10,10rg > platanitis_RNAseq_IFNG_2h_vs_UT/${i%.bed}_nearest_2h_vs_UT_IRG_top750.bed
done

# Get distance from B2 (defined using Platanitis et al. 2019 STAT1 ChIP-seq) to nearest Platanitis et al. 2019 IFNG 2h vs UT nonresponsive gene
for i in B2*plat_BMDM_WT_IFNG_1.5h_STAT1_intersected*bed
do
	bedtools closest -d -t first -a ${i} -b plat_BMDM_IFNG_2h_vs_UT_nonresponsive_padj0.90_baseMean100_log2FC0.10_random750.bed | sort -k10,10rg > platanitis_RNAseq_IFNG_2h_vs_UT/${i%.bed}_nearest_2h_vs_UT_nonresponsive_random750.bed
done