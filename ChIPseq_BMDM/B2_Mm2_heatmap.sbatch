#!/bin/bash

# Script for preparing a heatmap showing ChIP-seq signal over B2_Mm2
# These files are meant to be used with "B2_Mm2_heatmap.sbatch"

# Example usage:
# sbatch B2_Mm2_heatmap.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=24:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J heatmap
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# load modules
module purge
module load python/3.6.3

# Define temporary directory
tmpDir=ChIP_heatmap/tmp/${SLURM_ARRAY_JOB_ID}

# Make temporary directories
pwd; hostname; date

echo "Making temporary directory..."

mkdir -p ${tmpDir}

# Define temporary directory for deeptools
# Note that if you do not set this deeptools will write temporary files to /var and your jobs will likely be killed
TMPDIR="${tmpDir}"
export TMPDIR

# Use computeMatrix to generate data underlying heatmap
echo $(date +"[%b %d %H:%M:%S] Starting deeptools computeMatrix...")

computeMatrix reference-point \
--referencePoint center \
--scoreFileName picc_BMDM_WT_UT_H3K27ac.bw picc_BMDM_WT_IFNG_2h_H3K27ac.bw picc_BMDM_WT_UT_STAT1.bw picc_BMDM_WT_IFNG_2h_STAT1_R1.bw plat_BMDM_WT_UT_STAT1_R1.bw plat_BMDM_WT_IFNG_1.5h_STAT1_R1.bw gual_BMDM_shNT_UT_CTCF_R1.bw cuar_BMDM_WT_UT_RAD21_R1.bw GAS_fimo_out.bw CTCF_fimo_out.bw \
--regionsFileName ChIP_heatmap/B2_Mm2_picc_BMDM_WT_IFNG_2h_STAT1_intersected_only.bed ChIP_heatmap/B2_Mm2_picc_BMDM_WT_IFNG_2h_STAT1_intersected_and_CTCF.bed ChIP_heatmap/B2_Mm2_gual_BMDM_shNT_UT_CTCF_intersected_only.bed ChIP_heatmap/B2_Mm2_control_no_STAT1_or_CTCF_shuf1000.bed \
--beforeRegionStartLength 4000 --afterRegionStartLength 4000 \
--binSize 10 --missingDataAsZero \
-o ChIP_heatmap/B2_Mm2_heatmap.mat.gz -p 8

# Use "plotHeatmap" to create a pdf
echo $(date +"[%b %d %H:%M:%S] Starting deeptools plotHeatmap...")

plotHeatmap \
-m ChIP_heatmap/B2_Mm2_heatmap.mat.gz \
--outFileName ChIP_heatmap/B2_Mm2_heatmap.pdf \
--outFileSortedRegions ChIP_heatmap/B2_Mm2_heatmap.dt.bed \
--outFileNameMatrix ChIP_heatmap/B2_Mm2_heatmap.matrix.tab \
--sortRegions descend --sortUsingSamples 2 4 6 --colorMap Blues \
--zMin 0 --zMax "1 1 1 1 1 1 1 1 1 1" \
--yMin 0 --yMax "2 2 2 2 2 2 2 2 6 12" \
--regionsLabel "B2_Mm2_picc_BMDM_WT_IFNG_2h_STAT1_intersected_only B2_Mm2_picc_BMDM_WT_IFNG_2h_STAT1_intersected_and_CTCF B2_Mm2_gual_BMDM_shNT_UT_CTCF_intersected_only B2_Mm2_control_no_STAT1_or_CTCF_shuf1000" \
--samplesLabel "picc_BMDM_WT_UT_H3K27ac picc_BMDM_WT_IFNG_2h_H3K27ac picc_BMDM_WT_UT_STAT1 picc_BMDM_WT_IFNG_2h_STAT1_R1 plat_BMDM_WT_UT_STAT1_R1 plat_BMDM_WT_IFNG_1.5h_STAT1_R1 gual_BMDM_shNT_UT_CTCF_R1 cuar_BMDM_WT_UT_RAD21_R1 GAS CTCF" \
--plotTitle B2_Mm2_heatmap \
--heatmapWidth 10

# Remove intermediate files
echo $(date +"[%b %d %H:%M:%S] Removing intermediate files...")

rm -rf ${tmpDir}

echo $(date +"[%b %d %H:%M:%S] Done")