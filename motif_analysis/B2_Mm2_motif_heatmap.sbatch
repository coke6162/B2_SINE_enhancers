#!/bin/bash

# Script for preparing a heatmap showing predicted motif signal over B2_Mm2
# B2_Mm2 coordinates have been adjusted based on position relative to consensus according to RepeatMasker output 
# See run_createExpandedRepeatFile.sbatch

# Example usage:
# sbatch B2_Mm2_heatmap.sbatch

# General settings
#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=2-0:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J heatmap
#SBATCH -o /Users/%u/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/slurmErr/slurm-%j.err

# load modules
module purge
module load python/3.6.3

# Define temporary directory
tmpDir=tmp/${SLURM_ARRAY_JOB_ID}

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
--referencePoint TSS \
--scoreFileName GAS_fimo_out.bw CTCF_fimo_out.bw \
--regionsFileName B2_Mm2_expanded.bed \
--beforeRegionStartLength 100 --afterRegionStartLength 500 --binSize 1 \
--missingDataAsZero \
-o B2_Mm2_motif_heatmap.mat.gz -p 8

# Use "plotHeatmap" to create a pdf
echo $(date +"[%b %d %H:%M:%S] Starting deeptools plotHeatmap...")

plotHeatmap \
-m B2_Mm2_motif_heatmap.mat.gz \
--outFileName B2_Mm2_motif_heatmap.pdf \
--outFileSortedRegions B2_Mm2_motif_heatmap.dt.bed \
--outFileNameMatrix B2_Mm2_motif_heatmap.matrix.tab \
--sortRegions descend --sortUsing mean --sortUsingSamples 1 2 \
--colorMap Blues Greens \
--zMin 0 --zMax "8 8" \
--yMin 0 --yMax "8 8" \
--regionsLabel "B2_Mm2" \
--samplesLabel "GAS CTCF" \
--plotTitle B2_Mm2_motif_heatmap \
--heatmapWidth 10

# Remove intermediate files
echo $(date +"[%b %d %H:%M:%S] Removing intermediate files...")

rm -rf ${tmpDir}

echo $(date +"[%b %d %H:%M:%S] Done")