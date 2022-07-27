All bash scripts were written to be run using SLURM on the HPC cluster at the University of Colorado Boulder. Some scripts are written to be run in parallel as job arrays.

For simplicity all input and output files are written to the same directory.

Workflow:
1. [fimo.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/motif_analysis/fimo.sbatch)
2. [convert_fimo_txt_to_bw.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/motif_analysis/convert_fimo_txt_to_bw.sbatch)
3. [run_createExpandedRepeatFile.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/motif_analysis/run_createExpandedRepeatFile.sbatch)
4. [B2_Mm2_motif_heatmap.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/motif_analysis/B2_Mm2_motif_heatmap.sbatch), [B2_Mm1a_motif_heatmap.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/motif_analysis/B2_Mm1a_motif_heatmap.sbatch), [B2_Mm1t_motif_heatmap.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/motif_analysis/B2_Mm1t_motif_heatmap.sbatch)

Required packages:
* FIMO v5.4.1 (https://meme-suite.org/meme/)
* bedtools v2.28.0 (http://bedtools.readthedocs.io/en/latest/)
* Samtools v1.10 (http://www.htslib.org/)
* bedGraphToBigWig v4 (http://hgdownload.soe.ucsc.edu/downloads.html#source_downloads)
* deepTools v3.5.1 (https://deeptools.readthedocs.io/en/develop/index.html)

Note that the soft masked mm10 assembly file may be downloaded via:
* mm10.fa - http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz

The RepeatMasker output file for mm10 using Dfam 2.0 annotation may be downloaded via:
* mm10.fa.out (referenced in scripts as mm10_dfam.out) - https://www.repeatmasker.org/genomes/mm10/RepeatMasker-rm406-dfam2.0/mm10.fa.out.gz

Other files referenced in these scripts that have been provided through this repository:
* [MA0137.3.meme](https://github.com/coke6162/B2_SINE_enhancers/blob/main/motif_analysis/MA0137.3.meme), [MA0139.1.meme](https://github.com/coke6162/B2_SINE_enhancers/blob/main/motif_analysis/MA0139.1.meme), [MA0517.1.meme](https://github.com/coke6162/B2_SINE_enhancers/blob/main/motif_analysis/MA0517.1.meme)
* [B2_Mm2_expanded.bed.gz](https://github.com/coke6162/B2_SINE_enhancers/blob/main/motif_analysis/B2_Mm2_expanded.bed.gz), [B2_Mm1a_expanded.bed.gz](https://github.com/coke6162/B2_SINE_enhancers/blob/main/motif_analysis/B2_Mm1a_expanded.bed.gz), [B2_Mm1t_expanded.bed.gz](https://github.com/coke6162/B2_SINE_enhancers/blob/main/motif_analysis/B2_Mm1t_expanded.bed.gz)
* [repeats_turnRMFile_to_UCSC.py](https://github.com/coke6162/B2_SINE_enhancers/blob/main/motif_analysis/repeats_turnRMFile_to_UCSC.py)
* [createExpandedRepeatFile.py](https://github.com/coke6162/B2_SINE_enhancers/blob/main/motif_analysis/createExpandedRepeatFile.py)