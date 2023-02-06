All bash scripts were written to be run using SLURM on the HPC cluster at the University of Colorado Boulder. Some scripts are written to be run in parallel as job arrays.

The purpose of this analysis is to predict enhancer-gene contacts in IFNG-stimulated BMDMs using the Activity-by-Contact (ABC) Model, which incorporates both epigenomic signal and Hi-C 3D interaction data. As input into the ABC model, we used publicly available ATAC-seq and Hi-C data (PRJNA694816) from murine BMDMs stimulated with IFNG along with H3K27ac ChIP-seq (SRR3929173) data. Put explicity, this analysis requires the following files:
* ATAC-seq alignment & peak files (plat_BMDM_IFNG_2h_ATAC_R1.sorted.dedup.bam, plat_BMDM_IFNG_2h_ATAC_R2.sorted.dedup.bam, plat_BMDM_IFNG_2h_ATAC_R3.sorted.dedup.bam, plat_BMDM_IFNG_2h_ATAC_R1_peaks.narrowPeak)
* H3K27ac alignment file (picc_BMDM_WT_IFNG_2h_H3K27ac.sorted.bam)
* Directory containing KR-normalized Hi-C data (plat_BMDM_IFNG_2h_HiC) parsed by chromosome, one per subdirectory (from [juicebox_dump.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/HiC_BMDM/juicebox_dump.sbatch))

Workflow:
1. [call_candidate_regions.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ABC_BMDM/call_candidate_regions.sbatch)
2. [collapse_gene_boundaries.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ABC_BMDM/collapse_gene_boundaries.sbatch)
3. [find_neighborhoods.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ABC_BMDM/find_neighborhoods.sbatch)
4. [predict_enhancers.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ABC_BMDM/predict_enhancers.sbatch)
5. [subset_and_intersect_enhancers.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ABC_BMDM/subset_and_intersect_enhancers.sbatch)

Required packages:
* bedtools v2.28.0 (http://bedtools.readthedocs.io/en/latest/)
* Activity-by-Contact Model v0.2.2 (https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction)

Note that the soft masked mm10 assembly and Gencode vM18 gene annotation file may be downloaded via:
* [mm10-blocklist.v2.bed](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz)

Other files referenced in these scripts that have been provided through this repository:
* [mm10.chrom.sizes](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ABC_BMDM/mm10.chrom.sizes)
* [gencode.vM18.annotation.bed](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ABC_BMDM/gencode.vM18.annotation.bed)
* [gencode.vM18.annotation.boundaries.uniq.bed](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ABC_BMDM/gencode.vM18.annotation.boundaries.uniq.bed)
* [mm10_chrom_lengths.bed](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ABC_BMDM/mm10_chrom_lengths.bed)