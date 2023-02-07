All bash scripts were written to be run using SLURM on the HPC cluster at the University of Colorado Boulder. Some scripts are written to be run in parallel as job arrays.

Accession numbers and corresponding sample names are provided in [sample_names.txt](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/sample_names.txt). These scripts assume that all samples are named as described in [sample_names.txt](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/sample_names.txt). Samples may be downloaded from SRA using the provided SRR IDs using sra-tools (script not provided). For simplicity all input and output files are written to the same directory.

The purpose of this analysis is to determine whether B2_Mm2 shows evidence of broadly shaping mouse innate immunity by defining mouse-specific ISGs with a predicted B2_Mm2-derived enhancer. 

First, we re-analyzed publicly available RNA-seq data from human, IFNG-stimulated CD14+ monocytes to identify human ISGs. Note that, although classified as paired-end, these data when downloaded from SRA result in empty "R2" files (which can be confirmed by navigating to the samples on SRA). For this analysis these were deleted, and "R1" files were used as single-end data. Workflow:
1. [bbduk.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/bbduk.sbatch)
2. [fastqc.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/fastqc.sbatch)
3. [multiqc.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/multiqc.sbatch)
4. [hisat2.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/hisat2.sbatch)
5. [featureCounts.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/featureCounts.sbatch)
6. [DESeq2.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/DESeq2.R)

All human-to-mouse one-to-one orthologs identified with a confidence score of 1 were downloaded from BioMart and used to convert human ISG gene IDs to mouse gene IDs. Mouse and human-to-mouse gene IDs were collapsed into a single list and binned according to species specificity (i.e. induced in mouse, induced in human, induced in both). Note that any instance of "human" or "human-specific" is relative to mouse - human ISGs for which there were no one-to-one orthologs in mouse are not considered in this analysis. This assumes that you have a list of mouse ISGs, which here we defined using the Piccolo et al. 2017 IFNG 4h vs UT comparison generated from the [RNAseq_BMDM repos](https://github.com/coke6162/B2_SINE_enhancers/tree/main/RNAseq_BMDM). We then sought to identify (1) the nearest STAT1-bound B2_Mm2 element relative to each ISG TSS and (2) the top putative B2_Mm2 enhancer predicted to interact with an ISG. Workflow:
1. [generate_binary_matrix.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/generate_binary_matrix.sbatch)
2. [identify_nearest_STAT1_B2_Mm2.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/identify_nearest_STAT1_B2_Mm2.sbatch)
3. [identify_interacting_ABC_B2_Mm2.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/identify_interacting_ABC_B2_Mm2.sbatch)

Required packages:
* BBMap v38.05 (https://jgi.doe.gov/data-and-tools/bbtools/)
* FastQC v0.11.8 (https://github.com/s-andrews/FastQC)
* MultiQC v1.7 (https://github.com/ewels/MultiQC)
* HISAT2 v2.1.0 (https://github.com/DaehwanKimLab/hisat2)
* Samtools v1.14 (http://www.htslib.org/)
* Subread v1.6.2 (http://subread.sourceforge.net/)
* DESeq2 v1.26.0 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
* bedtools v2.28.0 (http://bedtools.readthedocs.io/en/latest/)

Note that the soft masked hg38 assembly and Gencode v39 gene annotation file may be downloaded via:
* [hg38.fa](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)
* [gencode.v39.annotation.gtf](https://www.gencodegenes.org/human/release_39.html) 

Other files referenced in these scripts that have been provided through this repository:
* [adapters.fa](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/adapters.fa)
* [bam_order.txt](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/bam_order.txt)
* [human_ensembl_v105_biomart_to_mouse_stable.txt](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/human_ensembl_v105_biomart_to_mouse_stable.txt)
* [B2_Mm2_picc_IFNG_2h_STAT1.bed](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/B2_Mm2_picc_IFNG_2h_STAT1.bed) - note that this includes B2_Mm2 bound by CTCF and is different from [the provided file to generate deepTools heatmaps](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/B2_Mm2_picc_BMDM_WT_IFNG_2h_STAT1_intersected_only.bed)
* [B2_Mm2.bed](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/B2_Mm2.bed)