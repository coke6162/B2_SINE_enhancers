All bash scripts were written to be run using SLURM on the HPC cluster at the University of Colorado Boulder. Some scripts are written to be run in parallel as job arrays.

Accession numbers and corresponding sample names are provided in [sample_names.txt](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ATACseq_BMDM/sample_names.txt). These scripts assume that all samples are named as described in [sample_names.txt](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ATACseq_BMDM/sample_names.txt). Samples may be downloaded from SRA using the provided SRR IDs using sra-tools (script not provided). For simplicity all input and output files are written to the same directory.

The purpose of this analysis is to generate a Hi-C matrix from IFNG-stimulated BMDMs ([Platanitis et al. iScience 2022](https://doi.org/10.1016/j.isci.2022.103840)) to predict enhancer-gene contacts using the Activity by Contact model. Note that adapter trimming was skipped as we observed no evidence of adapter contamination and low quality sequences after downloading from SRA. 

Workflow:
1. [bwa.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/HiC_BMDM/bwa.sbatch)
2. [pairtools_parse.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/HiC_BMDM/pairtools_parse.sbatch)
3. [pairtools_sort.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/HiC_BMDM/pairtools_sort.sbatch)
4. [pairtools_merge.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/HiC_BMDM/pairtools_merge.sbatch)
5. [pairtools_dedup.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/HiC_BMDM/pairtools_dedup.sbatch)
6. [pairtools_filter.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/HiC_BMDM/pairtools_filter.sbatch)
7. [run_generate_site_positions.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/HiC_BMDM/run_generate_site_positions.sbatch)
8. [run_addfrag2pairs.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/HiC_BMDM/run_addfrag2pairs.sbatch)
9. [juicer_pre.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/HiC_BMDM/juicer_pre.sbatch)
10. [juicebox_dump.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/HiC_BMDM/juicebox_dump.sbatch)

Required packages:
* BBMap v38.05 (https://jgi.doe.gov/data-and-tools/bbtools/)
* FastQC v0.11.8 (https://github.com/s-andrews/FastQC)
* MultiQC v1.7 (https://github.com/ewels/MultiQC)
* BWA v0.7.17 (https://github.com/lh3/bwa)
* Samtools v1.14 (http://www.htslib.org/)
* Pairtools v1.0.2 (https://pairtools.readthedocs.io/en/latest/)
* Pairix v0.3.7 (https://github.com/4dn-dcic/pairix)
* Juicer v1.6 (https://github.com/aidenlab/juicer)
* bedtools v2.28.0 (http://bedtools.readthedocs.io/en/latest/)

Note that the soft masked mm10 assembly and Gencode vM18 gene annotation file may be downloaded via:
* [mm10.fa](http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz)

Other files referenced in these scripts that have been provided through this repository:
* [mm10.chrom.sizes](https://github.com/coke6162/B2_SINE_enhancers/blob/main/HiC_BMDM/mm10.chrom.sizes)
* [sample_names.txt](https://github.com/coke6162/B2_SINE_enhancers/blob/main/HiC_BMDM/sample_names.txt)