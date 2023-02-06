All bash scripts were written to be run using SLURM on the HPC cluster at the University of Colorado Boulder. Some scripts are written to be run in parallel as job arrays.

Accession numbers and corresponding sample names are provided in [sample_names.txt](). These scripts assume that all samples are named as described in [sample_names.txt](). Samples may be download from SRA using the provided SRR IDs using sra-tools (script not provided). For simplicity all input and output files are written to the same directory.

The purpose of this analysis is to generate alignment and peak files using ATAC-seq data from IFNG-stimulated BMDMs [Platanitis et al. iScience 2022](https://doi.org/10.1016/j.isci.2022.103840). These files will be used in conjunction with H3K27ac ChIP, Hi-C, and RNA-seq data from IFNG-stimulated BMDMs to predict enhancer-gene contacts using the Activity by Contact model. 

A typical RNA-seq workflow to call differentially expressed genes looks like this:
1. [bbduk.sbatch]()
2. [fastqc.sbatch]()
3. [multiqc.sbatch]()
4. [bowtie2.sbatch]()
5. [remove_duplicates.sbatch]()
6. [shift_fragments.sbatch]()
7. [macs2.sbatch]()

Required packages:
* BBMap v38.05 (https://jgi.doe.gov/data-and-tools/bbtools/)
* FastQC v0.11.8 (https://github.com/s-andrews/FastQC)
* MultiQC v1.7 (https://github.com/ewels/MultiQC)
* Bowtie2 v2.2.9 (https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* Samtools v1.14 (http://www.htslib.org/)
* Picard v2.6.0 (https://broadinstitute.github.io/picard/)
* deepTools v3.5.1 (https://deeptools.readthedocs.io/en/develop/index.html)
* MACS2 v2.1.1 (https://pypi.org/project/MACS2/)

Note that the soft masked mm10 assembly and Gencode vM18 gene annotation file may be downloaded via:
* mm10.fa - http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz

Other files referenced in these scripts that have been provided through this repository:
* [adapters.fa]()
* [bam_order.txt]()